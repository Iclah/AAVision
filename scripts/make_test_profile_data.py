"""
Build the bundled test dataset used by `-profile test` (test_data/).

The full host genome (~3 GB) is too large for the repository, so the host
reference is cut down to the regions the test reads can reach:

  1. Align ALL reads to the full host (minimap2 -x map-ont, the pipeline's
     preset). All reads, not only the LQ/unmapped ones, so the regions do not
     depend on how the current classifier splits reads.
  2. Align the vector to the full host (minimap2 -x asm20 --secondary=yes, as
     the pipeline's homology step does), so vector-shared host regions are
     kept and the shared-region hybrid filter behaves as on the full host.
  3. Union every PAF target interval (primary, secondary, supplementary),
     pad by --pad bp on each side, clip to contig ends, and merge.
  4. Extract those regions with samtools faidx into host_subset.fasta.

The other references and the reads are copied unchanged. A samplesheet with
${projectDir}-relative paths, the region BED and a provenance README are
written alongside.

minimap2 and samtools run inside the pipeline container (--container), so
only Docker is needed on the host.

Coordinates: intervals are 0-based, half-open (PAF/BED) throughout. They are
converted to samtools' 1-based inclusive region strings only in
samtools_region(). Subset contig names use the BED coordinates:
<contig>_<start>_<end>.

Usage:
    python scripts/make_test_profile_data.py \\
        --reads reads.fastq.gz --vector vector.fasta --backbone backbone.fasta \\
        --helper helper.fasta --rep-cap rep_cap.fasta --host host.fasta \\
        --outdir test_data [--pad 5000] [--reads-paf r.paf --vector-paf v.paf]
"""

import argparse
import datetime
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterable, Iterator

DEFAULT_CONTAINER = "ghcr.io/iclah/aavision:1.0.2"
SAMPLE_NAME = "test_AAV9_CMV-GFP-Luc"


# ------------------------------------------------------------ pure helpers


def read_fai(path: str) -> dict[str, int]:
    """
    Read contig lengths from a samtools .fai index.

    Returns: {contig: length}. Raises ValueError if the index lists no contigs.
    """
    lengths = {}
    with open(path) as fh:
        for line in fh:
            if line.strip():
                name, length = line.split("\t")[:2]
                lengths[name] = int(length)
    if not lengths:
        raise ValueError(f"{path}: index lists no contigs")
    return lengths


def read_paf_intervals(paths: Iterable[str]) -> dict[str, list[tuple[int, int]]]:
    """
    Collect target intervals from one or more PAF files.

    Every record is kept (primary, secondary and supplementary), since any of
    them can become a read's best hit on a subset reference.

    Returns: {target contig: [(start, end), ...]}, 0-based half-open (PAF
    columns 8-9), in file order. Individual files may be empty; raises
    ValueError if all of them are, or if a line has fewer than 12 columns.
    """
    intervals: dict[str, list[tuple[int, int]]] = {}
    for path in paths:
        with open(path) as fh:
            for n, line in enumerate(fh, 1):
                if not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    raise ValueError(
                        f"{path}:{n}: not a PAF record ({len(cols)} columns)"
                    )
                intervals.setdefault(cols[5], []).append((int(cols[7]), int(cols[8])))
    if not intervals:
        raise ValueError(f"no alignments in {', '.join(paths)}; nothing to keep")
    return intervals


def merge_padded(
    intervals: dict[str, list[tuple[int, int]]],
    contig_lengths: dict[str, int],
    pad: int,
) -> list[tuple[str, int, int]]:
    """
    Pad each interval by `pad` bp per side, clip to [0, contig length), merge.

    Input and output are 0-based half-open. Overlapping and book-ended
    intervals are merged. Output is sorted by contig name, then start.
    Raises ValueError for a negative pad or a contig missing from
    contig_lengths (e.g. PAF and .fai from different references).
    """
    if pad < 0:
        raise ValueError(f"pad must be >= 0, got {pad}")
    merged = []
    for contig in sorted(intervals):
        if contig not in contig_lengths:
            raise ValueError(f"contig {contig} is not in the host index")
        length = contig_lengths[contig]
        padded = sorted(
            (max(0, s - pad), min(length, e + pad)) for s, e in intervals[contig]
        )
        cur_start, cur_end = padded[0]
        for start, end in padded[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                merged.append((contig, cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((contig, cur_start, cur_end))
    return merged


def samtools_region(contig: str, start: int, end: int) -> str:
    """
    Convert a 0-based half-open interval to a samtools region string
    (1-based inclusive). Contig names containing ':' are wrapped in braces,
    which samtools requires to parse them unambiguously.
    """
    name = f"{{{contig}}}" if ":" in contig else contig
    return f"{name}:{start + 1}-{end}"


def subset_contig_name(contig: str, start: int, end: int) -> str:
    """
    Name of a subset contig: <contig>_<start>_<end> with 0-based half-open
    (BED) coordinates. ':' is replaced by '_' so the name is a valid region
    prefix for samtools and IGV.
    """
    return f"{contig.replace(':', '_')}_{start}_{end}"


def rename_fasta_headers(
    lines: Iterable[str], regions: list[tuple[str, int, int]]
) -> Iterator[str]:
    """
    Rewrite samtools faidx headers (">chr1:1-10") to subset contig names.

    `regions` are the 0-based half-open intervals that were requested.
    Raises ValueError on a header that was not requested, or if any requested
    region is missing from the output.
    """
    names = {
        samtools_region(c, s, e): subset_contig_name(c, s, e) for c, s, e in regions
    }
    seen = set()
    for line in lines:
        if line.startswith(">"):
            header = line[1:].strip()
            if header not in names:
                raise ValueError(f"unexpected FASTA header from samtools: {header}")
            seen.add(header)
            yield f">{names[header]}\n"
        else:
            yield line
    missing = set(names) - seen
    if missing:
        raise ValueError(
            f"samtools output is missing {len(missing)} regions, e.g. {min(missing)}"
        )


# ------------------------------------------------------------ external tools


def run_in_container(image: str, mounts: list[str], workdir: str, cmd: str) -> str:
    """
    Run a shell command in the pipeline container as the current user.

    Each directory in `mounts` is bind-mounted at the same path; `workdir`
    (also the output location) is mounted read-write, the rest read-only.
    Returns stdout. Raises CalledProcessError (with stderr) on failure.
    """
    args = ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}"]
    for d in sorted(set(mounts)):
        if d != workdir:
            args += ["-v", f"{d}:{d}:ro"]
    args += ["-v", f"{workdir}:{workdir}", "-w", workdir, image, "sh", "-c", cmd]
    print(f"[docker] {cmd}", file=sys.stderr)
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout


def git_commit() -> str:
    """Current commit of the repository this script lives in, or 'unknown'."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=os.path.dirname(os.path.abspath(__file__)),
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return out.stdout.strip()


# ------------------------------------------------------------ main


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--reads", required=True, help="FASTQ(.gz) of test reads; copied unchanged"
    )
    p.add_argument("--vector", required=True)
    p.add_argument("--backbone", required=True)
    p.add_argument("--helper", required=True)
    p.add_argument("--rep-cap", required=True)
    p.add_argument(
        "--host", required=True, help="full host FASTA with a .fai index next to it"
    )
    p.add_argument("--outdir", required=True, help="output directory (e.g. test_data)")
    p.add_argument(
        "--pad",
        type=int,
        default=5000,
        help="bp kept on each side of every host hit (default 5000)",
    )
    p.add_argument("--threads", type=int, default=8)
    p.add_argument(
        "--split-index",
        default="4G",
        help="minimap2 -I; smaller fits less RAM (default 4G)",
    )
    p.add_argument("--container", default=DEFAULT_CONTAINER)
    p.add_argument(
        "--reads-paf", help="reuse an existing reads-vs-host PAF instead of aligning"
    )
    p.add_argument(
        "--vector-paf", help="reuse an existing vector-vs-host PAF instead of aligning"
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    inputs = {
        k: os.path.abspath(getattr(args, k))
        for k in ("reads", "vector", "backbone", "helper", "rep_cap", "host")
    }
    for name, path in inputs.items():
        if not os.path.isfile(path):
            sys.exit(f"--{name.replace('_', '-')}: {path} does not exist")
    fai = inputs["host"] + ".fai"
    if not os.path.isfile(fai):
        sys.exit(f"{fai} not found; run `samtools faidx {inputs['host']}` first")

    outdir = os.path.abspath(args.outdir)
    reads_dir = os.path.join(outdir, "reads")
    refs_dir = os.path.join(outdir, "references")
    os.makedirs(reads_dir, exist_ok=True)
    os.makedirs(refs_dir, exist_ok=True)

    work = tempfile.mkdtemp(prefix="test_profile_", dir=outdir)
    mounts = [os.path.dirname(p) for p in inputs.values()] + [work]
    host = inputs["host"]
    mm2 = (
        f"minimap2 -I {args.split_index} --split-prefix {work}/split -t {args.threads}"
    )
    # Commands as recorded in the README: basenames only, no local paths
    host_name = os.path.basename(host)
    reads_cmd = f"minimap2 -x map-ont -I {args.split_index} {host_name} {os.path.basename(inputs['reads'])}"
    vector_cmd = f"minimap2 -x asm20 --secondary=yes -I {args.split_index} {host_name} {os.path.basename(inputs['vector'])}"
    commands = [
        reads_cmd
        + (" (PAF from an earlier run of this command)" if args.reads_paf else ""),
        vector_cmd
        + (" (PAF from an earlier run of this command)" if args.vector_paf else ""),
    ]

    # 1-2. Host hits of the reads and the vector
    if args.reads_paf:
        reads_paf = os.path.abspath(args.reads_paf)
    else:
        reads_paf = os.path.join(work, "reads_vs_host.paf")
        cmd = f"{mm2} -x map-ont {host} {inputs['reads']} > {reads_paf}"
        run_in_container(args.container, mounts, work, cmd)
    if args.vector_paf:
        vector_paf = os.path.abspath(args.vector_paf)
    else:
        vector_paf = os.path.join(work, "vector_vs_host.paf")
        cmd = f"{mm2} -x asm20 --secondary=yes {host} {inputs['vector']} > {vector_paf}"
        run_in_container(args.container, mounts, work, cmd)

    # 3. Regions to keep
    intervals = read_paf_intervals([reads_paf, vector_paf])
    regions = merge_padded(intervals, read_fai(fai), args.pad)
    n_records = sum(len(v) for v in intervals.values())
    kept_bp = sum(e - s for _, s, e in regions)

    bed = os.path.join(refs_dir, "host_subset_regions.bed")
    with open(bed, "w") as fh:
        fh.writelines(
            f"{c}\t{s}\t{e}\t{subset_contig_name(c, s, e)}\n" for c, s, e in regions
        )

    # 4. Extract and rename
    region_file = os.path.join(work, "regions.txt")
    with open(region_file, "w") as fh:
        fh.writelines(samtools_region(c, s, e) + "\n" for c, s, e in regions)
    raw = os.path.join(work, "host_subset.raw.fasta")
    cmd = f"samtools faidx -r {region_file} {host} > {raw}"
    run_in_container(args.container, mounts, work, cmd)
    commands.append(f"samtools faidx -r regions.txt {host_name}")
    with (
        open(raw) as fin,
        open(os.path.join(refs_dir, "host_subset.fasta"), "w") as fout,
    ):
        fout.writelines(rename_fasta_headers(fin, regions))

    versions = run_in_container(
        args.container,
        mounts,
        work,
        "echo minimap2 $(minimap2 --version); samtools --version | head -1",
    ).strip()

    # Unchanged inputs
    shutil.copyfile(
        inputs["reads"], os.path.join(reads_dir, os.path.basename(inputs["reads"]))
    )
    for key, name in (
        ("vector", "vector"),
        ("backbone", "backbone"),
        ("helper", "helper"),
        ("rep_cap", "rep_cap"),
    ):
        shutil.copyfile(inputs[key], os.path.join(refs_dir, f"{name}.fasta"))

    ref = "${projectDir}/test_data/references"
    with open(os.path.join(outdir, "samplesheet.csv"), "w") as fh:
        fh.write("sample,fastq,vector,backbone,helper,host,rep_cap,spike_in\n")
        fh.write(
            f"{SAMPLE_NAME},${{projectDir}}/test_data/reads,{ref}/vector.fasta,{ref}/backbone.fasta,"
            f"{ref}/helper.fasta,{ref}/host_subset.fasta,{ref}/rep_cap.fasta,\n"
        )

    shutil.rmtree(work)

    contigs = sorted({c for c, _, _ in regions})
    summary = (
        f"{n_records} host alignment records on {len(contigs)} contigs -> "
        f"{len(regions)} regions, {kept_bp:,} bp (pad {args.pad} bp)"
    )
    write_readme(outdir, args, summary, commands, versions)

    print(summary)
    total = sum(
        os.path.getsize(os.path.join(r, f)) for r, _, fs in os.walk(outdir) for f in fs
    )
    print(f"{outdir}: {total / 1e6:.1f} MB")
    return 0


def write_readme(
    outdir: str,
    args: argparse.Namespace,
    summary: str,
    commands: list[str],
    versions: str,
) -> None:
    """Write test_data/README.md with the provenance of every file."""
    text = f"""# Test dataset for `-profile test`

Used by `nextflow run ... -profile test,docker` to check an installation.
It is a smoke test, not a validation: the host reference is a subset, so a
few reads can align differently than against the full genome, and the
numbers differ slightly from a full-host run.

| File | Content |
|------|---------|
| `samplesheet.csv` | one sample; paths use `${{projectDir}}`, which `main.nf` replaces with the pipeline directory |
| `reads/` | {os.path.basename(args.reads)}, copied unchanged |
| `references/vector,backbone,helper,rep_cap.fasta` | copied unchanged |
| `references/host_subset.fasta` | regions of GRCh38 that the reads or the vector align to |
| `references/host_subset_regions.bed` | those regions (0-based, half-open); column 4 is the contig name in `host_subset.fasta` |

## How the host subset was made

All reads and the vector were aligned to the full host; every alignment
(primary, secondary, supplementary) was padded by {args.pad} bp per side,
clipped to contig ends and merged.

{summary}

Commands (inside `{args.container}`):

```
{chr(10).join(commands)}
```

Tools: {versions.replace(chr(10), "; ")}
Generated {datetime.date.today().isoformat()} by `scripts/make_test_profile_data.py` at commit {git_commit()}.
"""
    with open(os.path.join(outdir, "README.md"), "w") as fh:
        fh.write(text)


if __name__ == "__main__":
    sys.exit(main())
