"""
Generate synthetic references and reads for the hybrid shared-region tests.

The vector shares one element with the backbone and one with the host, so
pure vector reads can align to those references. Writes into the output dir:
  vector.fa, backbone.fa, helper.fa, host.fa  - references
  reads.fastq        - reads with known structure, in their ORIGINAL orientation
  truth.tsv          - read_id, structure, strand, expected_category, read_length
  shared_truth.bed   - true shared regions on each reference (contig, start, end, ref_name)

Layout (0-based, half-open Python slices):
  vector    4,700 bp  V1[0,1000) E_bb[1000,1800) V2[1800,3000) E_host[3000,3500) V3[3500,4700)
  backbone  3,800 bp  B1[0,1500) E_bb[1500,2300) B2[2300,3800)
  helper    4,000 bp  random, nothing shared
  host      host_chr1 = H1[0,2000) E_host[2000,2500) H2[2500,4500); host_chr2 = 3,000 bp random

expected_category is one of:
  HQ               primary vector alignment covers >= 90% of the read (never enters the hybrid pass)
  LQ               LQ vector read with no eligible non-vector alignment (pure LQ)
  hybrid_<ref>     LQ vector read with >= 100 aligned bases on <ref> outside its shared regions

Reads labelled strand '-' are the reverse complement of the '+' construct.

Usage:
    python scripts/make_hybrid_testdata.py --outdir <dir> [--seed 20261017]
"""

import argparse
import os
import random

SUB_RATE = 0.01  # substitution-only errors keep read lengths exact
READS_PER_CLASS = 6

# Shared regions in reference coordinates (0-based, half-open)
BACKBONE_SHARED = (1500, 2300)
HOST_SHARED = (2000, 2500)

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    """Reverse complement of an uppercase ACGT sequence."""
    return seq[::-1].translate(_COMP)


def _random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(rng: random.Random, seq: str, rate: float) -> str:
    out = []
    for base in seq:
        if rng.random() < rate:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        else:
            out.append(base)
    return "".join(out)


def _references(rng: random.Random) -> dict:
    """Build the references; returns {name: {contig: sequence}}."""
    e_bb = _random_seq(rng, 800)
    e_host = _random_seq(rng, 500)
    vector = (
        _random_seq(rng, 1000)
        + e_bb
        + _random_seq(rng, 1200)
        + e_host
        + _random_seq(rng, 1200)
    )
    backbone = _random_seq(rng, 1500) + e_bb + _random_seq(rng, 1500)
    helper = _random_seq(rng, 4000)
    host1 = _random_seq(rng, 2000) + e_host + _random_seq(rng, 2000)
    host2 = _random_seq(rng, 3000)
    return {
        "vector": {"vector": vector},
        "backbone": {"backbone": backbone},
        "helper": {"helper": helper},
        "host": {"host_chr1": host1, "host_chr2": host2},
    }


def _constructs(rng: random.Random, refs: dict):
    """
    Yield (structure, expected_category, '+'-strand sequence) per read.

    Every non-HQ construct has a primary vector alignment well below 90% of
    the read, so it is LQ and enters the hybrid pass.
    """
    v = refs["vector"]["vector"]
    bb = refs["backbone"]["backbone"]
    hp = refs["helper"]["helper"]
    h1 = refs["host"]["host_chr1"]

    for _ in range(READS_PER_CLASS):
        j = rng.randint(0, 50)

        # Control: full-length genome
        yield "full_length", "HQ", v

        # Snapback covering E_bb: pure vector, aligns to backbone only inside
        # the shared region (the v1.0.1 bug: was called hybrid_backbone)
        a = v[800 + j : 2000]
        yield "snapback_shared_bb", "LQ", a + revcomp(a)

        # Snapback covering E_host: pure vector, aligns to host_chr1 only inside
        # the shared region
        a = v[2600 + j : 3800]
        yield "snapback_shared_host", "LQ", a + revcomp(a)

        # Vector V3 + 1.5 kb unique backbone: hybrid, nothing shared involved
        yield "vec_plus_bb_unique", "hybrid_backbone", v[3500 + j :] + bb[2300:3800]

        # Vector up to the end of E_bb, continuing into 300 bp of unique
        # backbone, then junk: the backbone alignment spans the shared/unique
        # boundary with ~300 unique bases -> hybrid
        yield (
            "vec_span_bb_300",
            "hybrid_backbone",
            v[j:1800] + bb[2300:2600] + _random_seq(rng, 1000),
        )

        # Same, but only 60 bp of unique backbone: below the 100 bp threshold
        yield "vec_span_bb_60", "LQ", v[j:1800] + bb[2300:2360] + _random_seq(rng, 1000)

        # Truncated vector covering E_bb + 1.5 kb foreign sequence: LQ (coverage
        # ~0.45, not a snapback), aligns to backbone only inside the shared
        # region -> pure LQ (the v1.0.1 bug, same-strand variant)
        yield "vec_shared_plus_junk", "LQ", v[800 + j : 2000] + _random_seq(rng, 1500)

        # Vector + unique helper: hybrid, helper has no shared region
        yield "vec_plus_helper", "hybrid_helper", v[j:1000] + hp[1000:2500]

        # Vector V3 + unique host_chr1 sequence: hybrid
        yield "vec_plus_host_unique", "hybrid_host", v[3500 + j :] + h1[2500:4000]


def _write_fasta(path: str, contigs: dict) -> None:
    with open(path, "w") as fa:
        fa.writelines(f">{name}\n{seq}\n" for name, seq in contigs.items())


def generate(outdir: str, seed: int = 20261017) -> dict:
    """
    Write the references, reads.fastq, truth.tsv and shared_truth.bed into outdir.

    Args:
        outdir: Output directory (created if missing)
        seed: Random seed; the same seed always gives identical files

    Returns:
        Dict mapping 'vector', 'backbone', 'helper', 'host', 'reads', 'truth'
        and 'shared_truth' to file paths
    """
    rng = random.Random(seed)
    os.makedirs(outdir, exist_ok=True)
    refs = _references(rng)

    paths = {name: os.path.join(outdir, f"{name}.fa") for name in refs}
    paths["reads"] = os.path.join(outdir, "reads.fastq")
    paths["truth"] = os.path.join(outdir, "truth.tsv")
    paths["shared_truth"] = os.path.join(outdir, "shared_truth.bed")

    for name, contigs in refs.items():
        _write_fasta(paths[name], contigs)

    with open(paths["shared_truth"], "w") as bed:
        bed.write(f"backbone\t{BACKBONE_SHARED[0]}\t{BACKBONE_SHARED[1]}\tbackbone\n")
        bed.write(f"host_chr1\t{HOST_SHARED[0]}\t{HOST_SHARED[1]}\thost\n")

    with open(paths["reads"], "w") as fq, open(paths["truth"], "w") as tt:
        tt.write("read_id\tstructure\tstrand\texpected_category\tread_length\n")
        n = 0
        for structure, expected, seq in _constructs(rng, refs):
            for strand in ("+", "-"):
                read = seq if strand == "+" else revcomp(seq)
                read = _mutate(rng, read, SUB_RATE)
                read_id = f"r{n:04d}_{structure}_{'fwd' if strand == '+' else 'rev'}"
                n += 1
                fq.write(f"@{read_id}\n{read}\n+\n{'I' * len(read)}\n")
                tt.write(f"{read_id}\t{structure}\t{strand}\t{expected}\t{len(read)}\n")

    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--seed", type=int, default=20261017, help="Random seed")
    args = parser.parse_args()
    paths = generate(args.outdir, args.seed)
    for name, path in paths.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
