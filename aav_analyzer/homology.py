"""
Regions of the non-vector references that share sequence with the vector.

LQ vector reads that align to another reference only inside these regions are
pure vector, not hybrids (docs/DECISIONS.md, 2026-10-16). The regions are
found by aligning vector.fasta to each reference with minimap2 and kept as
reference intervals.

All intervals in this module are 0-based, half-open [start, end), as in BED,
PAF and pysam. Strand is ignored: a shared element counts in either
orientation.

Shared-region BED format (one file per sample; per-reference files may be
concatenated):
    # searched: <ref_name> [<ref_name> ...]
    <contig>\t<start>\t<end>\t<ref_name>
The '# searched:' line names every reference the homology step was run on,
so a reference with no shared region is distinguishable from one that was
never searched.
"""

import logging
import os
import subprocess
from dataclasses import dataclass, field

from .utils import merge_intervals

Interval = tuple[int, int]
ContigIntervals = dict[str, list[Interval]]

_SEARCHED_PREFIX = "# searched:"

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SharedRegions:
    """
    Shared regions per reference.

    Attributes:
        searched: Reference names the homology step was run on
        intervals: ref_name -> contig -> merged, sorted 0-based half-open intervals.
            A searched reference with no shared region may be absent.
    """

    searched: frozenset
    intervals: dict[str, ContigIntervals] = field(default_factory=dict)


def read_paf_target_intervals(paf_path: str) -> ContigIntervals:
    """
    Target (reference) intervals of every hit in a PAF file, merged per contig.

    PAF columns 6, 8, 9 are target name, start and end (0-based half-open).
    Primary and secondary hits are both kept, so every copy of a shared
    element is included. Strand (column 5) is ignored.

    Args:
        paf_path: PAF file from minimap2 (vector as query, reference as target)

    Returns:
        contig -> merged 0-based half-open intervals; empty dict if no hits
    """
    with open(paf_path) as fh:
        return _paf_intervals(fh)


def _paf_intervals(lines) -> ContigIntervals:
    """Merged target intervals per contig from PAF lines (see read_paf_target_intervals)."""
    raw: dict[str, list[Interval]] = {}
    for line in lines:
        if not line.strip():
            continue
        cols = line.split("\t")
        raw.setdefault(cols[5], []).append((int(cols[7]), int(cols[8])))
    return {contig: merge_intervals(ivs) for contig, ivs in raw.items()}


def find_shared_regions(
    vector_fasta: str, ref_fasta: str, preset: str = "asm20", threads: int = 1
) -> ContigIntervals:
    """
    Find the regions of ref_fasta that share sequence with vector_fasta.

    Runs `minimap2 -x <preset> --secondary=yes -c <ref> <vector>`, so hits are
    reported directly in reference coordinates, including secondary hits for
    elements present more than once in the reference.

    Args:
        vector_fasta: Vector reference FASTA (query)
        ref_fasta: Other reference FASTA (target), e.g. backbone or host
        preset: minimap2 preset; asm20 tolerates ~5-10% divergence
        threads: minimap2 threads

    Returns:
        contig -> merged 0-based half-open intervals on ref_fasta; empty if none
    """
    cmd = [
        "minimap2",
        "-x",
        preset,
        "--secondary=yes",
        "-c",
        "-t",
        str(threads),
        ref_fasta,
        vector_fasta,
    ]
    logger.info(
        f"Finding vector-shared regions in {os.path.basename(ref_fasta)}: {' '.join(cmd)}"
    )
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)

    regions = _paf_intervals(result.stdout.splitlines())

    total = sum(e - s for ivs in regions.values() for s, e in ivs)
    logger.info(
        f"{os.path.basename(ref_fasta)}: {sum(len(v) for v in regions.values())} "
        f"vector-shared region(s), {total} bp"
    )
    for contig, ivs in regions.items():
        for s, e in ivs:
            logger.info(f"  shared: {contig}:{s}-{e} ({e - s} bp, 0-based half-open)")
    return regions


def compute_shared_regions(
    vector_fasta: str, refs: dict[str, str], preset: str = "asm20", threads: int = 1
) -> SharedRegions:
    """
    Run find_shared_regions for every reference.

    Args:
        vector_fasta: Vector reference FASTA
        refs: ref_name -> FASTA path (e.g. backbone, helper, host, rep_cap)
        preset: minimap2 preset (see find_shared_regions)
        threads: minimap2 threads

    Returns:
        SharedRegions with every name in refs marked as searched; references
        without a shared region have no entry in intervals
    """
    intervals = {}
    for ref_name, ref_fasta in refs.items():
        found = find_shared_regions(
            vector_fasta, ref_fasta, preset=preset, threads=threads
        )
        if found:
            intervals[ref_name] = found
    return SharedRegions(searched=frozenset(refs), intervals=intervals)


def write_shared_bed(regions: SharedRegions, bed_path: str) -> None:
    """
    Write shared regions as BED (0-based half-open) with a '# searched:' header.

    Args:
        regions: Shared regions to write
        bed_path: Output BED path
    """
    with open(bed_path, "w") as fh:
        fh.write(f"{_SEARCHED_PREFIX} {' '.join(sorted(regions.searched))}\n")
        for ref_name in sorted(regions.intervals):
            for contig in sorted(regions.intervals[ref_name]):
                fh.writelines(
                    f"{contig}\t{start}\t{end}\t{ref_name}\n"
                    for start, end in regions.intervals[ref_name][contig]
                )


def read_shared_bed(bed_path: str) -> SharedRegions:
    """
    Read a shared-region BED written by write_shared_bed (or several concatenated).

    Args:
        bed_path: BED path (0-based half-open)

    Returns:
        SharedRegions with intervals merged per (ref_name, contig)

    Raises:
        ValueError: on a malformed line (fewer than 4 columns, non-integer or
            negative coordinates, start >= end) or a region for a reference
            not named in any '# searched:' line
    """
    searched = set()
    raw: dict[str, dict[str, list[Interval]]] = {}
    with open(bed_path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith(_SEARCHED_PREFIX):
                searched.update(line[len(_SEARCHED_PREFIX) :].split())
                continue
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                raise ValueError(
                    f"{bed_path}:{lineno}: expected 4 columns, got {len(cols)}: {line!r}"
                )
            contig, start_s, end_s, ref_name = cols[:4]
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                raise ValueError(
                    f"{bed_path}:{lineno}: non-integer coordinates: {line!r}"
                ) from None
            if start < 0 or start >= end:
                raise ValueError(
                    f"{bed_path}:{lineno}: need 0 <= start < end: {line!r}"
                )
            raw.setdefault(ref_name, {}).setdefault(contig, []).append((start, end))

    unsearched = sorted(set(raw) - searched)
    if unsearched:
        raise ValueError(
            f"{bed_path}: regions for reference(s) {unsearched} not listed in a "
            f"'{_SEARCHED_PREFIX}' line"
        )
    intervals = {
        ref: {contig: merge_intervals(ivs) for contig, ivs in contigs.items()}
        for ref, contigs in raw.items()
    }
    return SharedRegions(searched=frozenset(searched), intervals=intervals)


def unique_aligned_bases(segment, shared: ContigIntervals) -> int:
    """
    Reference bases aligned by one record that lie outside the shared regions.

    Aligned reference bases are the M/=/X blocks (pysam get_blocks, 0-based
    half-open); deletions, skips, insertions and clips do not count. Strand
    does not matter: blocks are in forward reference coordinates.

    Args:
        segment: Mapped pysam AlignedSegment
        shared: contig -> merged 0-based half-open intervals for this reference

    Returns:
        Number of aligned reference bases not in a shared interval
    """
    intervals = shared.get(segment.reference_name, [])
    unique = 0
    for b_start, b_end in segment.get_blocks():
        overlap = 0
        for s, e in intervals:
            if s >= b_end:
                break
            overlap += max(0, min(b_end, e) - max(b_start, s))
        unique += (b_end - b_start) - overlap
    return unique
