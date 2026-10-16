"""
Tests for the read-side check in hybrid calls (docs/DECISIONS.md).

A reference is eligible for an LQ read only if at least min_unique_bases
distinct READ positions are aligned to it outside its vector-shared regions
AND are not covered by any primary/supplementary vector record of the read.

Every read here is 1000 bases long. Read intervals are 0-based half-open in
original read orientation; reference intervals are 0-based half-open.
Backbone contig 'bb' (1000 bp) has shared region [0, 500).
"""

import argparse
import csv
import json
import sys

import pysam
import pytest
from test_hybrid_homology import (
    _categories,
    _header,
    _segment,
    _write_bam,
    read_fastq_names,
    write_fastq,
)

from aav_analyzer import alignment, classification, cli, homology

SHARED = homology.SharedRegions(
    searched=frozenset({"backbone"}),
    intervals={"backbone": {"bb": [(0, 500)]}},
)
# Every read: vector covers read [0, 600)
VECTOR_IV = [(0, 600)]


# ----------------------------------------------------------------------
# homology.novel_read_positions
# ----------------------------------------------------------------------


def test_novel_positions_plus_strand_outside_vector():
    # read [600, 900) -> ref [500, 800), outside shared and outside vector
    seg = _segment(_header(), "r", 500, "600S300M100S")
    pos = homology.novel_read_positions(seg, {"bb": [(0, 500)]}, VECTOR_IV)
    assert pos == set(range(600, 900))


def test_novel_positions_inside_vector_part_are_not_novel():
    # read [100, 400) is already explained by the vector
    seg = _segment(_header(), "r", 500, "100S300M600S")
    assert homology.novel_read_positions(seg, {"bb": [(0, 500)]}, VECTOR_IV) == set()


def test_novel_positions_minus_strand_flipped():
    # Same CIGAR as above but minus strand: SEQ [100, 400) is read [600, 900)
    seg = _segment(_header(), "r", 500, "100S300M600S", flag=16)
    pos = homology.novel_read_positions(seg, {"bb": [(0, 500)]}, VECTOR_IV)
    assert pos == set(range(600, 900))


def test_novel_positions_hard_clip_offsets_read_coordinates():
    # Leading 600H: aligned bases are read [600, 900) even though SEQ starts at 0
    seg = _segment(_header(), "r", 500, "600H300M100H")
    pos = homology.novel_read_positions(seg, {"bb": [(0, 500)]}, VECTOR_IV)
    assert pos == set(range(600, 900))


def test_novel_positions_shared_reference_part_excluded():
    # read [600, 900) -> ref [350, 650): ref [350, 500) is shared
    seg = _segment(_header(), "r", 350, "600S300M100S")
    pos = homology.novel_read_positions(seg, {"bb": [(0, 500)]}, VECTOR_IV)
    assert pos == set(range(750, 900))


def test_novel_positions_insertions_and_deletions():
    # read [600,650) M, [650,700) I (no ref base: not counted), [700,750) M after a 20D
    seg = _segment(_header(), "r", 500, "600S50M50I20D50M250S")
    pos = homology.novel_read_positions(seg, {}, VECTOR_IV)
    assert pos == set(range(600, 650)) | set(range(700, 750))


def test_novel_positions_no_vector_intervals():
    seg = _segment(_header(), "r", 500, "100S300M600S")
    assert len(homology.novel_read_positions(seg, {}, [])) == 300


# ----------------------------------------------------------------------
# alignment.vector_read_intervals
# ----------------------------------------------------------------------


def _write_vector_bam(path, records):
    """records: list of (name, ref_start, cigar, flag). Vector contig 'vec'."""
    header = _header(("vec",))
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, start, cigar, flag in records:
            if flag & 4:
                seg = pysam.AlignedSegment(header)
                seg.query_name, seg.flag = name, flag
                seg.query_sequence = "A" * 1000
                seg.query_qualities = pysam.qualitystring_to_array("I" * 1000)
                out.write(seg)
            else:
                out.write(_segment(header, name, start, cigar, flag=flag))
    return str(path)


def test_vector_read_intervals_primary_supp_merged_secondary_skipped(tmp_path):
    bam = _write_vector_bam(
        tmp_path / "vector.bam",
        [
            ("a", 0, "500M500S", 0),  # read [0, 500)
            ("a", 0, "450S100M450S", 2048),  # read [450, 550): overlaps, merged
            ("a", 0, "900S100M", 256),  # secondary: ignored
            ("b", 0, "100S300M600S", 16),  # minus: read [600, 900)
            ("c", 0, "", 4),  # unmapped: no entry
        ],
    )
    got = alignment.vector_read_intervals(bam)
    assert got == {"a": [(0, 550)], "b": [(600, 900)]}


def test_vector_read_intervals_restricted_to_names(tmp_path):
    bam = _write_vector_bam(
        tmp_path / "vector.bam", [("a", 0, "500M500S", 0), ("b", 0, "500M500S", 0)]
    )
    assert alignment.vector_read_intervals(bam, read_names={"b"}) == {"b": [(0, 500)]}


# ----------------------------------------------------------------------
# classify_reads_from_bams with the read-side check
# ----------------------------------------------------------------------


@pytest.fixture
def case(tmp_path):
    """
    All reads: vector covers read [0, 600). Backbone ref bb shared [0, 500).
    Expected with min_unique_bases=100 and the read-side check on:
      junction       read [600,900) -> bb [500,800)               -> backbone (300 novel)
      vector_part    read [100,400) -> bb [500,800)               -> unmapped (0 novel;
                                                                     ref-side alone: 300)
      minus_junction minus, SEQ [100,400) = read [600,900)        -> backbone
      fragments      5 x 40M on read [0,600), bb outside shared   -> unmapped (ref-side: 200)
      double_count   primary + supp both on read [600,660)        -> unmapped (60 distinct;
                                                                     ref-side: 120)
      shared_novel   read [600,900) -> bb [100,400) (shared)      -> unmapped
      at_threshold   read [600,700) -> bb [500,600)               -> backbone (exactly 100)
      below          read [600,699) -> bb [500,599)               -> unmapped (99)
      no_alignment   FASTQ only                                   -> unmapped
    """
    recs = [
        ("junction", 500, "600S300M100S", 0, 300),
        ("vector_part", 500, "100S300M600S", 0, 300),
        ("minus_junction", 500, "100S300M600S", 16, 300),
        ("double_count", 500, "600S60M340S", 0, 60),
        ("double_count", 700, "600S60M340S", 2048, 60),
        ("shared_novel", 100, "600S300M100S", 0, 300),
        ("at_threshold", 500, "600S100M300S", 0, 100),
        ("below", 500, "600S99M301S", 0, 99),
    ]
    for i in range(5):
        flag = 0 if i == 0 else 2048
        recs.append(
            ("fragments", 500 + 100 * i, f"{100 * i}S40M{960 - 100 * i}S", flag, 40)
        )
    bam = _write_bam(tmp_path / "alignments_hybrid_backbone.bam", "bb", recs)
    names = sorted({r[0] for r in recs} | {"no_alignment"})
    fastq = tmp_path / "vector_LQ.fastq"
    write_fastq(fastq, names, length=1000)
    vec = _write_vector_bam(
        tmp_path / "vector.bam", [(n, 0, "600M400S", 0) for n in names]
    )
    return {
        "bams": {"backbone": bam},
        "fastq": str(fastq),
        "prefix": str(tmp_path / "hybrid"),
        "vector_bam": vec,
        "vector_iv": {n: list(VECTOR_IV) for n in names},
    }


def _classify(c, vector_intervals="from_case", min_unique=100):
    vi = c["vector_iv"] if vector_intervals == "from_case" else vector_intervals
    return classification.classify_reads_from_bams(
        c["bams"],
        c["fastq"],
        c["prefix"],
        "hybrid",
        shared_regions=SHARED,
        min_unique_bases=min_unique,
        vector_intervals=vi,
    )


def test_read_side_check_calls(case):
    got = _categories(_classify(case))
    assert got["backbone"] == {"junction", "minus_junction", "at_threshold"}
    assert got["unmapped"] == {
        "vector_part",
        "fragments",
        "double_count",
        "shared_novel",
        "below",
        "no_alignment",
    }


def test_without_read_side_check_is_reference_side_only(case):
    # vector_intervals=None keeps the 29fc7e6 behaviour
    got = _categories(_classify(case, vector_intervals=None))
    assert {"vector_part", "fragments", "double_count"} <= got["backbone"]
    assert "shared_novel" in got["unmapped"]


def test_read_side_report_columns_and_stats(case):
    _classify(case)
    with open(case["prefix"] + "_shared_only.tsv") as fh:
        rows = {r["read_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    assert rows["vector_part"]["unique_bases"] == "300"  # reference side
    assert rows["vector_part"]["novel_read_bases"] == "0"  # read side
    assert rows["double_count"]["unique_bases"] == "120"
    assert rows["double_count"]["novel_read_bases"] == "60"
    assert rows["below"]["novel_read_bases"] == "99"
    with open(case["prefix"] + "_filter_stats.json") as fh:
        stats = json.load(fh)
    assert stats["read_side_check"] is True
    assert stats["rescued_to_lq"] == 5


def test_report_without_read_side_check_has_empty_novel_column(case):
    _classify(case, vector_intervals=None)
    with open(case["prefix"] + "_shared_only.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows and all(r["novel_read_bases"] == "" for r in rows)
    with open(case["prefix"] + "_filter_stats.json") as fh:
        assert json.load(fh)["read_side_check"] is False


def test_read_with_hits_but_no_vector_record_fails(case):
    vi = dict(case["vector_iv"])
    del vi["junction"]
    with pytest.raises(ValueError, match="junction"):
        _classify(case, vector_intervals=vi)


def test_read_side_check_requires_shared_regions(case):
    with pytest.raises(ValueError, match="shared"):
        classification.classify_reads_from_bams(
            case["bams"],
            case["fastq"],
            case["prefix"],
            "hybrid",
            vector_intervals=case["vector_iv"],
        )


# ----------------------------------------------------------------------
# CLI: classify-best --vector-sam
# ----------------------------------------------------------------------


def _args(c, bed, vector_sam, mode="hybrid"):
    return argparse.Namespace(
        bams=[f"{n}={p}" for n, p in c["bams"].items()],
        source_fastq=c["fastq"],
        out_prefix=c["prefix"],
        mode=mode,
        shared_bed=bed,
        min_unique_bases=100,
        vector_sam=vector_sam,
    )


def test_cli_classify_best_with_vector_sam(case, tmp_path):
    bed = tmp_path / "shared.bed"
    homology.write_shared_bed(SHARED, str(bed))
    assert cli.cmd_classify_best(_args(case, str(bed), case["vector_bam"])) == 0
    assert read_fastq_names(case["prefix"] + "_backbone.fastq") == {
        "junction",
        "minus_junction",
        "at_threshold",
    }


def test_cli_vector_sam_rejected_for_impurity(case):
    with pytest.raises(ValueError, match="hybrid"):
        cli.cmd_classify_best(_args(case, None, case["vector_bam"], mode="impurity"))


def test_cli_vector_sam_requires_shared_bed(case):
    with pytest.raises(ValueError, match="shared-bed"):
        cli.cmd_classify_best(_args(case, None, case["vector_bam"]))


def test_classify_best_parser_has_vector_sam(monkeypatch):
    import aav_analyze

    base = [
        "aav_analyze.py",
        "classify-best",
        "--bams",
        "backbone=b.bam",
        "--source-fastq",
        "s.fq",
        "--out-prefix",
        "p",
    ]
    monkeypatch.setattr(sys, "argv", base)
    assert aav_analyze.parse_arguments().vector_sam is None
    monkeypatch.setattr(sys, "argv", base + ["--vector-sam", "v.sam"])
    assert aav_analyze.parse_arguments().vector_sam == "v.sam"
