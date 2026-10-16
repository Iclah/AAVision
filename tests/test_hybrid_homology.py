"""
Tests for excluding vector-shared reference regions from hybrid calls.

Fast tests build small PAF/BED/BAM files by hand. Slow tests run minimap2 on
the synthetic references and reads from scripts/make_hybrid_testdata.py and
compare to its truth table.

All intervals in this file are 0-based, half-open (BED/PAF/pysam convention).
"""

import argparse
import csv
import json
import os
import shutil
import sys

import pysam
import pytest

from aav_analyzer import alignment, classification, cli, homology

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
import make_hybrid_testdata

REF_LEN = 1000


def read_fastq_names(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return set()
    with pysam.FastxFile(path) as fq:
        return {entry.name for entry in fq}


def write_fastq(path, names, length=300):
    with open(path, "w") as fh:
        fh.writelines(f"@{n}\n{'A' * length}\n+\n{'I' * length}\n" for n in names)


# ----------------------------------------------------------------------
# merge_intervals
# ----------------------------------------------------------------------


def test_merge_overlapping_and_unsorted():
    assert homology.merge_intervals([(50, 80), (0, 10), (5, 20)]) == [(0, 20), (50, 80)]


def test_merge_adjacent_half_open_intervals():
    # [0,10) and [10,20) touch with no gap: one interval
    assert homology.merge_intervals([(0, 10), (10, 20)]) == [(0, 20)]


def test_merge_keeps_one_base_gap():
    assert homology.merge_intervals([(0, 10), (11, 20)]) == [(0, 10), (11, 20)]


def test_merge_empty():
    assert homology.merge_intervals([]) == []


# ----------------------------------------------------------------------
# PAF parsing
# ----------------------------------------------------------------------


def _paf_line(tname, tstart, tend, strand="+", tp="P"):
    # qname qlen qstart qend strand tname tlen tstart tend nmatch alen mapq tags
    return (
        f"vector\t4700\t0\t{tend - tstart}\t{strand}\t{tname}\t5000\t{tstart}\t{tend}"
        f"\t{tend - tstart}\t{tend - tstart}\t60\ttp:A:{tp}\n"
    )


def test_paf_target_intervals_both_strands_secondary_and_contigs(tmp_path):
    paf = tmp_path / "hits.paf"
    paf.write_text(
        _paf_line("chr1", 100, 200)
        + _paf_line("chr1", 150, 300, strand="-")  # overlaps, minus strand
        + _paf_line("chr1", 1000, 1100, tp="S")  # secondary copy: kept
        + _paf_line("chr2", 0, 50)
    )
    assert homology.read_paf_target_intervals(str(paf)) == {
        "chr1": [(100, 300), (1000, 1100)],
        "chr2": [(0, 50)],
    }


def test_paf_empty_file(tmp_path):
    paf = tmp_path / "empty.paf"
    paf.write_text("")
    assert homology.read_paf_target_intervals(str(paf)) == {}


# ----------------------------------------------------------------------
# Shared-region BED round trip and validation
# ----------------------------------------------------------------------


def test_bed_round_trip(tmp_path):
    regions = homology.SharedRegions(
        searched=frozenset({"backbone", "helper", "host"}),
        intervals={
            "backbone": {"backbone": [(1500, 2300)]},
            "host": {"host_chr1": [(2000, 2500)], "host_chr2": [(0, 10)]},
        },
    )
    bed = tmp_path / "shared.bed"
    homology.write_shared_bed(regions, str(bed))

    lines = bed.read_text().splitlines()
    assert lines[0] == "# searched: backbone helper host"
    assert "backbone\t1500\t2300\tbackbone" in lines

    assert homology.read_shared_bed(str(bed)) == regions


def test_bed_concatenated_files_merge(tmp_path):
    """Per-reference BEDs concatenated (as the Nextflow fan-in does) read as one."""
    a = homology.SharedRegions(frozenset({"backbone"}), {"backbone": {"bb": [(0, 10)]}})
    b = homology.SharedRegions(frozenset({"helper"}), {})
    pa, pb = tmp_path / "a.bed", tmp_path / "b.bed"
    homology.write_shared_bed(a, str(pa))
    homology.write_shared_bed(b, str(pb))
    combined = tmp_path / "all.bed"
    combined.write_text(pa.read_text() + pb.read_text())

    got = homology.read_shared_bed(str(combined))
    assert got.searched == frozenset({"backbone", "helper"})
    assert got.intervals == {"backbone": {"bb": [(0, 10)]}}


@pytest.mark.parametrize(
    "line",
    [
        "bb\t100\t100\tbackbone",  # empty interval
        "bb\t200\t100\tbackbone",  # start > end
        "bb\t-1\t100\tbackbone",  # negative
        "bb\tx\t100\tbackbone",  # not an integer
        "bb\t0\t100",  # missing ref_name column
    ],
)
def test_bed_malformed_line_fails(tmp_path, line):
    bed = tmp_path / "bad.bed"
    bed.write_text(f"# searched: backbone\n{line}\n")
    with pytest.raises(ValueError):
        homology.read_shared_bed(str(bed))


def test_bed_region_for_unsearched_ref_fails(tmp_path):
    bed = tmp_path / "bad.bed"
    bed.write_text("# searched: backbone\nh\t0\t10\thelper\n")
    with pytest.raises(ValueError):
        homology.read_shared_bed(str(bed))


# ----------------------------------------------------------------------
# unique_aligned_bases
# ----------------------------------------------------------------------


def _header(contigs=("bb",)):
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": REF_LEN} for c in contigs]}
    )


def _segment(header, name, ref_start, cigar, flag=0, contig_id=0, score=100):
    seg = pysam.AlignedSegment(header)
    seg.query_name = name
    seg.flag = flag
    seg.reference_id = contig_id
    seg.reference_start = ref_start
    seg.mapping_quality = 60
    seg.cigarstring = cigar
    qlen = sum(n for op, n in seg.cigartuples if op in (0, 1, 4, 7, 8))
    seg.query_sequence = "A" * qlen
    seg.query_qualities = pysam.qualitystring_to_array("I" * qlen)
    seg.set_tag("AS", score)
    return seg


def test_unique_bases_excludes_shared_and_deletions():
    # Blocks [100,150) and [160,210); shared [120,170) overlaps 30 + 10 bases
    seg = _segment(_header(), "q", 100, "20S50M10D50M")
    assert homology.unique_aligned_bases(seg, {"bb": [(120, 170)]}) == 60


def test_unique_bases_is_strand_independent():
    plus = _segment(_header(), "q", 100, "50M10D50M")
    minus = _segment(_header(), "q", 100, "50M10D50M", flag=16)
    shared = {"bb": [(120, 170)]}
    assert homology.unique_aligned_bases(plus, shared) == 60
    assert homology.unique_aligned_bases(minus, shared) == 60


def test_unique_bases_insertions_do_not_count():
    # Insertions consume no reference: 100M + 50I aligns 100 reference bases
    seg = _segment(_header(), "q", 0, "50M50I50M")
    assert homology.unique_aligned_bases(seg, {}) == 100


def test_unique_bases_shared_on_other_contig_ignored():
    seg = _segment(_header(("c1", "c2")), "q", 0, "100M", contig_id=1)
    assert homology.unique_aligned_bases(seg, {"c1": [(0, 100)]}) == 100


def test_unique_bases_fully_shared():
    seg = _segment(_header(), "q", 100, "100M")
    assert homology.unique_aligned_bases(seg, {"bb": [(0, 500)]}) == 0


# ----------------------------------------------------------------------
# classify_reads_from_bams with shared regions
# ----------------------------------------------------------------------


def _write_bam(path, contig, records):
    """records: list of (name, ref_start, cigar, flag, AS). Sorted and indexed."""
    header = _header((contig,))
    unsorted = str(path) + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for name, start, cigar, flag, score in records:
            out.write(_segment(header, name, start, cigar, flag=flag, score=score))
    pysam.sort("-o", str(path), unsorted)
    pysam.index(str(path))
    os.remove(unsorted)
    return str(path)


@pytest.fixture
def bam_case(tmp_path):
    """
    backbone contig 'bb' has shared region [0,500); helper has none.

    Reads and expected category with min_unique_bases=100:
      shared_only     bb [0,300)                          -> unmapped (rescued to LQ)
      spans           bb [400,700)                        -> backbone (200 unique)
      at_threshold    bb [400,600)                        -> backbone (exactly 100, inclusive)
      below           bb [399,599)                        -> unmapped (99 unique)
      supp_sum        bb [440,520) + supp [600,690)       -> backbone (20 + 90 = 110)
      minus_spans     bb [400,700) on - strand            -> backbone
      shared_beats    bb [0,500) AS 500; helper AS 200    -> helper (only eligible ref)
      no_alignment    in FASTQ only                       -> unmapped
    """
    bb = _write_bam(
        tmp_path / "alignments_hybrid_backbone.bam",
        "bb",
        [
            ("shared_only", 0, "300M", 0, 300),
            ("spans", 400, "300M", 0, 300),
            ("at_threshold", 400, "200M", 0, 200),
            ("below", 399, "200M", 0, 200),
            ("supp_sum", 440, "80M90S", 0, 80),
            ("supp_sum", 600, "80S90M", 2048, 90),
            ("minus_spans", 400, "300M", 16, 300),
            ("shared_beats", 0, "500M", 0, 500),
        ],
    )
    hp = _write_bam(
        tmp_path / "alignments_hybrid_helper.bam",
        "hp",
        [("shared_beats", 100, "200M", 0, 200)],
    )
    names = [
        "shared_only",
        "spans",
        "at_threshold",
        "below",
        "supp_sum",
        "minus_spans",
        "shared_beats",
        "no_alignment",
    ]
    fastq = tmp_path / "vector_LQ.fastq"
    write_fastq(fastq, names)
    shared = homology.SharedRegions(
        searched=frozenset({"backbone", "helper"}),
        intervals={"backbone": {"bb": [(0, 500)]}},
    )
    return {
        "bams": {"backbone": bb, "helper": hp},
        "fastq": str(fastq),
        "prefix": str(tmp_path / "hybrid"),
        "shared": shared,
    }


def _categories(output_files):
    return {cat: read_fastq_names(path) for cat, path in output_files.items()}


def test_shared_regions_filter_hybrids(bam_case):
    c = bam_case
    out = classification.classify_reads_from_bams(
        c["bams"],
        c["fastq"],
        c["prefix"],
        "hybrid",
        shared_regions=c["shared"],
        min_unique_bases=100,
    )
    got = _categories(out)
    assert got["backbone"] == {"spans", "at_threshold", "supp_sum", "minus_spans"}
    assert got["helper"] == {"shared_beats"}
    assert got["unmapped"] == {"shared_only", "below", "no_alignment"}


def test_shared_only_report_and_stats(bam_case):
    c = bam_case
    classification.classify_reads_from_bams(
        c["bams"],
        c["fastq"],
        c["prefix"],
        "hybrid",
        shared_regions=c["shared"],
        min_unique_bases=100,
    )
    with open(c["prefix"] + "_shared_only.tsv") as fh:
        rows = {
            (r["read_id"], r["ref_name"]): r for r in csv.DictReader(fh, delimiter="\t")
        }
    # Every rejected (read, ref) pair is listed, including shared_beats whose
    # backbone alignment was rejected although the read went to helper
    assert set(rows) == {
        ("shared_only", "backbone"),
        ("below", "backbone"),
        ("shared_beats", "backbone"),
    }
    assert rows[("below", "backbone")]["unique_bases"] == "99"
    assert rows[("below", "backbone")]["best_AS"] == "200"
    assert rows[("shared_beats", "backbone")]["final_category"] == "helper"
    assert rows[("shared_only", "backbone")]["final_category"] == "unmapped"

    with open(c["prefix"] + "_filter_stats.json") as fh:
        stats = json.load(fh)
    assert stats["min_unique_bases"] == 100
    assert stats["rejected_alignments_by_ref"] == {"backbone": 3, "helper": 0}
    # Reads that became pure LQ only because of the filter, by the ref of
    # their best rejected alignment
    assert stats["rescued_to_lq_by_ref"] == {"backbone": 2, "helper": 0}
    assert stats["rescued_to_lq"] == 2


def test_without_shared_regions_matches_v1_0_1(bam_case):
    c = bam_case
    out = classification.classify_reads_from_bams(
        c["bams"], c["fastq"], c["prefix"], "hybrid"
    )
    got = _categories(out)
    assert got["backbone"] == {
        "shared_only",
        "spans",
        "at_threshold",
        "below",
        "supp_sum",
        "minus_spans",
        "shared_beats",
    }
    assert got["helper"] == set()
    assert got["unmapped"] == {"no_alignment"}
    # No filter outputs when the filter is off
    assert not os.path.exists(c["prefix"] + "_shared_only.tsv")
    assert not os.path.exists(c["prefix"] + "_filter_stats.json")


def test_bam_ref_not_searched_fails(bam_case):
    c = bam_case
    partial = homology.SharedRegions(frozenset({"backbone"}), c["shared"].intervals)
    with pytest.raises(ValueError, match="helper"):
        classification.classify_reads_from_bams(
            c["bams"], c["fastq"], c["prefix"], "hybrid", shared_regions=partial
        )


def test_empty_inputs_with_filter(tmp_path):
    bam = _write_bam(tmp_path / "alignments_hybrid_backbone.bam", "bb", [])
    fastq = tmp_path / "empty.fastq"
    fastq.write_text("")
    prefix = str(tmp_path / "hybrid")
    shared = homology.SharedRegions(
        frozenset({"backbone"}), {"backbone": {"bb": [(0, 10)]}}
    )
    out = classification.classify_reads_from_bams(
        {"backbone": bam}, str(fastq), prefix, "hybrid", shared_regions=shared
    )
    assert all(os.path.getsize(p) == 0 for p in out.values())
    with open(prefix + "_filter_stats.json") as fh:
        stats = json.load(fh)
    assert stats["rescued_to_lq"] == 0


# ----------------------------------------------------------------------
# End to end with minimap2 on synthetic references and reads
# ----------------------------------------------------------------------

_TOOLS = ("minimap2", "samtools", "seqkit")
_needs_tools = pytest.mark.skipif(
    any(shutil.which(t) is None for t in _TOOLS),
    reason=f"needs {', '.join(_TOOLS)} on PATH",
)
_OTHER_REFS = ("backbone", "helper", "host")


def _truth_bed(path):
    truth = {}
    with open(path) as fh:
        for line in fh:
            contig, start, end, ref = line.rstrip("\n").split("\t")
            truth.setdefault(ref, {}).setdefault(contig, []).append(
                (int(start), int(end))
            )
    return truth


@pytest.mark.slow
@_needs_tools
def test_find_shared_regions_matches_truth(tmp_path):
    paths = make_hybrid_testdata.generate(str(tmp_path / "data"))
    truth = _truth_bed(paths["shared_truth"])
    for ref in _OTHER_REFS:
        found = homology.find_shared_regions(
            paths["vector"], paths[ref], preset="asm20"
        )
        expected = truth.get(ref, {})
        assert set(found) == set(expected), ref
        for contig, intervals in expected.items():
            assert len(found[contig]) == len(intervals), (ref, contig)
            for (fs, fe), (ts, te) in zip(found[contig], intervals):
                assert abs(fs - ts) <= 20 and abs(fe - te) <= 20, (ref, contig, fs, fe)


@pytest.mark.slow
@_needs_tools
def test_end_to_end_hybrids_match_truth_table(tmp_path):
    paths = make_hybrid_testdata.generate(str(tmp_path / "data"))
    with open(paths["truth"]) as fh:
        truth = {row["read_id"]: row for row in csv.DictReader(fh, delimiter="\t")}

    sam = str(tmp_path / "vector.sam")
    vector_fastq = str(tmp_path / "vector.fastq")
    alignment.align_reads(paths["reads"], paths["vector"], sam, threads=2)
    alignment.extract_mapped_reads(sam, vector_fastq)
    _, lq_path, _, _ = classification.classify_vector_reads(
        sam, vector_fastq, str(tmp_path), 0.9
    )

    bams = {}
    for ref in _OTHER_REFS:
        ref_sam = str(tmp_path / f"{ref}.sam")
        bams[ref] = str(tmp_path / f"alignments_hybrid_{ref}.bam")
        alignment.align_reads(lq_path, paths[ref], ref_sam, threads=2)
        alignment.sam_to_sorted_bam(ref_sam, bams[ref])

    shared = homology.SharedRegions(
        searched=frozenset(_OTHER_REFS),
        intervals={
            ref: homology.find_shared_regions(
                paths["vector"], paths[ref], preset="asm20"
            )
            for ref in _OTHER_REFS
        },
    )

    def classify(prefix, **kw):
        out = classification.classify_reads_from_bams(
            bams, lq_path, str(tmp_path / prefix), "hybrid", **kw
        )
        got = {}
        for cat, names in _categories(out).items():
            for name in names:
                got[name] = "LQ" if cat == "unmapped" else f"hybrid_{cat}"
        return got

    got = classify("hybrid", shared_regions=shared, min_unique_bases=100)
    expected = {
        r: t["expected_category"]
        for r, t in truth.items()
        if t["expected_category"] != "HQ"
    }
    wrong = sorted(
        f"{r}: {got.get(r)} != {e}" for r, e in expected.items() if got.get(r) != e
    )
    assert not wrong, f"{len(wrong)} misclassified, e.g. {wrong[:5]}"
    assert set(got) == set(expected)

    # Regression proof: without the filter, the shared-only LQ reads are hybrids
    unfiltered = classify("unfiltered")
    pure_lq = {r for r, e in expected.items() if e == "LQ"}
    assert all(unfiltered[r].startswith("hybrid_") for r in pure_lq)


# ----------------------------------------------------------------------
# CLI wiring: ref-homology, classify-best --shared-bed, run options
# ----------------------------------------------------------------------


def _classify_best_args(c, bed=None, mode="hybrid", min_unique=100):
    return argparse.Namespace(
        bams=[f"{name}={path}" for name, path in c["bams"].items()],
        source_fastq=c["fastq"],
        out_prefix=c["prefix"],
        mode=mode,
        shared_bed=bed,
        min_unique_bases=min_unique,
    )


def test_cli_classify_best_with_shared_bed(bam_case, tmp_path):
    c = bam_case
    bed = tmp_path / "shared.bed"
    homology.write_shared_bed(c["shared"], str(bed))
    assert cli.cmd_classify_best(_classify_best_args(c, str(bed))) == 0
    assert read_fastq_names(c["prefix"] + "_backbone.fastq") == {
        "spans",
        "at_threshold",
        "supp_sum",
        "minus_spans",
    }
    assert os.path.exists(c["prefix"] + "_filter_stats.json")


def test_cli_classify_best_without_shared_bed_is_unfiltered(bam_case):
    c = bam_case
    assert cli.cmd_classify_best(_classify_best_args(c)) == 0
    assert "shared_only" in read_fastq_names(c["prefix"] + "_backbone.fastq")
    assert not os.path.exists(c["prefix"] + "_filter_stats.json")


def test_cli_classify_best_shared_bed_rejected_for_impurity(bam_case, tmp_path):
    c = bam_case
    bed = tmp_path / "shared.bed"
    homology.write_shared_bed(c["shared"], str(bed))
    with pytest.raises(ValueError, match="hybrid"):
        cli.cmd_classify_best(_classify_best_args(c, str(bed), mode="impurity"))


def test_run_and_step_parsers_have_hybrid_options(monkeypatch):
    import aav_analyze

    monkeypatch.setattr(sys, "argv", ["aav_analyze.py", "run", "--run-dir", "x"])
    args = aav_analyze.parse_arguments()
    assert (args.hybrid_min_unique, args.homology_preset) == (100, "asm20")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "aav_analyze.py",
            "ref-homology",
            "--vector",
            "v.fa",
            "--ref",
            "backbone=b.fa",
            "--out-bed",
            "o.bed",
        ],
    )
    args = aav_analyze.parse_arguments()
    assert (args.preset, args.ref) == ("asm20", ["backbone=b.fa"])

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "aav_analyze.py",
            "classify-best",
            "--bams",
            "backbone=b.bam",
            "--source-fastq",
            "s.fq",
            "--out-prefix",
            "p",
        ],
    )
    args = aav_analyze.parse_arguments()
    assert (args.shared_bed, args.min_unique_bases) == (None, 100)


@pytest.mark.slow
@_needs_tools
def test_compute_shared_regions_and_cli_bed(tmp_path):
    paths = make_hybrid_testdata.generate(str(tmp_path / "data"))
    refs = {ref: paths[ref] for ref in _OTHER_REFS}
    regions = homology.compute_shared_regions(paths["vector"], refs, preset="asm20")
    assert regions.searched == frozenset(_OTHER_REFS)
    assert set(regions.intervals) == {"backbone", "host"}  # helper shares nothing

    bed = tmp_path / "shared.bed"
    args = argparse.Namespace(
        vector=paths["vector"],
        ref=[f"{name}={path}" for name, path in refs.items()],
        preset="asm20",
        threads=1,
        out_bed=str(bed),
    )
    assert cli.cmd_ref_homology(args) == 0
    assert homology.read_shared_bed(str(bed)) == regions


# ----------------------------------------------------------------------
# SUMMARIZE: shared-region filter counts in the classification stats
# ----------------------------------------------------------------------


def _fastq_dir_with_filter_stats(tmp_path, with_stats=True):
    """Classified FASTQs (2 LQ, 1 of them a backbone hybrid) plus filter stats."""
    d = tmp_path / "fastq"
    d.mkdir()
    write_fastq(d / "vector_HQ.fastq", ["hq1"])
    write_fastq(d / "vector_LQ.fastq", ["lq1", "lq2"])
    write_fastq(d / "hybrid_backbone.fastq", ["lq1"])
    for name in (
        "hybrid_helper",
        "hybrid_host",
        "impurity_backbone",
        "impurity_helper",
        "impurity_host",
    ):
        (d / f"{name}.fastq").write_text("")
    if with_stats:
        (d / "hybrid_filter_stats.json").write_text(
            json.dumps(
                {
                    "min_unique_bases": 100,
                    "rejected_alignments_by_ref": {
                        "backbone": 2,
                        "helper": 0,
                        "host": 1,
                    },
                    "rescued_to_lq_by_ref": {"backbone": 1, "helper": 0, "host": 0},
                    "rescued_to_lq": 1,
                }
            )
        )
    return str(d)


def test_summarize_includes_shared_filter_counts(tmp_path):
    stats = classification.summarize_classification(
        _fastq_dir_with_filter_stats(tmp_path)
    )
    assert stats["lq_shared_only_backbone"] == 1
    assert stats["lq_shared_only_host"] == 0
    assert stats["lq_shared_only_total"] == 1
    assert stats["shared_filter_min_unique_bases"] == 100
    # The new keys must not leak into the substring-based totals
    assert stats["total_hybrid"] == 1
    assert stats["vector_lq_pure"] == 1
    assert stats["total_reads"] == 3  # HQ 1 + pure LQ 1 + hybrid 1


def test_summarize_without_filter_stats_has_no_shared_keys(tmp_path):
    stats = classification.summarize_classification(
        _fastq_dir_with_filter_stats(tmp_path, with_stats=False)
    )
    assert not any(k.startswith(("lq_shared_only", "shared_filter")) for k in stats)


def test_summary_txt_reports_shared_only_reads(tmp_path):
    stats = classification.summarize_classification(
        _fastq_dir_with_filter_stats(tmp_path)
    )
    out = tmp_path / "summary.txt"
    cli._write_summary_txt(stats, str(out))
    assert "aligned only in vector-shared regions: 1" in out.read_text()
