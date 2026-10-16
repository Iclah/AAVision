"""
Tests for HQ/LQ vector-read classification and the alignment details it
depends on (supplementary/secondary handling, hard clips, strand).

Fast tests build small SAM/BAM files by hand. Slow tests run minimap2 on the
synthetic reads from scripts/make_hq_lq_testdata.py and compare to its truth
table.

SAM text in this file uses 1-based POS; pysam attributes are 0-based.
"""

import csv
import os
import random
import shutil
import sys

import pysam
import pytest

from aav_analyzer import alignment, classification, truncation

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
import make_hq_lq_testdata

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq[::-1].translate(_COMP)


def read_fastq(path):
    """Return {read_id: (sequence, quality)} from a 4-line FASTQ."""
    records = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return records
    with pysam.FastxFile(path) as fq:
        for entry in fq:
            records[entry.name] = (entry.sequence, entry.quality)
    return records


def write_fastq(path, records):
    with open(path, "w") as fh:
        fh.writelines(
            f"@{name}\n{seq}\n+\n{qual}\n" for name, (seq, qual) in records.items()
        )


# ----------------------------------------------------------------------
# Hand-built SAM fixture
# ----------------------------------------------------------------------

REF_LEN = 1000


@pytest.fixture
def ref_seq():
    rng = random.Random(7)
    return "".join(rng.choice("ACGT") for _ in range(REF_LEN))


@pytest.fixture
def sam_case(tmp_path, ref_seq):
    """
    A vector SAM plus the matching vector.fastq (original read orientation,
    as `samtools fastq -F 2308` produces it).

    Reads and expected classes (threshold 0.9):
      fwd_full        100M on + strand                       -> HQ
      rev_full        100M on - strand, asymmetric quality   -> HQ (written in original orientation)
      thresh          90M10S, coverage exactly 0.90          -> HQ (threshold is inclusive)
      sec             primary 40M60S + full-length secondary -> LQ
      unmapped        flag 4, not in vector.fastq            -> neither
    """
    r = ref_seq
    junk = "T" * 60
    originals = {
        "fwd_full": (r[0:100], "I" * 100),
        "rev_full": (revcomp(r[100:200]), "!" * 50 + "I" * 50),
        "thresh": (r[500:590] + junk[:10], "I" * 100),
        "sec": (r[400:440] + junk, "I" * 100),
    }

    def rec(name, flag, pos, cigar, seq, qual):
        return f"{name}\t{flag}\tvector\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tAS:i:100"

    o = originals
    lines = [
        "@HD\tVN:1.6\tSO:unsorted",
        f"@SQ\tSN:vector\tLN:{REF_LEN}",
        rec("fwd_full", 0, 1, "100M", *o["fwd_full"]),
        # Reverse-strand primary: SAM stores the reverse complement and reversed quality
        rec("rev_full", 16, 101, "100M", r[100:200], o["rev_full"][1][::-1]),
        rec("thresh", 0, 501, "90M10S", *o["thresh"]),
        rec("sec", 0, 401, "40M60S", *o["sec"]),
        # minimap2 writes secondary records with SEQ/QUAL '*'
        "sec\t256\tvector\t601\t0\t100M\t*\t0\t0\t*\t*\tAS:i:100",
        "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t" + "A" * 100 + "\t" + "I" * 100,
    ]
    sam = tmp_path / "vector.sam"
    sam.write_text("\n".join(lines) + "\n")

    vector_fastq = tmp_path / "vector.fastq"
    write_fastq(vector_fastq, originals)

    outdir = tmp_path / "out"
    outdir.mkdir()
    return {
        "sam": str(sam),
        "vector_fastq": str(vector_fastq),
        "outdir": str(outdir),
        "originals": originals,
        "expected_hq": {"fwd_full", "rev_full", "thresh"},
        "expected_lq": {"sec"},
    }


def _classify(case, min_coverage=0.9):
    hq_path, lq_path, hq_n, lq_n = classification.classify_vector_reads(
        case["sam"], case["vector_fastq"], case["outdir"], min_coverage
    )
    return read_fastq(hq_path), read_fastq(lq_path), hq_n, lq_n


# ----------------------------------------------------------------------
# read_interval: the part of the ORIGINAL read one record covers
# (0-based half-open, read coordinates in FASTQ orientation)
# ----------------------------------------------------------------------


def _segment(cigar, seq_len, flag=0):
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "vector", "LN": REF_LEN}]})
    seg = pysam.AlignedSegment(header)
    seg.query_name = "q"
    seg.flag = flag
    seg.reference_id = 0
    seg.reference_start = 0
    seg.cigarstring = cigar
    seg.query_sequence = "A" * seq_len
    seg.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
    return seg


def test_read_interval_soft_clip():
    assert alignment.read_interval(_segment("90M10S", 100)) == (0, 90)


def test_read_interval_hard_clips_count_in_read_coordinates():
    # 50 aligned bases of a 100 bp read; the hard-clipped half is not in SEQ
    assert alignment.read_interval(_segment("50H50M", 50, flag=2048)) == (50, 100)


def test_read_interval_mixed_soft_and_hard_clips():
    assert alignment.read_interval(_segment("20H60M20S", 80, flag=2048)) == (20, 80)


def test_read_interval_ignores_reference_deletions():
    # Expected behaviour (docs/DECISIONS.md): HQ measures read-side clipping,
    # so an internal deletion does not shrink the covered read interval.
    assert alignment.read_interval(_segment("50M500D50M", 100)) == (0, 100)


def test_read_interval_insertions_count():
    assert alignment.read_interval(_segment("10S40M20I30M", 100)) == (10, 100)


def test_read_interval_reverse_strand_is_flipped():
    # SAM stores a minus-strand record reverse-complemented: 30S70M covers
    # record positions [30,100), i.e. original read positions [0,70)
    assert alignment.read_interval(_segment("30S70M", 100, flag=16)) == (0, 70)
    assert alignment.read_interval(_segment("30H70M", 70, flag=2064)) == (0, 70)


# ----------------------------------------------------------------------
# is_snapback: (contig, ref_start, ref_end, is_reverse) per record,
# reference spans 0-based half-open
# ----------------------------------------------------------------------


def test_snapback_opposite_strands_overlapping():
    spans = [("vector", 0, 2000, False), ("vector", 1800, 3000, True)]
    assert alignment.is_snapback(spans, min_overlap=200)


def test_snapback_overlap_threshold_inclusive():
    assert alignment.is_snapback(
        [("vector", 0, 1000, False), ("vector", 800, 1500, True)], min_overlap=200
    )
    assert not alignment.is_snapback(
        [("vector", 0, 1000, False), ("vector", 801, 1500, True)], min_overlap=200
    )


def test_snapback_needs_opposite_strands():
    spans = [("vector", 0, 2000, False), ("vector", 0, 2000, False)]
    assert not alignment.is_snapback(spans, min_overlap=200)


def test_snapback_needs_same_contig():
    spans = [("vector", 0, 2000, False), ("other", 0, 2000, True)]
    assert not alignment.is_snapback(spans, min_overlap=200)


def test_snapback_single_record():
    assert not alignment.is_snapback([("vector", 0, 2000, True)], min_overlap=200)


# ----------------------------------------------------------------------
# classify_vector_reads: merged primary + supplementary coverage
# ----------------------------------------------------------------------

STRUCT_REF_LEN = 10000

# name -> (read length, records, expected class)
# record: (read_start, read_end, ref_start, strand, is_supplementary)
# read_start/read_end are positions in the ORIGINAL read (0-based half-open);
# ref_start is 0-based. The first record is the primary unless reordered.
STRUCTURES = {
    # 3,000 vector + 1,000 human, one alignment: 0.75
    "single_75": (4000, [(0, 3000, 0, "+", False)], "LQ"),
    # 2,000 vector + 1,000 human + 1,000 vector: 0.75 merged
    "split_75": (4000, [(0, 2000, 0, "+", False), (3000, 4000, 3000, "+", True)], "LQ"),
    # 2,000 vector forward + the same 2,000 reversed: coverage 1.0 but snapback
    "snapback": (4000, [(0, 2000, 0, "+", False), (2000, 4000, 0, "-", True)], "LQ"),
    # rcAAV-like: ITR + 4,010 foreign + ITR: 290 / 4,300
    "rcaav": (4300, [(0, 145, 0, "+", False), (4155, 4300, 4555, "-", True)], "LQ"),
    # single alignment, 0.95
    "single_95": (4000, [(0, 3800, 0, "+", False)], "HQ"),
    # all vector, split over distant non-overlapping regions: promoted by supplementary
    "split_all_vector": (
        4000,
        [(0, 2500, 2000, "+", False), (2500, 4000, 0, "+", True)],
        "HQ",
    ),
    # read pieces [0,2000) and [1000,3000) overlap: 3,000 merged = 0.75 (sum would be 1.0)
    "overlap_counted_once": (
        4000,
        [(0, 2000, 0, "+", False), (1000, 3000, 5000, "+", True)],
        "LQ",
    ),
    # full genome + 145 bp ITR supplementary on the other strand: overlap < 200
    "full_itr_supp": (4000, [(0, 4000, 0, "+", False), (0, 145, 0, "-", True)], "HQ"),
    # foldback with exactly 200 bp reference overlap: snapback (inclusive)
    "snap_overlap_200": (
        4000,
        [(0, 3800, 0, "+", False), (3800, 4000, 3600, "-", True)],
        "LQ",
    ),
    # 199 bp reference overlap: not a snapback, coverage 1.0
    "snap_overlap_199": (
        4000,
        [(0, 3800, 0, "+", False), (3800, 4000, 3601, "-", True)],
        "HQ",
    ),
    # minus-strand primary covering the whole read
    "rev_full": (4000, [(0, 4000, 1000, "-", False)], "HQ"),
}
# Same structures with the supplementary record written before the primary
REORDERED = {"split_all_vector", "snapback"}


def _structure_case(tmp_path, clip):
    """
    Build vector.sam and vector.fastq for STRUCTURES.

    clip: 'H' writes supplementary records hard-clipped (minimap2 default),
    'S' soft-clipped (minimap2 -Y). Primary records are always soft-clipped.
    """
    rng = random.Random(11)
    originals, lines, expected = {}, [], {}
    for name, (length, records, label) in STRUCTURES.items():
        for suffix in ("", "_suppfirst") if name in REORDERED else ("",):
            read_name = name + suffix
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            qual = "".join(rng.choice("5?I") for _ in range(length))
            originals[read_name] = (seq, qual)
            expected[read_name] = label

            recs = []
            for r_start, r_end, ref_start, strand, is_supp in records:
                rev = strand == "-"
                rec_seq = revcomp(seq) if rev else seq
                rec_qual = qual[::-1] if rev else qual
                # Interval in record (SAM) orientation
                s0, e0 = (length - r_end, length - r_start) if rev else (r_start, r_end)
                op = clip if is_supp else "S"
                cigar = (f"{s0}{op}" if s0 else "") + f"{e0 - s0}M"
                cigar += f"{length - e0}{op}" if e0 < length else ""
                if op == "H":
                    rec_seq, rec_qual = rec_seq[s0:e0], rec_qual[s0:e0]
                flag = (16 if rev else 0) | (2048 if is_supp else 0)
                recs.append(
                    f"{read_name}\t{flag}\tvector\t{ref_start + 1}\t60\t{cigar}"
                    f"\t*\t0\t0\t{rec_seq}\t{rec_qual}\tAS:i:100"
                )
            if suffix:
                recs.reverse()
            lines.extend(recs)

    sam = tmp_path / f"vector_{clip}.sam"
    sam.write_text(
        f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:vector\tLN:{STRUCT_REF_LEN}\n"
        + "\n".join(lines)
        + "\n"
    )
    vector_fastq = tmp_path / f"vector_{clip}.fastq"
    write_fastq(vector_fastq, originals)
    outdir = tmp_path / f"out_{clip}"
    outdir.mkdir()
    return {
        "sam": str(sam),
        "vector_fastq": str(vector_fastq),
        "outdir": str(outdir),
        "originals": originals,
        "expected": expected,
    }


@pytest.mark.parametrize("clip", ["H", "S"])
def test_structures_classified(tmp_path, clip):
    case = _structure_case(tmp_path, clip)
    hq, lq, hq_n, lq_n = _classify(case)
    got = {**{n: "HQ" for n in hq}, **{n: "LQ" for n in lq}}
    wrong = {n: (got.get(n), e) for n, e in case["expected"].items() if got.get(n) != e}
    assert not wrong, f"(got, expected): {wrong}"
    assert (hq_n, lq_n) == (len(hq), len(lq))


@pytest.mark.parametrize("clip", ["H", "S"])
def test_structures_written_as_full_original_reads(tmp_path, clip):
    case = _structure_case(tmp_path, clip)
    hq, lq, _, _ = _classify(case)
    for name, record in {**hq, **lq}.items():
        assert record == case["originals"][name], name


def test_hard_and_soft_clipped_supplementary_agree(tmp_path):
    hard, soft = _structure_case(tmp_path, "H"), _structure_case(tmp_path, "S")
    hq_h, lq_h, _, _ = _classify(hard)
    hq_s, lq_s, _, _ = _classify(soft)
    assert set(hq_h) == set(hq_s)
    assert set(lq_h) == set(lq_s)


def test_snapback_overlap_is_a_parameter(tmp_path):
    case = _structure_case(tmp_path, "H")
    hq_path, _, _, _ = classification.classify_vector_reads(
        case["sam"], case["vector_fastq"], case["outdir"], 0.9, min_snapback_overlap=100
    )
    hq = read_fastq(hq_path)
    # full_itr_supp has a 145 bp opposite-strand overlap: a snapback at 100 bp
    assert "full_itr_supp" not in hq
    assert "snap_overlap_199" not in hq


def _write_sam(path, records):
    path.write_text(
        f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:vector\tLN:{REF_LEN}\n"
        + "\n".join(records)
        + "\n"
    )


def test_inconsistent_read_length_fails(tmp_path):
    sam = tmp_path / "vector.sam"
    _write_sam(
        sam,
        [
            "q\t0\tvector\t1\t60\t100M\t*\t0\t0\t" + "A" * 100 + "\t" + "I" * 100,
            # 50M60H implies a 110 bp read
            "q\t2048\tvector\t201\t60\t50M60H\t*\t0\t0\t" + "A" * 50 + "\t" + "I" * 50,
        ],
    )
    fq = tmp_path / "vector.fastq"
    write_fastq(fq, {"q": ("A" * 100, "I" * 100)})
    with pytest.raises(ValueError, match="read length"):
        alignment.filter_reads_by_coverage(
            str(sam), str(fq), str(tmp_path / "hq.fastq")
        )


def test_supplementary_without_primary_fails(tmp_path):
    sam = tmp_path / "vector.sam"
    _write_sam(
        sam,
        ["q\t2048\tvector\t1\t60\t50M50H\t*\t0\t0\t" + "A" * 50 + "\t" + "I" * 50],
    )
    fq = tmp_path / "vector.fastq"
    fq.write_text("")
    with pytest.raises(ValueError, match="primary"):
        alignment.filter_reads_by_coverage(
            str(sam), str(fq), str(tmp_path / "hq.fastq")
        )


def test_secondary_alignment_is_ignored(sam_case):
    hq, lq, _, _ = _classify(sam_case)
    assert "sec" not in hq
    assert "sec" in lq


def test_threshold_is_inclusive(sam_case):
    hq, _, _, _ = _classify(sam_case)
    assert "thresh" in hq


def test_hq_and_lq_partition_mapped_reads(sam_case):
    hq, lq, hq_n, lq_n = _classify(sam_case)
    assert set(hq) == sam_case["expected_hq"]
    assert set(lq) == sam_case["expected_lq"]
    assert set(hq).isdisjoint(lq)
    assert set(hq) | set(lq) == set(sam_case["originals"])
    assert (hq_n, lq_n) == (len(hq), len(lq))


def test_hq_reads_written_as_original_reads(sam_case):
    # Minus-strand reads must not be reverse-complemented, and reads must be
    # the full original read, not one alignment record's portion of it.
    hq, lq, _, _ = _classify(sam_case)
    for name, record in {**hq, **lq}.items():
        assert record == sam_case["originals"][name], name


def test_empty_input(tmp_path):
    sam = tmp_path / "vector.sam"
    sam.write_text(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:vector\tLN:{REF_LEN}\n")
    vector_fastq = tmp_path / "vector.fastq"
    vector_fastq.write_text("")
    hq_path, lq_path, hq_n, lq_n = classification.classify_vector_reads(
        str(sam), str(vector_fastq), str(tmp_path), 0.9
    )
    assert (hq_n, lq_n) == (0, 0)
    assert os.path.exists(hq_path) and os.path.exists(lq_path)


# ----------------------------------------------------------------------
# align_reads
# ----------------------------------------------------------------------


def test_align_reads_soft_clips_supplementary(tmp_path, monkeypatch):
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)

    monkeypatch.setattr(alignment.subprocess, "run", fake_run)
    alignment.align_reads("reads.fq", "ref.fa", str(tmp_path / "out.sam"), threads=2)
    assert calls and calls[0][0] == "minimap2"
    assert "-Y" in calls[0]


# ----------------------------------------------------------------------
# compute_truncation_hotspots
# ----------------------------------------------------------------------


def test_truncation_uses_primary_records_only(tmp_path):
    """
    One 900 bp read: primary 600M300S at ref [0, 600), supplementary
    600H300M at ref [700, 1000). Plus a secondary at ref [300, 400).
    Only the primary's termini (0-based 0 and 599) and coverage may count.
    """
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "vector", "LN": REF_LEN}]})
    bam_path = str(tmp_path / "vector.bam")

    def seg(flag, start, cigar, seq_len):
        s = pysam.AlignedSegment(header)
        s.query_name = "chimera"
        s.flag = flag
        s.reference_id = 0
        s.reference_start = start
        s.mapping_quality = 60
        s.cigarstring = cigar
        if seq_len:
            s.query_sequence = "A" * seq_len
            s.query_qualities = pysam.qualitystring_to_array("I" * seq_len)
        return s

    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        out.write(seg(0, 0, "600M300S", 900))
        out.write(seg(256, 300, "100M", 0))
        out.write(seg(2048, 700, "600H300M", 300))
    pysam.index(bam_path)

    df = truncation.compute_truncation_hotspots(bam_path)
    termini = dict(zip(df["pos"], df["termini"]))  # pos is 1-based
    cov = dict(zip(df["pos"], df["coverage"]))

    assert termini[1] == 1 and termini[600] == 1
    assert termini[301] == 0 and termini[400] == 0  # secondary
    assert termini[701] == 0 and termini[1000] == 0  # supplementary
    assert cov[350] == 1
    assert cov[800] == 0


# ----------------------------------------------------------------------
# End to end with minimap2 on synthetic reads
# ----------------------------------------------------------------------

_TOOLS = ("minimap2", "samtools", "seqkit")


@pytest.mark.slow
@pytest.mark.skipif(
    any(shutil.which(t) is None for t in _TOOLS),
    reason=f"needs {', '.join(_TOOLS)} on PATH",
)
def test_end_to_end_matches_truth_table(tmp_path):
    paths = make_hq_lq_testdata.generate(str(tmp_path / "data"))
    with open(paths["truth"]) as fh:
        truth = {row["read_id"]: row for row in csv.DictReader(fh, delimiter="\t")}
    inputs = read_fastq(paths["reads"])

    sam = str(tmp_path / "vector.sam")
    vector_fastq = str(tmp_path / "vector.fastq")
    alignment.align_reads(paths["reads"], paths["vector"], sam, threads=2)
    alignment.extract_mapped_reads(sam, vector_fastq)
    hq_path, lq_path, hq_n, lq_n = classification.classify_vector_reads(
        sam, vector_fastq, str(tmp_path), 0.9
    )
    hq, lq = read_fastq(hq_path), read_fastq(lq_path)

    expected = {
        c: {r for r, t in truth.items() if t["expected_class"] == c}
        for c in ("HQ", "LQ")
    }
    wrong_hq = sorted(set(hq) - expected["HQ"])
    wrong_lq = sorted(set(lq) - expected["LQ"])
    assert not wrong_hq, f"{len(wrong_hq)} reads wrongly HQ, e.g. {wrong_hq[:5]}"
    assert not wrong_lq, f"{len(wrong_lq)} reads wrongly LQ, e.g. {wrong_lq[:5]}"
    assert set(hq) == expected["HQ"]
    assert set(lq) == expected["LQ"]
    assert (hq_n, lq_n) == (len(expected["HQ"]), len(expected["LQ"]))

    # Output reads are the input reads, byte for byte (orientation, length, quality)
    for name, record in {**hq, **lq}.items():
        assert record == inputs[name], name
