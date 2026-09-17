"""
Fine-grained subcommand handlers for the AAV analysis pipeline.

Each handler is a thin wrapper around the core module functions
(alignment, classification, truncation, taxonomy, report). They take
explicit input/output file paths so that each pipeline step can run
independently — this is what makes the Nextflow port a near 1:1 mapping.

The handlers here do NOT do sample discovery, merging, or reference
resolution; that orchestration lives in aav_analyze.py's `run` command
(and will be replicated by the Nextflow workflow layer). These handlers
operate on already-resolved single files.

Each `cmd_*` function takes an argparse.Namespace and returns an int
exit code (0 = success).
"""

import os
import json
import logging

from . import alignment
from . import classification
from . import truncation
from . import taxonomy
from . import report
from . import utils


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def _parse_ref_pairs(pairs):
    """
    Parse a list of 'name=path' or 'name:path' strings into a dict.

    Used for subcommands that accept multiple references or multiple
    named FASTQ files on the command line.
    """
    result = {}
    if not pairs:
        return result
    for item in pairs:
        if '=' in item:
            name, path = item.split('=', 1)
        elif ':' in item:
            name, path = item.split(':', 1)
        else:
            raise ValueError(
                f"Expected 'name=path' (or 'name:path'), got: {item}"
            )
        result[name.strip()] = path.strip()
    return result


def _ensure_parent_dir(filepath):
    """Create the parent directory of a file path if it doesn't exist."""
    parent = os.path.dirname(os.path.abspath(filepath))
    if parent:
        os.makedirs(parent, exist_ok=True)


def _same_file(a, b):
    """
    True if two paths refer to the same file on disk.

    Uses os.path.samefile when both exist (handles symlinks, ./ prefixes,
    relative-vs-absolute), and falls back to normalized-path comparison
    when one side doesn't exist yet.
    """
    try:
        return os.path.exists(a) and os.path.exists(b) and os.path.samefile(a, b)
    except OSError:
        return os.path.abspath(a) == os.path.abspath(b)


# ----------------------------------------------------------------------
# Subcommand: prepare-fastq
# ----------------------------------------------------------------------

def cmd_prepare_fastq(args):
    """
    Merge a sample's FASTQ files (and optionally length-filter them).

    Wraps utils.prepare_sample_fastq, which handles single-file,
    multi-file merge, and the optional length filter. Prints the path
    of the prepared FASTQ to stdout so a caller (e.g. Nextflow) can
    capture it.
    """
    result = utils.prepare_sample_fastq(
        args.sample_dir,
        args.sample_name,
        min_length=args.min_length,
        work_dir=getattr(args, 'work_dir', None)
    )
    if result is None:
        logging.error("No FASTQ files found in sample directory")
        return 1

    # If --out is given, copy the prepared FASTQ to that explicit path.
    # This is what workflow managers (Nextflow) use so the output lands in
    # the task's work directory rather than inside the (possibly read-only,
    # symlink-staged) sample directory. Also copy the length-filter stats
    # alongside if they were produced.
    if getattr(args, 'out', None):
        import shutil as _shutil
        _ensure_parent_dir(args.out)
        # copyfile (not copy) — copies bytes only, no chmod. chmod fails on
        # Windows-mounted filesystems (e.g. /mnt/c under WSL) where the
        # container runs as the host UID, and preserving the source mode for
        # a freshly-produced output is unnecessary.
        # Skip the copy when source and destination are already the same file
        # (happens when --work-dir places the output exactly where --out wants
        # it, e.g. under Nextflow with --work-dir . and a relative --out).
        if not _same_file(result, args.out):
            _shutil.copyfile(result, args.out)

        if args.out_filter_stats:
            # Stats live in the work dir (or sample dir if no work dir)
            search_dir = getattr(args, 'work_dir', None) or args.sample_dir
            stats_src = os.path.join(search_dir, 'length_filter_stats.json')
            if os.path.exists(stats_src):
                _ensure_parent_dir(args.out_filter_stats)
                if not _same_file(stats_src, args.out_filter_stats):
                    _shutil.copyfile(stats_src, args.out_filter_stats)

        print(os.path.abspath(args.out))
        return 0

    # Otherwise emit the resolved in-sample-dir path (standalone behavior)
    print(result)
    return 0


# ----------------------------------------------------------------------
# Subcommand: align-vector
# ----------------------------------------------------------------------

def cmd_align_vector(args):
    """Align reads to the vector reference (minimap2 -> SAM)."""
    _ensure_parent_dir(args.out_sam)
    alignment.align_reads(args.reads, args.vector, args.out_sam, args.threads)
    return 0


# ----------------------------------------------------------------------
# Subcommand: extract-vector
# ----------------------------------------------------------------------

def cmd_extract_vector(args):
    """
    From a vector SAM:
      - extract primary mapped reads to FASTQ
      - build a sorted, indexed vector BAM (primary + supplementary)
      - extract the unmapped (non-vector) reads to a separate FASTQ
    """
    _ensure_parent_dir(args.out_fastq)
    _ensure_parent_dir(args.out_bam)
    _ensure_parent_dir(args.out_unmap)

    alignment.extract_mapped_reads(args.sam, args.out_fastq)
    alignment.create_bam_for_mapped_reads(args.sam, args.out_bam)
    # Unmapped = all input reads minus the vector-mapped reads
    alignment.create_fastq_difference(args.reads, args.out_fastq, args.out_unmap)
    return 0


# ----------------------------------------------------------------------
# Subcommand: classify-hq-lq
# ----------------------------------------------------------------------

def cmd_classify_hq_lq(args):
    """
    Split vector reads into high-quality / low-quality by alignment
    coverage. Writes vector_HQ.fastq and vector_LQ.fastq.
    """
    _ensure_parent_dir(args.out_hq)
    _ensure_parent_dir(args.out_lq)

    # filter_reads_by_coverage writes HQ; create_fastq_difference makes LQ
    hq_count = alignment.filter_reads_by_coverage(
        args.sam, args.out_hq, args.min_coverage
    )
    lq_count = alignment.create_fastq_difference(
        args.vector_fastq, args.out_hq, args.out_lq
    )
    logging.info(f"HQ reads: {hq_count}, LQ reads: {lq_count}")
    return 0


# ----------------------------------------------------------------------
# Subcommand: align-ref  (generic; used for hybrid and impurity fan-out)
# ----------------------------------------------------------------------

def cmd_align_ref(args):
    """
    Align a set of reads to a single reference and produce a sorted BAM.

    This is the generic per-reference aligner used for both the hybrid
    fan-out (LQ reads vs each impurity ref) and the impurity fan-out
    (unmapped reads vs each impurity ref). The caller decides which
    reads and which ref; this command just does one alignment.
    """
    _ensure_parent_dir(args.out_bam)

    # Align to SAM (temp), convert to sorted+indexed BAM, clean up SAM
    sam_path = args.out_bam.rsplit('.bam', 1)[0] + '.sam'
    alignment.align_reads(args.reads, args.ref, sam_path, args.threads)
    alignment.sam_to_sorted_bam(sam_path, args.out_bam)
    if os.path.exists(sam_path):
        os.remove(sam_path)
    return 0


# ----------------------------------------------------------------------
# Subcommand: classify-best  (generic; hybrid or impurity)
# ----------------------------------------------------------------------

def cmd_classify_best(args):
    """
    Assign each read to its best-matching reference across multiple BAMs.

    Wraps classification.classify_reads_from_bams. The --bams arguments
    are 'refname=path' pairs. Output FASTQs are written as
    {out_prefix}_{refname}.fastq plus {out_prefix}_unmapped.fastq.
    """
    bam_files = _parse_ref_pairs(args.bams)
    if not bam_files:
        logging.error("No BAM files provided to classify-best")
        return 1

    _ensure_parent_dir(args.out_prefix + "_x")  # ensure prefix parent exists

    classification.classify_reads_from_bams(
        bam_files,
        args.source_fastq,
        args.out_prefix,
        args.mode
    )
    return 0


# ----------------------------------------------------------------------
# Subcommand: compute-coverage
# ----------------------------------------------------------------------

def cmd_compute_coverage(args):
    """
    Compute per-position coverage depth for one BAM (one category) and
    write it as JSON: {position: depth}.
    """
    _ensure_parent_dir(args.out_json)
    coverage = truncation.compute_coverage(args.bam, max_pos=args.max_pos)
    with open(args.out_json, 'w') as f:
        json.dump(coverage, f, indent=2)
    logging.info(f"Coverage for '{args.category}' written to {args.out_json}")
    return 0


# ----------------------------------------------------------------------
# Subcommand: compute-truncation
# ----------------------------------------------------------------------

def cmd_compute_truncation(args):
    """
    Compute truncation hotspots for the vector (rAAV) BAM. Writes both
    a full CSV and a position->enrichment JSON for the report.
    """
    _ensure_parent_dir(args.out_csv)
    _ensure_parent_dir(args.out_json)

    trunc_df = truncation.compute_truncation_hotspots(args.bam)
    trunc_df.to_csv(args.out_csv, index=False)

    trunc_data = truncation.get_truncation_data_dict(trunc_df)
    with open(args.out_json, 'w') as f:
        json.dump(trunc_data, f, indent=2)
    logging.info(f"Truncation data written to {args.out_csv} and {args.out_json}")
    return 0


# ----------------------------------------------------------------------
# Subcommand: compute-lengths
# ----------------------------------------------------------------------

def cmd_compute_lengths(args):
    """
    Compute read-length distributions for a set of named FASTQ files.

    --fastqs arguments are 'name=path' pairs (e.g. vector=vector.fastq
    helper=impurity_helper.fastq ...). Writes one JSON mapping each
    name to its binned length distribution.
    """
    fastq_map = _parse_ref_pairs(args.fastqs)
    if not fastq_map:
        logging.error("No FASTQ files provided to compute-lengths")
        return 1

    # truncation.analyze_read_lengths expects a list of (name, path) tuples
    fastq_list = [(name, path) for name, path in fastq_map.items()]

    _ensure_parent_dir(args.out_json)
    distributions = truncation.analyze_read_lengths(fastq_list)
    with open(args.out_json, 'w') as f:
        json.dump(distributions, f, indent=2)
    logging.info(f"Length distributions written to {args.out_json}")
    return 0


# ----------------------------------------------------------------------
# Subcommand: summarize
# ----------------------------------------------------------------------

def cmd_summarize(args):
    """
    Aggregate classification statistics across all categorized FASTQ
    files in a directory. Writes both a JSON (classification_stats.json)
    and a human-readable summary.txt.
    """
    optional_refs = args.optional_refs or []
    stats = classification.summarize_classification(args.fastq_dir, optional_refs)

    _ensure_parent_dir(args.out_json)
    with open(args.out_json, 'w') as f:
        json.dump(stats, f, indent=2)

    if args.out_txt:
        _ensure_parent_dir(args.out_txt)
        _write_summary_txt(stats, args.out_txt)

    logging.info(f"Statistics written to {args.out_json}")
    return 0


def _write_summary_txt(stats, out_txt):
    """Write the human-readable summary.txt (mirrors the run orchestrator)."""
    with open(out_txt, 'w') as f:
        f.write("AAV Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total reads: {stats['total_reads']:,}\n")
        f.write(f"  Vector reads: {stats['total_vector']:,} "
                f"({stats['vector_percent']:.2f}%)\n")
        f.write(f"    - High quality: {stats['vector_hq']:,}\n")
        f.write(f"    - Low quality: {stats.get('vector_lq_pure', stats.get('vector_lq', 0)):,}\n")
        f.write(f"  Hybrid reads: {stats['total_hybrid']:,} "
                f"({stats['hybrid_percent']:.2f}%)\n")
        for key in ['backbone', 'helper', 'host', 'rep_cap']:
            hybrid_key = f'{key}_hybrid'
            if hybrid_key in stats:
                f.write(f"    - {key.capitalize()}: {stats[hybrid_key]:,}\n")
        f.write(f"  Impurity reads: {stats['total_impurity']:,} "
                f"({stats['impurity_percent']:.2f}%)\n")
        for key in ['backbone', 'helper', 'host', 'rep_cap', 'spike_in']:
            impurity_key = f'{key}_impurity'
            if impurity_key in stats:
                f.write(f"    - {key.capitalize()}: {stats[impurity_key]:,}\n")


# ----------------------------------------------------------------------
# Subcommand: taxonomy  (combined kraken + bracken + summary)
# ----------------------------------------------------------------------

def cmd_taxonomy(args):
    """
    Run the full taxonomy unit (Kraken2 + optional Bracken) and write the
    taxonomy_summary.json that the report consumes. This mirrors the
    standalone's run_taxonomy_analysis end-to-end so the Nextflow output
    matches the control. Kraken2 needs the whole DB in RAM, so this step is
    intentionally heavyweight and opt-in (gated on --db at the workflow level).

    Fails LOUDLY if kraken2 is missing or the DB is invalid — taxonomy was
    explicitly requested, so silently skipping would hide a real problem.
    Bracken remains a graceful internal skip if the bracken tool/DB files
    are absent (mirrors standalone behavior).
    """
    if not taxonomy.check_kraken_available():
        logging.error("kraken2 not found on PATH but taxonomy was requested")
        return 1
    if not taxonomy.validate_kraken_db(args.db):
        logging.error(f"Invalid Kraken2 database: {args.db}")
        return 1

    os.makedirs(args.out_dir, exist_ok=True)

    taxonomy_data = taxonomy.run_taxonomy_analysis(
        fastq_file=args.reads,
        db_path=args.db,
        output_dir=args.out_dir,
        threads=args.threads,
        confidence=args.confidence,
        bracken_level=args.level,
        bracken_threshold=args.threshold,
        bracken_read_length=args.read_length,
    )

    _ensure_parent_dir(args.out_summary)
    with open(args.out_summary, 'w') as f:
        json.dump(taxonomy_data, f, indent=2)
    logging.info(f"Taxonomy summary written to {args.out_summary}")
    return 0


# ----------------------------------------------------------------------
# Subcommand: taxonomy-kraken
# ----------------------------------------------------------------------

def cmd_taxonomy_kraken(args):
    """Run Kraken2 classification on a FASTQ of reads."""
    if not taxonomy.check_kraken_available():
        logging.error("kraken2 not found on PATH")
        return 1
    if not taxonomy.validate_kraken_db(args.db):
        logging.error(f"Invalid Kraken2 database: {args.db}")
        return 1

    os.makedirs(args.out_dir, exist_ok=True)
    taxonomy.run_kraken2(
        args.reads, args.db, args.out_dir,
        threads=args.threads, confidence=args.confidence
    )
    return 0


# ----------------------------------------------------------------------
# Subcommand: taxonomy-bracken
# ----------------------------------------------------------------------

def cmd_taxonomy_bracken(args):
    """Run Bracken abundance re-estimation on a Kraken2 report."""
    os.makedirs(args.out_dir, exist_ok=True)
    result = taxonomy.run_bracken(
        args.kraken_report, args.db, args.out_dir,
        level=args.level, threshold=args.threshold
    )
    if result is None:
        logging.warning("Bracken did not produce output (may be unavailable)")
        return 0  # non-fatal, mirrors pipeline behavior
    return 0


# ----------------------------------------------------------------------
# Subcommand: generate-report
# ----------------------------------------------------------------------

def cmd_generate_report(args):
    """
    Assemble the final HTML report from previously-computed JSON inputs.

    Inputs are paths to JSON files produced by earlier steps:
      --stats        classification_stats.json
      --truncation   truncation JSON (position -> enrichment)
      --coverage     coverage JSON (category -> {pos: depth})
      --lengths      read length distributions JSON
      --taxonomy     (optional) taxonomy_summary.json
    """
    with open(args.stats) as f:
        stats = json.load(f)
    with open(args.truncation) as f:
        truncation_data = json.load(f)
    with open(args.coverage) as f:
        coverage_data = json.load(f)
    with open(args.lengths) as f:
        length_distributions = json.load(f)

    taxonomy_data = None
    if args.taxonomy and os.path.exists(args.taxonomy):
        with open(args.taxonomy) as f:
            taxonomy_data = json.load(f)

    _ensure_parent_dir(args.out)
    report.generate_html_report(
        stats,
        truncation_data,
        length_distributions,
        args.out,
        coverage_data=coverage_data,
        taxonomy_data=taxonomy_data
    )
    return 0


# ----------------------------------------------------------------------
# Registry: maps subcommand name -> handler function
# ----------------------------------------------------------------------

HANDLERS = {
    'prepare-fastq': cmd_prepare_fastq,
    'align-vector': cmd_align_vector,
    'extract-vector': cmd_extract_vector,
    'classify-hq-lq': cmd_classify_hq_lq,
    'align-ref': cmd_align_ref,
    'classify-best': cmd_classify_best,
    'compute-coverage': cmd_compute_coverage,
    'compute-truncation': cmd_compute_truncation,
    'compute-lengths': cmd_compute_lengths,
    'summarize': cmd_summarize,
    'taxonomy': cmd_taxonomy,
    'taxonomy-kraken': cmd_taxonomy_kraken,
    'taxonomy-bracken': cmd_taxonomy_bracken,
    'generate-report': cmd_generate_report,
}
