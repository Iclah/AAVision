#!/usr/bin/env python3
"""
AAVision — AAV analysis pipeline

Main command-line tool for analyzing AAV sequencing data from Nanopore runs.
Performs alignment, quality filtering, read classification, and generates
interactive HTML reports with visualization.

Usage:
    python aav_analyze.py --run-dir my_run -o output [-t 8] [-q 90] [-v]
                          [--kraken-db PATH] [--force]

See `aav_analyze.py --help` for the full list of options and the expected
run directory layout.
"""

import os
import sys
import time
import json
import shutil
import argparse
import logging
from pathlib import Path
from typing import List, Dict, Optional

# Import modules (assuming they're in aav_analyzer package)
try:
    from aav_analyzer import utils, alignment, classification, truncation, report, taxonomy
except ImportError:
    # Fallback for standalone script
    import utils
    import alignment
    import classification
    import truncation
    import report
    import taxonomy


def _add_run_subparser(subparsers):
    """Add the 'run' subcommand — the full multi-sample pipeline orchestrator."""
    p = subparsers.add_parser(
        'run',
        help='Run the full pipeline on all samples in a run directory',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze all samples in a run directory
  python aav_analyze.py run --run-dir my_run -o output

  # With Kraken2 taxonomic classification
  python aav_analyze.py run --run-dir my_run -o output --kraken-db /path/to/k2_standard

  # Verbose, 8 threads, custom coverage threshold, force overwrite, length filter
  python aav_analyze.py run --run-dir my_run -o output -t 8 -q 85 -v --force --min-length 500

Expected run directory layout:
  my_run/
  ├── references/                  # run-level references (shared across samples)
  │   ├── vector.fasta             # required
  │   ├── backbone.fasta           # required
  │   ├── helper.fasta             # required
  │   ├── host.fasta               # required
  │   ├── rep_cap.fasta            # optional
  │   └── spike_in.fasta           # optional
  └── samples/
      ├── sample_A/
      │   ├── fastq_pass/          # FASTQ files (single dir or one level deep)
      │   │   ├── chunk_0.fastq.gz #   tool finds them automatically and merges
      │   │   └── chunk_1.fastq.gz #   originals are preserved
      │   └── references/          # OPTIONAL — overrides run-level refs
      │       └── vector.fasta     #   only files that differ need to be present
      └── sample_B/
          └── reads.fastq.gz       # single pre-merged file is also fine
        """
    )
    p.add_argument('--run-dir', required=True,
                   help='Path to the run directory (must contain samples/; '
                        'optionally contains references/ with shared refs)')
    p.add_argument('-q', '--quality', type=float, default=90.0,
                   help='Minimum alignment coverage for HQ reads (default: 90%%)')
    p.add_argument('--min-length', type=int, default=0,
                   help='Minimum read length filter in bp (default: 0 = off). '
                        'Filtered reads go to merged.filtered.fastq.gz; '
                        'originals preserved.')
    p.add_argument('-o', '--output', default='analysis',
                   help='Output directory (default: analysis). One subdir per sample.')
    p.add_argument('-t', '--threads', type=int, default=4,
                   help='Number of threads for alignment (default: 4)')
    p.add_argument('-v', '--verbose', action='store_true',
                   help='Enable verbose logging')
    p.add_argument('--force', action='store_true',
                   help='Overwrite the output directory if it already exists')
    p.add_argument('--kraken-db', default=None,
                   help='Path to Kraken2 database. If provided, runs Kraken2/'
                        'Bracken on unmapped impurity reads.')
    p.add_argument('--kraken-confidence', type=float, default=0.05,
                   help='Kraken2 confidence threshold (default: 0.05, ONT-tuned)')
    p.add_argument('--bracken-level', default='S',
                   choices=['D', 'P', 'C', 'O', 'F', 'G', 'S'],
                   help='Bracken taxonomic level (default: S)')
    p.add_argument('--bracken-threshold', type=int, default=10,
                   help='Bracken minimum read count threshold (default: 10)')


def _add_step_subparsers(subparsers):
    """Add the fine-grained per-step subcommands (one per Nextflow process)."""

    # prepare-fastq
    p = subparsers.add_parser('prepare-fastq',
                              help='Merge (and optionally length-filter) a sample\'s FASTQ files')
    p.add_argument('--sample-dir', required=True)
    p.add_argument('--sample-name', required=True)
    p.add_argument('--min-length', type=int, default=0)
    p.add_argument('--work-dir', default=None,
                   help='Directory for merged/filtered intermediates (default: '
                        'sample-dir). Workflow managers pass a writable task dir '
                        'so the input directory is never written to.')
    p.add_argument('--out', default=None,
                   help='Explicit output path for the prepared FASTQ. If given, '
                        'the prepared file is copied here (used by Nextflow).')
    p.add_argument('--out-filter-stats', default=None,
                   help='Explicit output path for length_filter_stats.json '
                        '(only relevant when --min-length > 0).')

    # align-vector
    p = subparsers.add_parser('align-vector', help='Align reads to the vector reference')
    p.add_argument('--reads', required=True)
    p.add_argument('--vector', required=True)
    p.add_argument('--out-sam', required=True)
    p.add_argument('--threads', type=int, default=4)

    # extract-vector
    p = subparsers.add_parser('extract-vector',
                              help='Extract mapped/unmapped reads and build vector BAM')
    p.add_argument('--sam', required=True)
    p.add_argument('--reads', required=True)
    p.add_argument('--out-fastq', required=True)
    p.add_argument('--out-bam', required=True)
    p.add_argument('--out-unmap', required=True)

    # classify-hq-lq
    p = subparsers.add_parser('classify-hq-lq',
                              help='Split vector reads into HQ/LQ by coverage')
    p.add_argument('--sam', required=True)
    p.add_argument('--vector-fastq', required=True)
    p.add_argument('--out-hq', required=True)
    p.add_argument('--out-lq', required=True)
    p.add_argument('--min-coverage', type=float, default=0.9)

    # align-ref (generic)
    p = subparsers.add_parser('align-ref',
                              help='Align reads to a single reference -> sorted BAM')
    p.add_argument('--reads', required=True)
    p.add_argument('--ref', required=True)
    p.add_argument('--ref-name', required=True)
    p.add_argument('--out-bam', required=True)
    p.add_argument('--threads', type=int, default=4)

    # classify-best (generic)
    p = subparsers.add_parser('classify-best',
                              help='Assign reads to best-matching reference across BAMs')
    p.add_argument('--bams', nargs='+', required=True,
                   help="One or more 'refname=path' pairs")
    p.add_argument('--source-fastq', required=True)
    p.add_argument('--out-prefix', required=True)
    p.add_argument('--mode', default='hybrid', choices=['hybrid', 'impurity'])

    # compute-coverage
    p = subparsers.add_parser('compute-coverage',
                              help='Compute per-position coverage for one BAM')
    p.add_argument('--bam', required=True)
    p.add_argument('--category', required=True)
    p.add_argument('--out-json', required=True)
    p.add_argument('--max-pos', type=int, default=6000)

    # compute-truncation
    p = subparsers.add_parser('compute-truncation',
                              help='Compute truncation hotspots for the vector BAM')
    p.add_argument('--bam', required=True)
    p.add_argument('--out-csv', required=True)
    p.add_argument('--out-json', required=True)

    # compute-lengths
    p = subparsers.add_parser('compute-lengths',
                              help='Compute read-length distributions for named FASTQs')
    p.add_argument('--fastqs', nargs='+', required=True,
                   help="One or more 'name=path' pairs")
    p.add_argument('--out-json', required=True)

    # summarize
    p = subparsers.add_parser('summarize',
                              help='Aggregate classification statistics')
    p.add_argument('--fastq-dir', required=True)
    p.add_argument('--optional-refs', nargs='*', default=[],
                   help='Optional reference names present (e.g. rep_cap spike_in)')
    p.add_argument('--out-json', required=True)
    p.add_argument('--out-txt', default=None)

    # taxonomy (combined kraken + bracken + summary)
    p = subparsers.add_parser('taxonomy',
                              help='Run full taxonomy (Kraken2 + Bracken) and write summary JSON')
    p.add_argument('--reads', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--out-summary', required=True,
                   help='Path for taxonomy_summary.json (consumed by the report)')
    p.add_argument('--confidence', type=float, default=0.05)
    p.add_argument('--level', default='S', choices=['D', 'P', 'C', 'O', 'F', 'G', 'S'])
    p.add_argument('--threshold', type=int, default=10)
    p.add_argument('--read-length', type=int, default=None,
                   help='Bracken read length (-r) to match the data read-length '
                        'distribution. Must be a length the Bracken DB was built '
                        'for. Default: none (uses the DB maximum).')
    p.add_argument('--threads', type=int, default=4)

    # taxonomy-kraken
    p = subparsers.add_parser('taxonomy-kraken', help='Run Kraken2 classification')
    p.add_argument('--reads', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--confidence', type=float, default=0.05)
    p.add_argument('--threads', type=int, default=4)

    # taxonomy-bracken
    p = subparsers.add_parser('taxonomy-bracken', help='Run Bracken re-estimation')
    p.add_argument('--kraken-report', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--out-dir', required=True)
    p.add_argument('--level', default='S', choices=['D', 'P', 'C', 'O', 'F', 'G', 'S'])
    p.add_argument('--threshold', type=int, default=10)

    # generate-report
    p = subparsers.add_parser('generate-report', help='Assemble the HTML report')
    p.add_argument('--stats', required=True)
    p.add_argument('--truncation', required=True)
    p.add_argument('--coverage', required=True)
    p.add_argument('--lengths', required=True)
    p.add_argument('--taxonomy', default=None)
    p.add_argument('--out', required=True)


def parse_arguments():
    """Parse command-line arguments using subcommands."""
    parser = argparse.ArgumentParser(
        description='AAVision — AAV analysis pipeline for Nanopore sequencing data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Subcommands:
  run               Full multi-sample pipeline (most users want this)

  Fine-grained steps (for debugging / Nextflow orchestration):
  prepare-fastq     Merge + optional length filter for one sample
  align-vector      Align reads to the vector reference
  extract-vector    Extract mapped/unmapped reads, build vector BAM
  classify-hq-lq    Split vector reads into HQ/LQ by coverage
  align-ref         Align reads to a single reference -> sorted BAM
  classify-best     Assign reads to best-matching reference across BAMs
  compute-coverage  Per-position coverage for one BAM
  compute-truncation  Truncation hotspots for the vector BAM
  compute-lengths   Read-length distributions for named FASTQs
  summarize         Aggregate classification statistics
  taxonomy-kraken   Kraken2 classification of reads
  taxonomy-bracken  Bracken abundance re-estimation
  generate-report   Assemble the final HTML report

Run `aav_analyze.py <subcommand> --help` for per-subcommand options.
        """
    )

    subparsers = parser.add_subparsers(dest='command', metavar='<subcommand>')

    _add_run_subparser(subparsers)
    _add_step_subparsers(subparsers)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        parser.exit(1)

    return args


def analyze_single_sample(fastq_path: str, sample_name: str, output_dir: str,
                          refs: dict, quality: float, threads: int, 
                          verbose: bool, kraken_db: str = None,
                          kraken_confidence: float = 0.05,
                          bracken_level: str = 'S',
                          bracken_threshold: int = 10) -> dict:
    """
    Run full analysis pipeline for a single sample.
    
    Args:
        fastq_path: Path to input FASTQ file
        sample_name: Name for this sample
        output_dir: Output directory for this sample
        refs: Dictionary of reference file paths
        quality: Minimum coverage threshold (0-100)
        threads: Number of threads for alignment
        verbose: Enable verbose logging
        kraken_db: Path to Kraken2 database (optional, skips taxonomy if None)
        kraken_confidence: Kraken2 confidence threshold (default: 0.05)
        bracken_level: Bracken taxonomic level (default: S)
        bracken_threshold: Bracken minimum read count (default: 10)
        
    Returns:
        Dictionary with analysis statistics
    """
    log = logging.getLogger('aav_analyzer')
    
    # Create output directory structure
    dirs = utils.create_output_structure(output_dir)
    log.info(f"Created output directory: {output_dir}")
    
    # Step 1: Align reads to vector reference
    log.info("="*60)
    log.info("STEP 1: Aligning reads to vector reference")
    log.info("="*60)
    
    vector_sam = os.path.join(dirs['bam'], 'vector.sam')
    vector_fastq = os.path.join(dirs['fastq'], 'vector.fastq')
    vector_bam = os.path.join(dirs['bam'], 'vector.bam')
    
    alignment.align_reads(fastq_path, refs['vector'], vector_sam, threads)
    alignment.extract_mapped_reads(vector_sam, vector_fastq)
    alignment.create_bam_for_mapped_reads(vector_sam, vector_bam)
    
    # Step 2: Classify vector reads by quality
    log.info("="*60)
    log.info("STEP 2: Classifying vector reads by quality")
    log.info("="*60)
    
    min_coverage = quality / 100.0
    hq_fastq, lq_fastq, hq_count, lq_count = classification.classify_vector_reads(
        vector_sam, vector_fastq, dirs['fastq'], min_coverage
    )
    
    # Extract unmapped reads
    unmap_fastq = os.path.join(dirs['fastq'], 'vector_unmap.fastq')
    alignment.create_fastq_difference(fastq_path, vector_fastq, unmap_fastq)
    
    # Clean up SAM file after classification is complete
    os.remove(vector_sam)
    
    # Step 3: Process hybrid reads
    log.info("="*60)
    log.info("STEP 3: Processing hybrid reads (low-quality vector reads)")
    log.info("="*60)
    
    hybrid_refs = {k: v for k, v in refs.items() 
                  if k in ['backbone', 'helper', 'host', 'rep_cap']}
    
    # Create hybrid output directly in fastq directory
    hybrid_bam_prefix = os.path.join(dirs['bam'], 'alignments_hybrid')
    hybrid_fastq_prefix = os.path.join(dirs['fastq'], 'hybrid')
    
    # Align to multiple references
    hybrid_bams = alignment.align_to_multiple_references(lq_fastq, hybrid_refs, hybrid_bam_prefix, threads)
    
    # Classify reads directly to fastq directory
    hybrid_files = classification.classify_reads_from_bams(hybrid_bams, lq_fastq, hybrid_fastq_prefix, "hybrid")
    
    log.info(f"Hybrid classification complete: {len(hybrid_files)} categories")
    
    # Step 4: Process impurity reads
    log.info("="*60)
    log.info("STEP 4: Processing impurity reads (unmapped reads)")
    log.info("="*60)
    
    impurity_refs = {k: v for k, v in refs.items() 
                    if k in ['backbone', 'helper', 'host', 'rep_cap', 'spike_in']}
    
    # Create impurity output directly in fastq directory
    impurity_bam_prefix = os.path.join(dirs['bam'], 'alignments_impurity')
    impurity_fastq_prefix = os.path.join(dirs['fastq'], 'impurity')
    
    # Align to multiple references
    impurity_bams = alignment.align_to_multiple_references(unmap_fastq, impurity_refs, impurity_bam_prefix, threads)
    
    # Classify reads directly to fastq directory
    impurity_files = classification.classify_reads_from_bams(impurity_bams, unmap_fastq, impurity_fastq_prefix, "impurity")
    
    # Step 5: Generate statistics
    log.info("="*60)
    log.info("STEP 5: Generating statistics")
    log.info("="*60)
    
    optional_refs = [k for k in refs.keys() if k in ['rep_cap', 'spike_in']]
    stats = classification.summarize_classification(dirs['fastq'], optional_refs)
    
    # Save statistics to file
    stats_file = os.path.join(dirs['stats'], 'classification_stats.json')
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    log.info(f"Statistics saved to {stats_file}")
    
    # Also save human-readable summary
    summary_file = os.path.join(dirs['stats'], 'summary.txt')
    with open(summary_file, 'w') as f:
        f.write("AAV Analysis Summary\n")
        f.write("="*60 + "\n\n")
        f.write(f"Sample: {sample_name}\n\n")
        f.write(f"Total reads: {stats['total_reads']:,}\n")
        f.write(f"  Vector reads: {stats['total_vector']:,} ({stats['vector_percent']:.2f}%)\n")
        f.write(f"    - High quality: {stats['vector_hq']:,}\n")
        f.write(f"    - Low quality: {stats['vector_lq_pure']:,}\n")
        f.write(f"  Hybrid reads: {stats['total_hybrid']:,} ({stats['hybrid_percent']:.2f}%)\n")
        
        # Add hybrid breakdown
        for key in ['backbone', 'helper', 'host', 'rep_cap']:
            hybrid_key = f'{key}_hybrid'
            if hybrid_key in stats:
                f.write(f"    - {key.capitalize()}: {stats[hybrid_key]:,}\n")
        
        f.write(f"  Impurity reads: {stats['total_impurity']:,} ({stats['impurity_percent']:.2f}%)\n")
        
        # Add impurity breakdown
        for key in ['backbone', 'helper', 'host', 'rep_cap', 'spike_in']:
            impurity_key = f'{key}_impurity'
            if impurity_key in stats:
                f.write(f"    - {key.capitalize()}: {stats[impurity_key]:,}\n")
    
    log.info(f"Summary saved to {summary_file}")
    
    log.info(f"Total reads: {stats['total_reads']:,}")
    log.info(f"  Vector reads: {stats['total_vector']:,} ({stats['vector_percent']:.2f}%)")
    log.info(f"    - High quality: {stats['vector_hq']:,}")
    log.info(f"    - Low quality: {stats['vector_lq_pure']:,}")
    log.info(f"  Hybrid reads: {stats['total_hybrid']:,} ({stats['hybrid_percent']:.2f}%)")
    log.info(f"  Impurity reads: {stats['total_impurity']:,} ({stats['impurity_percent']:.2f}%)")
    
    # Step 6: Analyze truncation hotspots and coverage
    log.info("="*60)
    log.info("STEP 6: Analyzing truncation hotspots and coverage")
    log.info("="*60)
    
    # Compute truncation for vector (rAAV) only
    trunc_df = truncation.compute_truncation_hotspots(vector_bam)
    trunc_data = truncation.get_truncation_data_dict(trunc_df)
    
    # Save truncation data
    trunc_file = os.path.join(dirs['stats'], 'truncation_hotspots.csv')
    trunc_df.to_csv(trunc_file, index=False)
    log.info(f"Truncation data saved to {trunc_file}")
    
    # Save truncation JSON for plotting
    trunc_json_file = os.path.join(dirs['stats'], 'truncation_hotspots.json')
    with open(trunc_json_file, 'w') as f:
        json.dump(trunc_data, f, indent=2)
    
    # Compute coverage for rAAV and impurities (not host genome)
    all_coverage_data = {}
    
    # Vector (rAAV) coverage
    all_coverage_data['vector'] = truncation.compute_coverage(vector_bam)
    log.info("Computed coverage for rAAV")
    
    # Helper impurity coverage
    helper_bam = os.path.join(dirs['bam'], 'alignments_impurity_helper.bam')
    if os.path.exists(helper_bam):
        try:
            all_coverage_data['helper'] = truncation.compute_coverage(helper_bam)
            log.info("Computed coverage for helper impurity")
        except Exception as e:
            log.warning(f"Could not compute coverage for helper: {e}")
    
    # Rep/Cap impurity coverage
    if 'rep_cap' in optional_refs:
        rep_cap_bam = os.path.join(dirs['bam'], 'alignments_impurity_rep_cap.bam')
        if os.path.exists(rep_cap_bam):
            try:
                all_coverage_data['rep_cap'] = truncation.compute_coverage(rep_cap_bam)
                log.info("Computed coverage for rep_cap impurity")
            except Exception as e:
                log.warning(f"Could not compute coverage for rep_cap: {e}")
    
    # Spike-in impurity coverage
    if 'spike_in' in optional_refs:
        spike_in_bam = os.path.join(dirs['bam'], 'alignments_impurity_spike_in.bam')
        if os.path.exists(spike_in_bam):
            try:
                all_coverage_data['spike_in'] = truncation.compute_coverage(spike_in_bam)
                log.info("Computed coverage for spike_in impurity")
            except Exception as e:
                log.warning(f"Could not compute coverage for spike_in: {e}")
    
    # Save coverage data as JSON
    coverage_json_file = os.path.join(dirs['stats'], 'coverage.json')
    with open(coverage_json_file, 'w') as f:
        json.dump(all_coverage_data, f, indent=2)
    log.info(f"Coverage data saved to {coverage_json_file}")
    
    # Step 7: Analyze read length distributions
    log.info("="*60)
    log.info("STEP 7: Analyzing read length distributions")
    log.info("="*60)
    
    # Build list of FASTQ files for length analysis
    fastq_files = [('vector', vector_fastq)]
    
    # Add helper impurity
    helper_impurity = os.path.join(dirs['fastq'], 'impurity_helper.fastq')
    if os.path.exists(helper_impurity) and os.path.getsize(helper_impurity) > 0:
        fastq_files.append(('helper', helper_impurity))
    
    # Add rep_cap impurity if present
    if 'rep_cap' in optional_refs:
        rep_cap_impurity = os.path.join(dirs['fastq'], 'impurity_rep_cap.fastq')
        if os.path.exists(rep_cap_impurity) and os.path.getsize(rep_cap_impurity) > 0:
            fastq_files.append(('rep_cap', rep_cap_impurity))
    
    # Add host genome impurity
    host_impurity = os.path.join(dirs['fastq'], 'impurity_host.fastq')
    if os.path.exists(host_impurity) and os.path.getsize(host_impurity) > 0:
        fastq_files.append(('host', host_impurity))
    
    # Add spike_in impurity if present
    if 'spike_in' in optional_refs:
        spike_in_impurity = os.path.join(dirs['fastq'], 'impurity_spike_in.fastq')
        if os.path.exists(spike_in_impurity) and os.path.getsize(spike_in_impurity) > 0:
            fastq_files.append(('spike_in', spike_in_impurity))
    
    # Add unmapped impurity reads (reads that didn't map to any reference)
    impurity_unmapped = os.path.join(dirs['fastq'], 'impurity_unmapped.fastq')
    if os.path.exists(impurity_unmapped) and os.path.getsize(impurity_unmapped) > 0:
        fastq_files.append(('unmapped', impurity_unmapped))
    
    length_distributions = truncation.analyze_read_lengths(fastq_files)
    
    # Save length distributions
    length_file = os.path.join(dirs['stats'], 'read_length_distributions.json')
    with open(length_file, 'w') as f:
        json.dump(length_distributions, f, indent=2)
    log.info(f"Read length distributions saved to {length_file}")
    
    # Step 8: Taxonomic classification of unmapped impurity reads (optional)
    taxonomy_data = None
    if kraken_db:
        log.info("="*60)
        log.info("STEP 8: Taxonomic classification of unmapped impurity reads")
        log.info("="*60)
        
        if not taxonomy.check_kraken_available():
            log.warning("kraken2 not found on PATH — skipping taxonomy step")
        elif not taxonomy.validate_kraken_db(kraken_db):
            log.warning(f"Invalid Kraken2 database at {kraken_db} — skipping taxonomy step")
        else:
            impurity_unmapped = os.path.join(dirs['fastq'], 'impurity_unmapped.fastq')
            
            if not os.path.exists(impurity_unmapped) or os.path.getsize(impurity_unmapped) == 0:
                log.warning("No impurity_unmapped.fastq reads available — skipping taxonomy step")
            else:
                taxonomy_dir = os.path.join(dirs['stats'], 'taxonomy')
                os.makedirs(taxonomy_dir, exist_ok=True)
                
                try:
                    taxonomy_data = taxonomy.run_taxonomy_analysis(
                        fastq_file=impurity_unmapped,
                        db_path=kraken_db,
                        output_dir=taxonomy_dir,
                        threads=threads,
                        confidence=kraken_confidence,
                        bracken_level=bracken_level,
                        bracken_threshold=bracken_threshold
                    )
                    
                    # Save taxonomy summary as JSON
                    tax_summary_file = os.path.join(taxonomy_dir, 'taxonomy_summary.json')
                    with open(tax_summary_file, 'w') as f:
                        json.dump(taxonomy_data, f, indent=2)
                    log.info(f"Taxonomy summary saved to {tax_summary_file}")
                    
                    stats_line = taxonomy_data.get('kraken_stats', {})
                    log.info(f"Kraken2: {stats_line.get('classified', 0):,} classified "
                             f"({stats_line.get('classified_percent', 0):.2f}%), "
                             f"{stats_line.get('unclassified', 0):,} unclassified "
                             f"({stats_line.get('unclassified_percent', 0):.2f}%)")
                    
                    if taxonomy_data.get('bracken_results'):
                        log.info(f"Bracken: top taxa at {bracken_level}-level available in report")
                    
                except Exception as e:
                    log.warning(f"Taxonomy analysis failed (non-fatal): {e}")
                    taxonomy_data = None
    
    # Step 9: Generate HTML report
    log.info("="*60)
    log.info("STEP 9: Generating HTML report")
    log.info("="*60)
    
    report_file = os.path.join(dirs['base'], 'analysis_report.html')
    report.generate_html_report(
        stats,
        trunc_data,
        length_distributions,
        report_file,
        coverage_data=all_coverage_data,
        taxonomy_data=taxonomy_data
    )
    
    return {
        'sample_name': sample_name,
        'stats': stats,
        'report_file': report_file,
        'dirs': dirs
    }


def main():
    """Entry point — parse args and dispatch to the right subcommand."""
    args = parse_arguments()

    # Configure logging early (verbose only meaningful for 'run', but harmless)
    verbose = getattr(args, 'verbose', False)
    utils.setup_logging(verbose)

    if args.command == 'run':
        return run_pipeline(args)

    # Fine-grained subcommands dispatch through the cli handler registry
    from aav_analyzer import cli
    handler = cli.HANDLERS.get(args.command)
    if handler is None:
        print(f"Error: unknown subcommand '{args.command}'", file=sys.stderr)
        return 1
    return handler(args)


def run_pipeline(args):
    """Run the full multi-sample pipeline (the 'run' subcommand)."""

    # ------------------------------------------------------------------
    # Validate run directory exists
    # ------------------------------------------------------------------
    if not os.path.isdir(args.run_dir):
        print(f"Error: Run directory not found: {args.run_dir}", file=sys.stderr)
        sys.exit(1)

    run_dir = os.path.abspath(args.run_dir)

    # ------------------------------------------------------------------
    # Handle output directory (--force vs error on existing)
    # ------------------------------------------------------------------
    if os.path.exists(args.output):
        if args.force:
            print(f"Removing existing output directory: {args.output}")
            shutil.rmtree(args.output)
        else:
            print(f"Error: Output directory already exists: {args.output}",
                  file=sys.stderr)
            print("Use --force to overwrite.", file=sys.stderr)
            sys.exit(1)

    os.makedirs(args.output, exist_ok=True)

    # Setup logging (run-level log captures everything)
    log_file = os.path.join(args.output, 'pipeline.log')
    log = utils.setup_logging(args.verbose, log_file)

    # ------------------------------------------------------------------
    # Discover samples
    # ------------------------------------------------------------------
    try:
        samples = utils.discover_samples(run_dir)
    except (FileNotFoundError, ValueError) as e:
        log.error(str(e))
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if not samples:
        log.error(f"No sample directories found in {run_dir}/samples/")
        print(f"Error: No sample directories found in {run_dir}/samples/",
              file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------
    # Print run banner
    # ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("AAVision — AAV analysis pipeline")
    print("=" * 60)
    print(f"Run directory:       {run_dir}")
    print(f"Output directory:    {args.output}")
    print(f"Number of samples:   {len(samples)}")
    print(f"Coverage threshold:  {args.quality}%")
    print(f"Threads:             {args.threads}")
    if args.kraken_db:
        print(f"Kraken2 database:    {args.kraken_db}")
    print("=" * 60 + "\n")

    log.info(f"Run directory: {run_dir}")
    log.info(f"Discovered {len(samples)} sample(s): "
             f"{', '.join(s['name'] for s in samples)}")

    # ------------------------------------------------------------------
    # Check required tools once for the whole run
    # ------------------------------------------------------------------
    log.info("Checking dependencies...")
    if not utils.check_required_tools():
        sys.exit(1)

    # ------------------------------------------------------------------
    # Process each sample
    # ------------------------------------------------------------------
    total_start_time = time.time()
    succeeded = []
    failed = []

    for i, sample in enumerate(samples, start=1):
        sample_name = sample['name']
        sample_dir = sample['dir']

        print(f"\n{'=' * 60}")
        print(f"Sample {i}/{len(samples)}: {sample_name}")
        print(f"{'=' * 60}")

        log.info(f"=== Processing sample {i}/{len(samples)}: {sample_name} ===")
        log.info(f"Sample directory: {sample_dir}")

        sample_start_time = time.time()

        try:
            # Resolve references (per-file fallback: sample → run-level)
            log.info(f"[{sample_name}] Resolving references:")
            refs = utils.resolve_sample_references(sample_dir, run_dir, sample_name)
            if refs is None:
                failed.append({
                    'name': sample_name,
                    'reason': 'Missing required references'
                })
                continue

            # Prepare FASTQ (single file or merge multiple)
            log.info(f"[{sample_name}] Preparing FASTQ input...")
            fastq_path = utils.prepare_sample_fastq(sample_dir, sample_name,
                                                    min_length=args.min_length)
            if fastq_path is None:
                failed.append({
                    'name': sample_name,
                    'reason': 'No FASTQ files found'
                })
                continue

            if not utils.validate_fastq(fastq_path):
                failed.append({
                    'name': sample_name,
                    'reason': f'Invalid FASTQ format: {fastq_path}'
                })
                continue

            # Per-sample output directory
            sample_output_dir = os.path.join(args.output, sample_name)

            # Run analysis
            result = analyze_single_sample(
                fastq_path=fastq_path,
                sample_name=sample_name,
                output_dir=sample_output_dir,
                refs=refs,
                quality=args.quality,
                threads=args.threads,
                verbose=args.verbose,
                kraken_db=args.kraken_db,
                kraken_confidence=args.kraken_confidence,
                bracken_level=args.bracken_level,
                bracken_threshold=args.bracken_threshold
            )

            sample_elapsed = time.time() - sample_start_time
            result['elapsed_seconds'] = sample_elapsed

            log.info(f"[{sample_name}] Completed in {sample_elapsed:.1f} seconds")
            print(f"\n[{sample_name}] Done in {sample_elapsed:.1f}s")
            print(f"  Report: {result['report_file']}")

            succeeded.append(result)

        except KeyboardInterrupt:
            log.error("Analysis interrupted by user")
            print("\nInterrupted by user", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            log.error(f"[{sample_name}] Failed: {e}", exc_info=args.verbose)
            print(f"[{sample_name}] FAILED: {e}", file=sys.stderr)
            failed.append({'name': sample_name, 'reason': str(e)})
            continue

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    total_elapsed = time.time() - total_start_time

    summary_path = os.path.join(args.output, 'run_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("AAVision — Run Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Run directory:    {run_dir}\n")
        f.write(f"Output directory: {os.path.abspath(args.output)}\n")
        f.write(f"Total samples:    {len(samples)}\n")
        f.write(f"Succeeded:        {len(succeeded)}\n")
        f.write(f"Failed:           {len(failed)}\n")
        f.write(f"Total time:       {total_elapsed:.1f} seconds\n")
        f.write("\n")

        if succeeded:
            f.write("Succeeded samples\n")
            f.write("-" * 60 + "\n")
            for r in succeeded:
                stats = r.get('stats', {})
                f.write(f"  {r['sample_name']}  "
                        f"({r.get('elapsed_seconds', 0):.1f}s)\n")
                if stats:
                    f.write(f"    total reads:  {stats.get('total_reads', 0):,}\n")
                    f.write(f"    rAAV:         "
                            f"{stats.get('vector_percent', 0):.2f}%\n")
                    f.write(f"    hybrid:       "
                            f"{stats.get('hybrid_percent', 0):.2f}%\n")
                    f.write(f"    impurity:     "
                            f"{stats.get('impurity_percent', 0):.2f}%\n")
                f.write(f"    report:       {r['report_file']}\n\n")

        if failed:
            f.write("Failed samples\n")
            f.write("-" * 60 + "\n")
            for entry in failed:
                f.write(f"  {entry['name']}: {entry['reason']}\n")
            f.write("\n")

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"Total time:        {total_elapsed:.1f} seconds")
    print(f"Samples succeeded: {len(succeeded)}/{len(samples)}")
    if failed:
        print(f"Samples failed:    {len(failed)}")
        for entry in failed:
            print(f"  - {entry['name']}: {entry['reason']}")
    print(f"\nRun summary:  {summary_path}")
    print(f"Pipeline log: {log_file}")
    print("=" * 60 + "\n")

    log.info(f"Run complete: {len(succeeded)}/{len(samples)} succeeded "
             f"in {total_elapsed:.1f}s")

    # Non-zero exit if anything failed (useful for scripting)
    sys.exit(0 if not failed else 1)


if __name__ == '__main__':
    sys.exit(main())
