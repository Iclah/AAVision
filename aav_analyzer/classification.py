"""
Read classification functions for AAV analysis.
Classifies hybrid and impurity reads based on alignment scores.
"""

import os
import json
import logging
import pysam
from typing import Dict, List
from collections import defaultdict

from .homology import SharedRegions, unique_aligned_bases

logger = logging.getLogger(__name__)


def classify_reads_from_bams(bam_files: Dict[str, str], source_fastq: str,
                             output_prefix: str, read_type: str = "hybrid",
                             shared_regions: SharedRegions | None = None,
                             min_unique_bases: int = 100) -> Dict[str, str]:
    """
    Classify reads based on best alignment score across multiple BAM files.

    Without shared_regions every mapped record is eligible (v1.0.1 behaviour).
    With shared_regions (hybrid pass, docs/DECISIONS.md), a reference is
    eligible for a read only if the read's records on it align at least
    min_unique_bases reference bases outside that reference's vector-shared
    regions, summed over its primary and supplementary records
    (homology.unique_aligned_bases; 0-based half-open intervals, strand
    ignored). The best eligible reference by AS wins; reads with no eligible
    reference go to unmapped (pure LQ in the hybrid pass). Rejected
    (read, reference) pairs are written to {output_prefix}_shared_only.tsv and
    counts to {output_prefix}_filter_stats.json.

    Args:
        bam_files: Dict mapping reference names to BAM file paths
        source_fastq: Source FASTQ file containing reads to classify
        output_prefix: Prefix for output FASTQ files
        read_type: Type of reads being classified ("hybrid" or "impurity")
        shared_regions: Vector-shared regions per reference, or None to disable
            the filter. Must have searched every reference in bam_files.
        min_unique_bases: Minimum aligned reference bases outside shared
            regions for a reference to be eligible (inclusive)

    Returns:
        Dict mapping categories to output FASTQ file paths

    Raises:
        ValueError: if shared_regions was not searched for a reference in bam_files
    """
    logging.info(f"Classifying {read_type} reads")

    if shared_regions is not None:
        unsearched = sorted(set(bam_files) - shared_regions.searched)
        if unsearched:
            raise ValueError(
                f"Shared regions were not computed for reference(s) {unsearched}; "
                f"searched: {sorted(shared_regions.searched)}"
            )
        logger.info(
            f"Excluding vector-shared regions: a reference counts only with "
            f">= {min_unique_bases} aligned bases outside them"
        )

    # Process BAMs in a fixed canonical reference order so that tie-breaking
    # (equal alignment scores) is DETERMINISTIC and independent of the order
    # the BAMs happen to be passed in. The selection loop below keeps the first
    # reference seen on a score tie (uses strict <), so the iteration order
    # decides tied reads. Without this sort, the standalone pipeline (fixed
    # dict order) and Nextflow (groupTuple order) could assign tied reads to
    # different references, producing different per-category counts for the
    # same total. The canonical order matches the standalone's reference order.
    _CANONICAL_ORDER = ['backbone', 'helper', 'host', 'rep_cap', 'spike_in']

    def _rank(name):
        return _CANONICAL_ORDER.index(name) if name in _CANONICAL_ORDER else len(_CANONICAL_ORDER)

    ordered_items = sorted(bam_files.items(), key=lambda kv: _rank(kv[0]))

    # read name -> ref name -> [best AS, aligned bases outside shared regions]
    hits = {}
    for ref_name, bam_file in ordered_items:
        logging.debug(f"Parsing alignments from {ref_name}")
        shared = shared_regions.intervals.get(ref_name, {}) if shared_regions else None
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                if read.is_unmapped:
                    continue
                try:
                    score = read.get_tag("AS")  # Alignment score
                except KeyError:
                    score = 0
                hit = hits.setdefault(read.query_name, {}).setdefault(ref_name, [score, 0])
                hit[0] = max(hit[0], score)
                if shared is not None:
                    hit[1] += unique_aligned_bases(read, shared)

    # Best eligible reference per read. Strict < means the FIRST reference
    # (in canonical order) wins on a tie.
    read_scores = {}
    rejected = []  # (read, ref, best AS, unique bases)
    for read_name, refs in hits.items():
        for ref_name, _ in ordered_items:
            if ref_name not in refs:
                continue
            score, unique = refs[ref_name]
            if shared_regions is not None and unique < min_unique_bases:
                rejected.append((read_name, ref_name, score, unique))
                continue
            if read_name not in read_scores or read_scores[read_name][1] < score:
                read_scores[read_name] = (ref_name, score)

    if shared_regions is not None:
        _write_shared_filter_report(
            rejected, read_scores, list(bam_files), min_unique_bases, output_prefix
        )

    # Initialize read collections
    categories = list(bam_files.keys()) + ['unmapped']
    reads_classified = {cat: [] for cat in categories}
    
    # Classify reads from source FASTQ
    with pysam.FastxFile(source_fastq) as fq:
        for entry in fq:
            if entry.name in read_scores:
                best_ref = read_scores[entry.name][0]
                fastq_entry = f"@{entry.name}\n{entry.sequence}\n+\n{entry.quality}"
                reads_classified[best_ref].append(fastq_entry)
            else:
                fastq_entry = f"@{entry.name}\n{entry.sequence}\n+\n{entry.quality}"
                reads_classified['unmapped'].append(fastq_entry)
    
    # Write classified reads to separate FASTQ files
    output_files = {}
    
    for category, reads in reads_classified.items():
        if category == 'unmapped':
            output_file = f"{output_prefix}_unmapped.fastq"
        else:
            output_file = f"{output_prefix}_{category}.fastq"
        
        with open(output_file, 'w') as f:
            if reads:
                # Join records with newlines AND terminate the final record with
                # one, so the file ends in a newline. Without the trailing '\n'
                # the last record has no line terminator, which malforms the
                # FASTQ (breaks `wc -l`/4 counts and naive concatenation).
                f.write('\n'.join(reads) + '\n')
            # Empty categories still produce an empty file (no stray newline).
        
        output_files[category] = output_file
        logging.debug(f"Wrote {len(reads)} reads to {category}")
    
    logging.info(f"Classification complete: {len(read_scores)} reads assigned")
    return output_files


def _write_shared_filter_report(rejected: list, read_scores: dict, ref_names: list[str],
                                min_unique_bases: int, output_prefix: str) -> None:
    """
    Write {output_prefix}_shared_only.tsv (one row per rejected read/reference
    pair) and {output_prefix}_filter_stats.json, and log the counts.

    A read is "rescued to LQ" when it had alignments but every one was
    rejected; it is attributed to the reference of its best rejected AS
    (canonical order on ties, as rejected is built in that order).
    """
    tsv_path = f"{output_prefix}_shared_only.tsv"
    with open(tsv_path, 'w') as fh:
        fh.write("read_id\tref_name\tbest_AS\tunique_bases\tfinal_category\n")
        for read_name, ref_name, score, unique in rejected:
            final = read_scores[read_name][0] if read_name in read_scores else 'unmapped'
            fh.write(f"{read_name}\t{ref_name}\t{score}\t{unique}\t{final}\n")

    rejected_by_ref = {ref: 0 for ref in ref_names}
    best_rejected = {}  # read -> (ref, AS) for reads with no eligible reference
    for read_name, ref_name, score, _ in rejected:
        rejected_by_ref[ref_name] += 1
        if read_name in read_scores:
            continue
        if read_name not in best_rejected or best_rejected[read_name][1] < score:
            best_rejected[read_name] = (ref_name, score)
    rescued_by_ref = {ref: 0 for ref in ref_names}
    for ref_name, _ in best_rejected.values():
        rescued_by_ref[ref_name] += 1

    stats = {
        'min_unique_bases': min_unique_bases,
        'rejected_alignments_by_ref': rejected_by_ref,
        'rescued_to_lq_by_ref': rescued_by_ref,
        'rescued_to_lq': len(best_rejected),
    }
    with open(f"{output_prefix}_filter_stats.json", 'w') as fh:
        json.dump(stats, fh, indent=2)

    logger.info(
        f"Shared-region filter: {len(best_rejected)} reads with alignments only "
        f"inside vector-shared regions stay unassigned ({rescued_by_ref}); "
        f"rejected (read, reference) pairs: {rejected_by_ref}; details in {tsv_path}"
    )


def get_classification_stats(output_files: Dict[str, str]) -> Dict[str, int]:
    """
    Get read counts for each classification category.
    
    Args:
        output_files: Dict mapping categories to FASTQ file paths
        
    Returns:
        Dict mapping categories to read counts
    """
    stats = {}
    
    for category, filepath in output_files.items():
        count = 0
        if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
            with open(filepath) as f:
                count = sum(1 for line in f if line.startswith('@'))
        stats[category] = count
    
    return stats


def classify_vector_reads(vector_sam: str, vector_fastq: str,
                         output_dir: str, min_coverage: float = 0.9,
                         min_snapback_overlap: int = 200) -> tuple:
    """
    Classify vector reads into high-quality and low-quality based on coverage.
    
    Args:
        vector_sam: SAM file from vector alignment
        vector_fastq: FASTQ file with mapped vector reads
        output_dir: Output directory for classified reads
        min_coverage: Minimum merged primary+supplementary read coverage for HQ
        min_snapback_overlap: Opposite-strand reference overlap (bp) that makes
            a read a snapback (LQ)
        
    Returns:
        Tuple of (HQ_fastq_path, LQ_fastq_path, HQ_count, LQ_count)
    """
    from .alignment import filter_reads_by_coverage, create_fastq_difference
    
    logging.info("Classifying vector reads by quality")
    
    hq_fastq = os.path.join(output_dir, 'vector_HQ.fastq')
    lq_fastq = os.path.join(output_dir, 'vector_LQ.fastq')
    
    # HQ: primary + supplementary cover >= min_coverage of the read, not a snapback
    hq_count = filter_reads_by_coverage(
        vector_sam, vector_fastq, hq_fastq, min_coverage, min_snapback_overlap
    )
    
    # Create low-quality reads (vector_mapped - HQ)
    lq_count = create_fastq_difference(vector_fastq, hq_fastq, lq_fastq)
    
    logging.info(f"Vector reads classified: {hq_count} HQ, {lq_count} LQ")
    
    return hq_fastq, lq_fastq, hq_count, lq_count


def process_hybrid_reads(lq_fastq: str, references: dict, output_dir: str, 
                        threads: int = 4) -> Dict[str, str]:
    """
    Process low-quality vector reads to identify hybrids.
    
    Args:
        lq_fastq: FASTQ file with low-quality vector reads
        references: Dict of reference names to FASTA paths
        output_dir: Output directory
        threads: Number of threads for alignment
        
    Returns:
        Dict mapping hybrid categories to FASTQ file paths
    """
    from .alignment import align_to_multiple_references
    
    logging.info("Processing hybrid reads")
    
    # Align to multiple references
    bam_prefix = os.path.join(output_dir, 'alignments_hybrid')
    bam_files = align_to_multiple_references(lq_fastq, references, bam_prefix, threads)
    
    # Classify reads
    output_prefix = os.path.join(output_dir, 'hybrid')
    output_files = classify_reads_from_bams(bam_files, lq_fastq, output_prefix, "hybrid")
    
    # Rename unmapped file
    os.rename(
        os.path.join(output_dir, 'hybrid_unmapped.fastq'),
        os.path.join(output_dir, 'hybrid_unmapped.fastq')
    )
    
    return output_files


def process_impurity_reads(unmap_fastq: str, references: dict, 
                          output_dir: str, threads: int = 4) -> Dict[str, str]:
    """
    Process unmapped reads to identify impurities.
    
    Args:
        unmap_fastq: FASTQ file with unmapped reads
        references: Dict of reference names to FASTA paths
        output_dir: Output directory
        threads: Number of threads for alignment
        
    Returns:
        Dict mapping impurity categories to FASTQ file paths
    """
    from .alignment import align_to_multiple_references
    
    logging.info("Processing impurity reads")
    
    # Align to multiple references
    bam_prefix = os.path.join(output_dir, 'alignments_impurity')
    bam_files = align_to_multiple_references(unmap_fastq, references, bam_prefix, threads)
    
    # Classify reads
    output_prefix = os.path.join(output_dir, 'impurity')
    output_files = classify_reads_from_bams(bam_files, unmap_fastq, output_prefix, "impurity")
    
    # Rename unmapped file
    os.rename(
        os.path.join(output_dir, 'impurity_unmapped.fastq'),
        os.path.join(output_dir, 'impurity_unmapped.fastq')
    )
    
    return output_files


def summarize_classification(fastq_dir: str, optional_refs: List[str] = None) -> dict:
    """
    Generate summary statistics for all classified reads.
    
    Args:
        fastq_dir: Directory containing classified FASTQ files
        optional_refs: List of optional reference names (e.g., ['rep_cap', 'spike_in'])
        
    Returns:
        Dict with classification statistics
    """
    from .utils import count_fastq_reads
    
    logging.info("Generating classification summary")
    
    # Count all FASTQ files
    counts = {
        'vector_hq': count_fastq_reads(os.path.join(fastq_dir, 'vector_HQ.fastq')),
        'vector_lq': count_fastq_reads(os.path.join(fastq_dir, 'vector_LQ.fastq')),
        'backbone_hybrid': count_fastq_reads(os.path.join(fastq_dir, 'hybrid_backbone.fastq')),
        'helper_hybrid': count_fastq_reads(os.path.join(fastq_dir, 'hybrid_helper.fastq')),
        'host_hybrid': count_fastq_reads(os.path.join(fastq_dir, 'hybrid_host.fastq')),
        'backbone_impurity': count_fastq_reads(os.path.join(fastq_dir, 'impurity_backbone.fastq')),
        'helper_impurity': count_fastq_reads(os.path.join(fastq_dir, 'impurity_helper.fastq')),
        'host_impurity': count_fastq_reads(os.path.join(fastq_dir, 'impurity_host.fastq')),
    }
    
    # Add optional references if present
    if optional_refs:
        for ref in optional_refs:
            hybrid_file = os.path.join(fastq_dir, f'hybrid_{ref}.fastq')
            impurity_file = os.path.join(fastq_dir, f'impurity_{ref}.fastq')
            
            if os.path.exists(hybrid_file):
                counts[f'{ref}_hybrid'] = count_fastq_reads(hybrid_file)
            if os.path.exists(impurity_file):
                counts[f'{ref}_impurity'] = count_fastq_reads(impurity_file)
    
    # Calculate totals
    counts['total_hybrid'] = sum(v for k, v in counts.items() if 'hybrid' in k)
    counts['vector_lq_pure'] = counts['vector_lq'] - counts['total_hybrid']  # LQ reads that are NOT hybrids
    counts['total_vector'] = counts['vector_hq'] + counts['vector_lq_pure']  # HQ + pure LQ
    counts['total_impurity'] = sum(v for k, v in counts.items() if 'impurity' in k)
    counts['total_reads'] = counts['total_vector'] + counts['total_hybrid'] + counts['total_impurity']
    
    # Calculate percentages
    if counts['total_reads'] > 0:
        counts['vector_percent'] = (counts['total_vector'] / counts['total_reads']) * 100
        counts['hybrid_percent'] = (counts['total_hybrid'] / counts['total_reads']) * 100
        counts['impurity_percent'] = (counts['total_impurity'] / counts['total_reads']) * 100
    
    return counts
