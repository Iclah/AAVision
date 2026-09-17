"""
Read classification functions for AAV analysis.
Classifies hybrid and impurity reads based on alignment scores.
"""

import os
import logging
import pysam
from typing import Dict, List
from collections import defaultdict


def classify_reads_from_bams(bam_files: Dict[str, str], source_fastq: str, 
                             output_prefix: str, read_type: str = "hybrid") -> Dict[str, str]:
    """
    Classify reads based on best alignment score across multiple BAM files.
    
    Args:
        bam_files: Dict mapping reference names to BAM file paths
        source_fastq: Source FASTQ file containing reads to classify
        output_prefix: Prefix for output FASTQ files
        read_type: Type of reads being classified ("hybrid" or "impurity")
        
    Returns:
        Dict mapping categories to output FASTQ file paths
    """
    logging.info(f"Classifying {read_type} reads")
    
    # Dictionary to store best alignment for each read
    read_scores = {}

    # Process BAMs in a fixed canonical reference order so that tie-breaking
    # (equal alignment scores) is DETERMINISTIC and independent of the order
    # the BAMs happen to be passed in. The scoring loop below keeps the first
    # reference seen on a score tie (uses strict <), so the iteration order
    # decides tied reads. Without this sort, the standalone pipeline (fixed
    # dict order) and Nextflow (groupTuple order) could assign tied reads to
    # different references, producing different per-category counts for the
    # same total. The canonical order matches the standalone's reference order.
    _CANONICAL_ORDER = ['backbone', 'helper', 'host', 'rep_cap', 'spike_in']

    def _rank(name):
        return _CANONICAL_ORDER.index(name) if name in _CANONICAL_ORDER else len(_CANONICAL_ORDER)

    ordered_items = sorted(bam_files.items(), key=lambda kv: _rank(kv[0]))

    # Parse each BAM file to get alignment scores
    for ref_name, bam_file in ordered_items:
        logging.debug(f"Parsing alignments from {ref_name}")
        bam = pysam.AlignmentFile(bam_file, "rb")

        for read in bam:
            if not read.is_unmapped:
                read_name = read.query_name
                try:
                    score = read.get_tag("AS")  # Alignment score
                except KeyError:
                    score = 0

                # Keep best score for each read. Strict < means the FIRST
                # reference (in canonical order) wins on a tie.
                if read_name not in read_scores or read_scores[read_name][1] < score:
                    read_scores[read_name] = (ref_name, score)

        bam.close()
    
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
                         output_dir: str, min_coverage: float = 0.9) -> tuple:
    """
    Classify vector reads into high-quality and low-quality based on coverage.
    
    Args:
        vector_sam: SAM file from vector alignment
        vector_fastq: FASTQ file with mapped vector reads
        output_dir: Output directory for classified reads
        min_coverage: Minimum coverage threshold for HQ reads
        
    Returns:
        Tuple of (HQ_fastq_path, LQ_fastq_path, HQ_count, LQ_count)
    """
    from .alignment import filter_reads_by_coverage, create_fastq_difference
    
    logging.info("Classifying vector reads by quality")
    
    hq_fastq = os.path.join(output_dir, 'vector_HQ.fastq')
    lq_fastq = os.path.join(output_dir, 'vector_LQ.fastq')
    
    # Extract high-quality reads (≥90% coverage)
    hq_count = filter_reads_by_coverage(vector_sam, hq_fastq, min_coverage)
    
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
