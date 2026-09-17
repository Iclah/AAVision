"""
Alignment and read filtering functions for AAV analysis.
Handles minimap2 alignment and quality-based filtering.

SAM/BAM Flag Strategy:
- Primary alignments (flag 0): Main alignment for each read â†’ ANALYZE
- Secondary alignments (flag 256): Same read mapping elsewhere â†’ EXCLUDE from all files
- Supplementary alignments (flag 2048): Chimeric/split alignments â†’ KEEP in BAM for visualization, EXCLUDE from analysis

Flags used:
- -F 260: Exclude unmapped(4) + secondary(256) â†’ Keep primary + supplementary
- -F 2308: Exclude unmapped(4) + secondary(256) + supplementary(2048) â†’ Keep primary only
- -F 256: Exclude secondary(256) only â†’ Keep primary + supplementary + unmapped
"""

import os
import logging
import pysam
import subprocess
from typing import Tuple, List
from collections import defaultdict


def align_reads(fastq_file: str, reference_fasta: str, output_sam: str, 
                threads: int = 4) -> str:
    """
    Align reads to reference using minimap2.
    
    Args:
        fastq_file: Input FASTQ file
        reference_fasta: Reference FASTA file
        output_sam: Output SAM file path
        threads: Number of threads for minimap2
        
    Returns:
        Path to output SAM file
    """
    cmd = [
        'minimap2',
        '-ax', 'map-ont',
        '-t', str(threads),
        reference_fasta,
        fastq_file
    ]
    
    logging.info(f"Aligning reads to {os.path.basename(reference_fasta)}")
    
    with open(output_sam, 'w') as outfile:
        subprocess.run(cmd, stdout=outfile, check=True, stderr=subprocess.PIPE)
    
    logging.info(f"Alignment complete: {output_sam}")
    return output_sam


def sam_to_sorted_bam(sam_file: str, output_bam: str, 
                     exclude_secondary: bool = True) -> str:
    """
    Convert SAM to sorted BAM and create index.
    
    Args:
        sam_file: Input SAM file
        output_bam: Output BAM file path
        exclude_secondary: If True, exclude secondary alignments (keep primary + supplementary)
        
    Returns:
        Path to output BAM file
    """
    logging.info(f"Converting {os.path.basename(sam_file)} to sorted BAM")
    
    # Create temporary unsorted BAM
    temp_bam = output_bam + ".tmp.bam"
    
    try:
        # Step 1: Convert SAM to BAM
        if exclude_secondary:
            # -F 256 excludes secondary alignments only
            # Keeps: primary + supplementary (for visualization)
            cmd1 = ['samtools', 'view', '-F', '256', '-b', '-o', temp_bam, sam_file]
        else:
            cmd1 = ['samtools', 'view', '-b', '-o', temp_bam, sam_file]
        
        result1 = subprocess.run(cmd1, check=True, capture_output=True, text=True)
        if result1.stderr:
            logging.debug(f"samtools view output: {result1.stderr}")
        
        # Step 2: Sort the BAM
        cmd2 = ['samtools', 'sort', '-o', output_bam, temp_bam]
        result2 = subprocess.run(cmd2, check=True, capture_output=True, text=True)
        if result2.stderr:
            logging.debug(f"samtools sort output: {result2.stderr}")
        
        # Step 3: Create index
        cmd3 = ['samtools', 'index', output_bam]
        result3 = subprocess.run(cmd3, check=True, capture_output=True, text=True)
        if result3.stderr:
            logging.debug(f"samtools index output: {result3.stderr}")
        
        logging.info(f"Created {os.path.basename(output_bam)} and index")
        
    except subprocess.CalledProcessError as e:
        logging.error(f"samtools command failed: {' '.join(e.cmd)}")
        logging.error(f"Return code: {e.returncode}")
        if e.stderr:
            logging.error(f"Error output: {e.stderr}")
        if e.stdout:
            logging.error(f"Standard output: {e.stdout}")
        raise
    
    finally:
        # Clean up temporary file
        if os.path.exists(temp_bam):
            os.remove(temp_bam)
    
    return output_bam


def extract_mapped_reads(sam_file: str, output_fastq: str) -> int:
    """
    Extract primary mapped reads from SAM file to FASTQ format using samtools.
    Excludes unmapped, secondary, and supplementary alignments.
    
    Args:
        sam_file: Input SAM file
        output_fastq: Output FASTQ file path
        
    Returns:
        Number of reads extracted
    """
    logging.info("Extracting primary mapped reads to FASTQ")
    
    # Use samtools fastq on primary mapped reads only
    # -F 2308 excludes: unmapped (4) + secondary (256) + supplementary (2048)
    cmd = ['samtools', 'fastq', '-F', '2308', sam_file, '-0', output_fastq]
    
    try:
        result = subprocess.run(
            cmd, 
            check=True, 
            capture_output=True,
            text=True
        )
        if result.stderr:
            logging.debug(f"samtools output: {result.stderr}")
    except subprocess.CalledProcessError as e:
        logging.error(f"samtools command failed: {' '.join(cmd)}")
        logging.error(f"Error output: {e.stderr}")
        raise
    
    # Count reads using seqkit (more reliable)
    if os.path.exists(output_fastq) and os.path.getsize(output_fastq) > 0:
        try:
            result = subprocess.run(
                ['seqkit', 'stats', '-T', output_fastq],
                capture_output=True,
                text=True,
                check=True
            )
            # Parse seqkit output (second line, num_seqs column)
            lines = result.stdout.strip().split('\n')
            if len(lines) > 1:
                fields = lines[1].split('\t')
                read_count = int(fields[3])  # num_seqs is 4th column
            else:
                read_count = 0
        except:
            # Fallback to simple counting
            with open(output_fastq) as f:
                read_count = sum(1 for line in f if line.startswith('@'))
    else:
        read_count = 0
        logging.warning(f"Output FASTQ is empty or not created: {output_fastq}")
    
    logging.info(f"Extracted {read_count} primary mapped reads")
    return read_count


def calculate_alignment_coverage(read) -> float:
    """
    Calculate alignment coverage from a pysam AlignedSegment.
    
    Coverage = (aligned_length_on_read / total_read_length)
    
    Aligned length = read length - soft clipped bases
    
    Args:
        read: pysam AlignedSegment object
        
    Returns:
        Coverage as a float (0.0 to 1.0)
    """
    if read.query_length is None or read.query_length == 0:
        return 0.0
    
    # Get total read length
    read_length = read.query_length
    
    # Count all soft-clipped bases from CIGAR
    soft_clipped = 0
    if read.cigartuples:
        for op, length in read.cigartuples:
            if op == 4:  # 4 = soft clip (S in CIGAR)
                soft_clipped += length
    
    # Calculate aligned length
    aligned_length = read_length - soft_clipped
    
    # Calculate coverage
    coverage = aligned_length / read_length
    
    return min(coverage, 1.0)  # Cap at 1.0


def filter_reads_by_coverage(sam_file: str, output_fastq: str, 
                            min_coverage: float = 0.9) -> int:
    """
    Extract high-quality reads based on alignment coverage threshold.
    Mimics sam2fastq behavior: processes primary + supplementary,
    but outputs only FIRST occurrence of each read name.
    
    Coverage = aligned_length / read_length
    Aligned length = read_length - all_soft_clipped_bases
    
    Args:
        sam_file: Input SAM file
        output_fastq: Output FASTQ file path
        min_coverage: Minimum alignment coverage (0.0 to 1.0)
        
    Returns:
        Number of reads passing filter
    """
    logging.info(f"Filtering reads with >={min_coverage*100:.0f}% coverage")
    
    samfile = pysam.AlignmentFile(sam_file, "r")
    read_count = 0
    seen_reads = {}  # Track which reads we've already written
    
    with open(output_fastq, 'w') as outfile:
        for read in samfile:
            # Skip unmapped reads
            if read.is_unmapped:
                continue
            
            # Skip secondary alignments (mimics sam2fastq checking tp:A:S)
            if read.is_secondary:
                continue
            
            # Supplementary alignments have tp:A:P, so they pass through
            # (just like in sam2fastq)
            
            # Check if we've already written this read name
            if read.query_name in seen_reads:
                continue  # Skip - already written (sam2fastq behavior)
            
            # Skip reads without sequence or quality
            if read.query_sequence is None or read.query_qualities is None:
                continue
            
            # Calculate coverage
            try:
                coverage = calculate_alignment_coverage(read)
            except Exception as e:
                logging.debug(f"Could not calculate coverage for {read.query_name}: {e}")
                continue
            
            # Write if passes threshold
            if coverage >= min_coverage:
                try:
                    quality_string = ''.join(chr(q + 33) for q in read.query_qualities)
                    outfile.write(f"@{read.query_name}\n")
                    outfile.write(f"{read.query_sequence}\n")
                    outfile.write("+\n")
                    outfile.write(f"{quality_string}\n")
                    seen_reads[read.query_name] = True  # Mark as written
                    read_count += 1
                except Exception as e:
                    logging.debug(f"Could not write read {read.query_name}: {e}")
                    continue
    
    samfile.close()
    
    logging.info(f"Extracted {read_count} high-coverage reads (unique reads, first occurrence)")
    return read_count


def create_fastq_difference(fastq_all: str, fastq_subset: str, 
                           output_fastq: str) -> int:
    """
    Create FASTQ file containing reads in fastq_all but not in fastq_subset.
    Uses seqkit for efficiency.
    
    Args:
        fastq_all: FASTQ file with all reads
        fastq_subset: FASTQ file with reads to exclude
        output_fastq: Output FASTQ file path
        
    Returns:
        Number of reads in output (validated count)
    """
    logging.info(f"Creating difference FASTQ: {os.path.basename(output_fastq)}")
    
    # Extract read names from subset
    cmd1 = f"seqkit seq -i --name {fastq_subset}"
    
    # Grep inverse to get reads NOT in subset
    cmd2 = f"seqkit grep -v -f - {fastq_all}"
    
    # Full pipeline
    cmd = f"{cmd1} | {cmd2} > {output_fastq}"
    
    subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE)
    
    # Count reads properly using seqkit (not grep!)
    from .utils import count_fastq_reads
    read_count = count_fastq_reads(output_fastq)
    
    logging.info(f"Created difference FASTQ with {read_count} reads")
    return read_count


def create_bam_for_mapped_reads(sam_file: str, output_bam: str) -> str:
    """
    Create sorted BAM file containing primary and supplementary mapped reads.
    Excludes unmapped and secondary alignments.
    
    Supplementary alignments are kept for visualization (IGV) but not used in analysis.
    
    Args:
        sam_file: Input SAM file
        output_bam: Output BAM file path
        
    Returns:
        Path to output BAM file
    """
    logging.info(f"Creating BAM for mapped reads (primary + supplementary)")
    
    # Create temporary unsorted BAM
    temp_bam = output_bam + ".tmp.bam"
    
    try:
        # Step 1: Convert SAM to BAM, filtering out unmapped (4) and secondary (256)
        # -F 260 excludes: unmapped (4) + secondary (256)
        # Keeps: primary + supplementary
        cmd1 = ['samtools', 'view', '-F', '260', '-b', '-o', temp_bam, sam_file]
        result1 = subprocess.run(cmd1, check=True, capture_output=True, text=True)
        if result1.stderr:
            logging.debug(f"samtools view output: {result1.stderr}")
        
        # Step 2: Sort the BAM
        cmd2 = ['samtools', 'sort', '-o', output_bam, temp_bam]
        result2 = subprocess.run(cmd2, check=True, capture_output=True, text=True)
        if result2.stderr:
            logging.debug(f"samtools sort output: {result2.stderr}")
        
        # Step 3: Create index
        cmd3 = ['samtools', 'index', output_bam]
        result3 = subprocess.run(cmd3, check=True, capture_output=True, text=True)
        if result3.stderr:
            logging.debug(f"samtools index output: {result3.stderr}")
        
        logging.info(f"Created {os.path.basename(output_bam)} and index")
        
    except subprocess.CalledProcessError as e:
        logging.error(f"samtools command failed: {' '.join(e.cmd)}")
        logging.error(f"Return code: {e.returncode}")
        if e.stderr:
            logging.error(f"Error output: {e.stderr}")
        if e.stdout:
            logging.error(f"Standard output: {e.stdout}")
        raise
    
    finally:
        # Clean up temporary file
        if os.path.exists(temp_bam):
            os.remove(temp_bam)
    
    return output_bam


def align_to_multiple_references(fastq_file: str, references: dict, 
                                 output_prefix: str, threads: int = 4) -> dict:
    """
    Align reads to multiple reference sequences.
    
    Args:
        fastq_file: Input FASTQ file
        references: Dict mapping reference names to FASTA paths
        output_prefix: Prefix for output files
        threads: Number of threads
        
    Returns:
        Dict mapping reference names to BAM file paths
    """
    bam_files = {}
    
    for ref_name, ref_path in references.items():
        sam_file = f"{output_prefix}_{ref_name}.sam"
        bam_file = f"{output_prefix}_{ref_name}.bam"
        
        # Align
        align_reads(fastq_file, ref_path, sam_file, threads)
        
        # Convert to BAM and sort
        sam_to_sorted_bam(sam_file, bam_file)
        
        # Clean up SAM
        os.remove(sam_file)
        
        bam_files[ref_name] = bam_file
    
    return bam_files