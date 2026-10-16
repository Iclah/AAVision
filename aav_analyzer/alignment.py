"""
Alignment and read filtering functions for AAV analysis.
Handles minimap2 alignment and quality-based filtering.

minimap2 runs with -Y, so supplementary records are soft-clipped and carry
the full read sequence.

Which records each step uses:
- vector.fastq (samtools fastq -F 2308): primary alignments only, written in
  the read's original orientation.
- HQ/LQ split (filter_reads_by_coverage): primary alignments only.
- vector.bam (-F 260): primary + supplementary, no secondary. Used for IGV,
  for the samtools-depth coverage plot (which counts supplementary records),
  and for truncation hotspots (which skip supplementary records).
- Hybrid/impurity BAMs (sam_to_sorted_bam, -F 256): primary + supplementary
  + unmapped. Best-reference assignment takes the highest AS over a read's
  primary and supplementary records.

samtools flags used:
- -F 260:  exclude unmapped (4) + secondary (256)
- -F 2308: exclude unmapped (4) + secondary (256) + supplementary (2048)
- -F 256:  exclude secondary (256) only
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
    # -Y: soft-clip supplementary alignments (minimap2 hard-clips them by
    # default). Supplementary records then carry the full read SEQ and length,
    # so IGV shows full read context and CIGAR-based read lengths are correct
    # for every record. HQ/LQ does not depend on this: it uses primaries only.
    cmd = [
        'minimap2',
        '-ax', 'map-ont',
        '-Y',
        '-t', str(threads),
        reference_fasta,
        fastq_file
    ]
    
    logging.info(f"Aligning reads to {os.path.basename(reference_fasta)}: {' '.join(cmd)}")
    
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
    Fraction of the full read that one alignment record covers.

    Coverage = read bases in M/I/=/X operations / full read length, where the
    full read length is taken from the CIGAR including hard clips
    (pysam infer_read_length). Including hard clips matters: minimap2 without
    -Y hard-clips supplementary records, so query_length alone would make
    every supplementary record look like it covers 100% of the read.

    Deletions (D/N) consume reference only and do not lower coverage: a read
    with a large internal deletion bridged by one alignment scores 1.0. This
    is intended (docs/DECISIONS.md). No positions are involved.

    Args:
        read: pysam AlignedSegment (mapped, with a CIGAR)

    Returns:
        Coverage in [0.0, 1.0]; 0.0 if the read has no CIGAR
    """
    read_length = read.infer_read_length()
    if not read_length or not read.cigartuples:
        return 0.0

    # 0=M, 1=I, 7='=', 8=X consume query bases that are aligned
    aligned_length = sum(length for op, length in read.cigartuples if op in (0, 1, 7, 8))
    return aligned_length / read_length


def filter_reads_by_coverage(sam_file: str, vector_fastq: str, output_fastq: str,
                             min_coverage: float = 0.9) -> int:
    """
    Write high-quality (HQ) vector reads: reads whose PRIMARY alignment
    covers >= min_coverage of the full read.

    The HQ decision uses the primary record only, once per read. Secondary
    and supplementary records are ignored, so a split (chimeric) read is
    judged on its primary alignment whatever order records appear in.

    HQ reads are copied from vector_fastq (as written by extract_mapped_reads),
    not rebuilt from SAM records, so they keep their original orientation,
    full length and qualities. Minus-strand SAM records store the reverse
    complement, and hard-clipped records store only part of the read.

    Args:
        sam_file: Vector alignment SAM/BAM (any record order)
        vector_fastq: FASTQ of all primary-mapped vector reads, original orientation
        output_fastq: Output FASTQ for HQ reads
        min_coverage: Minimum primary-alignment coverage, inclusive (0.0 to 1.0)

    Returns:
        Number of HQ reads written

    Raises:
        ValueError: if an HQ read from the SAM is missing from vector_fastq
            (inputs out of sync)
    """
    logging.info(f"Filtering reads with >={min_coverage*100:.0f}% primary-alignment coverage")

    hq_names = set()
    n_primary = 0
    n_lq_with_supp = 0  # LQ reads with an SA tag: promoted to HQ before this fix
    n_skipped = {'unmapped': 0, 'secondary': 0, 'supplementary': 0}

    with pysam.AlignmentFile(sam_file, "r") as samfile:
        for read in samfile:
            if read.is_unmapped:
                n_skipped['unmapped'] += 1
                continue
            if read.is_secondary:
                n_skipped['secondary'] += 1
                continue
            if read.is_supplementary:
                n_skipped['supplementary'] += 1
                continue

            n_primary += 1
            if calculate_alignment_coverage(read) >= min_coverage:
                hq_names.add(read.query_name)
            elif read.has_tag('SA'):
                n_lq_with_supp += 1

    read_count = 0
    with open(output_fastq, 'w') as outfile:
        if os.path.exists(vector_fastq) and os.path.getsize(vector_fastq) > 0:
            with pysam.FastxFile(vector_fastq) as fq:
                for entry in fq:
                    if entry.name in hq_names:
                        outfile.write(f"@{entry.name}\n{entry.sequence}\n+\n{entry.quality}\n")
                        read_count += 1

    if read_count != len(hq_names):
        raise ValueError(
            f"{len(hq_names) - read_count} HQ reads from {sam_file} are missing "
            f"from {vector_fastq}; were both produced from the same alignment?"
        )

    logging.info(
        f"Evaluated {n_primary} primary alignments: {read_count} HQ, "
        f"{n_primary - read_count} below threshold "
        f"({n_lq_with_supp} of these have supplementary alignments)"
    )
    logging.info(
        "Records not used for the HQ decision: "
        + ", ".join(f"{k}={v}" for k, v in n_skipped.items())
    )
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