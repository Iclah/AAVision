"""
Truncation hotspot analysis for AAV reads.
Identifies positions with enriched read termini relative to coverage.
"""

import os
import logging
import pysam
import pandas as pd
from typing import Dict


def compute_truncation_hotspots(bam_file: str, smoothing_window: int = 10, 
                               max_pos: int = 6000) -> pd.DataFrame:
    """
    Analyze BAM file to identify truncation hotspots.
    
    Counts read termini (start + end positions), normalizes by coverage,
    applies rolling smoothing, and converts to percentage enrichment scores.
    
    Args:
        bam_file: Path to sorted BAM file
        smoothing_window: Window size for rolling average smoothing
        max_pos: Maximum genomic position to analyze
        
    Returns:
        DataFrame with columns: chrom, pos, termini, coverage, enrichment, smoothed
    """
    logging.info(f"Computing truncation hotspots from {bam_file}")
    
    bam = pysam.AlignmentFile(bam_file, "rb")
    frames = []
    
    for chrom, chrom_len in zip(bam.references, bam.lengths):
        # Initialize arrays for this chromosome
        termini_counts = [0] * chrom_len
        coverage = [0] * chrom_len
        
        # Count termini and coverage
        for read in bam.fetch(chrom):
            if read.is_unmapped:
                continue
                
            start = read.reference_start
            end = read.reference_end - 1
            
            # Count read termini
            if 0 <= start < chrom_len:
                termini_counts[start] += 1
            if 0 <= end < chrom_len:
                termini_counts[end] += 1
            
            # Count coverage at each position
            for i in range(start, min(end + 1, chrom_len)):
                coverage[i] += 1
        
        # Create DataFrame for this chromosome
        df = pd.DataFrame({
            "chrom": chrom,
            "pos": range(1, chrom_len + 1),
            "termini": termini_counts,
            "coverage": coverage
        })
        
        # Normalize by coverage (enrichment score)
        df["enrichment"] = df.apply(
            lambda row: row["termini"] / row["coverage"] if row["coverage"] > 0 else 0,
            axis=1
        )
        
        # Apply rolling smoothing
        df["smoothed"] = df["enrichment"].rolling(
            window=smoothing_window,
            center=True,
            min_periods=1
        ).mean()
        
        # Clip to max_pos
        df = df[df["pos"] <= max_pos].copy()
        
        # Convert to percentage
        df["smoothed"] = df["smoothed"] * 100
        
        frames.append(df)
    
    bam.close()
    
    result = pd.concat(frames, ignore_index=True)
    logging.info(f"Truncation analysis complete: {len(result)} positions analyzed")
    
    return result


def get_truncation_data_dict(df: pd.DataFrame) -> Dict[int, float]:
    """
    Convert truncation DataFrame to dictionary for plotting.
    
    Args:
        df: DataFrame from compute_truncation_hotspots
        
    Returns:
        Dict mapping position to smoothed percentage
    """
    data = {}
    for _, row in df.iterrows():
        pos = int(row['pos'])
        percentage = float(row['smoothed'])
        if percentage > 0:  # Only include non-zero values
            data[pos] = percentage
    
    return data


def compute_coverage(bam_file: str, max_pos: int = 6000) -> Dict[int, int]:
    """
    Calculate per-position coverage from BAM file using samtools depth.
    This matches IGV's coverage calculation method.

    Args:
        bam_file: Path to sorted BAM file with AAV reads
        max_pos: Maximum genomic position to analyze

    Returns:
        Dict mapping genomic position to coverage depth
    """
    import subprocess

    logging.info(f"Computing coverage from {bam_file}")

    coverage_data = {}

    try:
        # Use samtools depth with flags that match IGV behavior:
        # -a: output all positions (including zero coverage)
        # -d 0: no max depth limit
        # -Q 0: no base quality filtering
        # -q 0: no mapping quality filtering
        cmd = ['samtools', 'depth', '-a', '-d', '0', '-Q', '0', '-q', '0', bam_file]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        # Parse output: format is "chrom\tpos\tdepth"
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) >= 3:
                pos = int(parts[1])
                depth = int(parts[2])

                if pos <= max_pos:
                    coverage_data[pos] = depth

    except subprocess.CalledProcessError as e:
        logging.error(f"samtools depth failed: {e.stderr}")
        raise
    except Exception as e:
        logging.error(f"Error computing coverage: {e}")
        raise

    logging.info(f"Coverage computed for {len(coverage_data)} positions")
    return coverage_data


def find_hotspot_peaks(df: pd.DataFrame, threshold: float = 2.0) -> pd.DataFrame:
    """
    Identify peak positions in truncation hotspots.
    
    Args:
        df: DataFrame from compute_truncation_hotspots
        threshold: Minimum enrichment percentage to consider as hotspot
        
    Returns:
        DataFrame with hotspot peak positions and their enrichment scores
    """
    # Filter positions above threshold
    hotspots = df[df["smoothed"] > threshold].copy()
    
    if len(hotspots) == 0:
        return pd.DataFrame(columns=["chrom", "pos", "enrichment"])
    
    # Find local maxima
    hotspots["is_peak"] = False
    
    for idx in hotspots.index:
        curr_val = hotspots.loc[idx, "smoothed"]
        curr_pos = hotspots.loc[idx, "pos"]
        
        # Check if this is a local maximum
        neighbors = df[
            (df["pos"] >= curr_pos - 5) & 
            (df["pos"] <= curr_pos + 5)
        ]
        
        if curr_val == neighbors["smoothed"].max():
            hotspots.loc[idx, "is_peak"] = True
    
    peaks = hotspots[hotspots["is_peak"]][["chrom", "pos", "smoothed"]].copy()
    peaks.rename(columns={"smoothed": "enrichment"}, inplace=True)
    
    logging.info(f"Found {len(peaks)} truncation hotspot peaks")
    
    return peaks


def analyze_read_lengths(fastq_files: list, bin_size: int = 100, 
                        max_length: int = 6000) -> Dict[str, Dict[str, float]]:
    """
    Analyze read length distributions from multiple FASTQ files.
    
    Args:
        fastq_files: List of tuples (name, filepath) for FASTQ files to analyze
        bin_size: Bin size for histogram (default 100bp)
        max_length: Maximum length to bin separately (longer reads go to "≥max_length")
        
    Returns:
        Dict mapping file names to binned length distributions (as percentages)
        All bins from 0 to max_length are included (even if empty) for consistent axes
    """
    logging.info("Analyzing read length distributions")
    
    results = {}
    
    # Pre-generate all bin keys for consistent axes across samples
    all_bin_keys = []
    for bin_start in range(0, max_length, bin_size):
        bin_end = bin_start + bin_size - 1
        all_bin_keys.append(f"{bin_start}-{bin_end}")
    overflow_key = f"≥{max_length}"
    all_bin_keys.append(overflow_key)
    
    for name, filepath in fastq_files:
        # Skip if file doesn't exist or is empty
        if not os.path.exists(filepath):
            logging.warning(f"File not found, skipping: {filepath}")
            continue
        
        if os.path.getsize(filepath) == 0:
            logging.warning(f"File is empty, skipping: {filepath}")
            continue
        
        # Initialize all bins to zero for consistent axes
        bins = {key: 0 for key in all_bin_keys}
        total_reads = 0
        
        try:
            with pysam.FastxFile(filepath) as fq:
                for entry in fq:
                    length = len(entry.sequence)
                    total_reads += 1
                    
                    # Determine bin
                    if length >= max_length:
                        bin_key = overflow_key
                    else:
                        bin_start = (length // bin_size) * bin_size
                        bin_end = bin_start + bin_size - 1
                        bin_key = f"{bin_start}-{bin_end}"
                    
                    bins[bin_key] = bins.get(bin_key, 0) + 1
            
            # Convert to percentages
            if total_reads > 0:
                bins = {k: (v / total_reads) * 100 for k, v in bins.items()}
            
            results[name] = bins
            logging.debug(f"Analyzed {total_reads} reads from {name}")
        except Exception as e:
            logging.warning(f"Error analyzing {filepath}: {e}")
            continue
    
    return results