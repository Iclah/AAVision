"""
Taxonomic classification of unclassified reads using Kraken2 and Bracken.

Provides an optional pipeline step to identify what unclassified reads are
by running Kraken2 for k-mer-based taxonomic assignment, followed by Bracken
for abundance re-estimation at a specified taxonomic level.

Requires:
    - kraken2 (on PATH)
    - bracken (on PATH, optional — skipped gracefully if absent)
    - A pre-built Kraken2 database (e.g. k2_standard_08gb)
"""

import os
import logging
import shutil
import subprocess
import csv
from typing import Dict, Optional, Tuple


def check_kraken_available() -> bool:
    """Check if kraken2 is installed and on PATH."""
    return shutil.which('kraken2') is not None


def check_bracken_available() -> bool:
    """Check if bracken is installed and on PATH."""
    return shutil.which('bracken') is not None


def validate_kraken_db(db_path: str) -> bool:
    """
    Validate that a Kraken2 database directory contains required files.
    
    Args:
        db_path: Path to Kraken2 database directory
        
    Returns:
        True if database looks valid
    """
    required_files = ['hash.k2d', 'opts.k2d', 'taxo.k2d']
    
    for fname in required_files:
        if not os.path.exists(os.path.join(db_path, fname)):
            logging.error(f"Kraken2 database missing required file: {fname}")
            logging.error(f"Database path: {db_path}")
            return False
    
    return True


def run_kraken2(fastq_file: str, db_path: str, output_dir: str,
                threads: int = 4, confidence: float = 0.05) -> Dict:
    """
    Run Kraken2 classification on a FASTQ file.
    
    Args:
        fastq_file: Input FASTQ file (unclassified reads)
        db_path: Path to Kraken2 database
        output_dir: Directory for output files
        threads: Number of threads
        confidence: Minimum confidence score for classification (0.0-1.0)
        
    Returns:
        Dict with paths to output files and summary stats
    """
    logging.info(f"Running Kraken2 classification (confidence={confidence})")
    
    output_file = os.path.join(output_dir, 'kraken2_output.txt')
    report_file = os.path.join(output_dir, 'kraken2_report.txt')
    
    cmd = [
        'kraken2',
        '--db', db_path,
        '--threads', str(threads),
        '--confidence', str(confidence),
        '--output', output_file,
        '--report', report_file,
        fastq_file
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse stderr for summary stats
        stats = _parse_kraken2_stderr(result.stderr)
        
        logging.info(f"Kraken2 complete: {stats.get('classified', 0)} classified, "
                     f"{stats.get('unclassified', 0)} unclassified")
        
        return {
            'output_file': output_file,
            'report_file': report_file,
            'stats': stats
        }
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Kraken2 failed: {e.stderr}")
        raise


def _parse_kraken2_stderr(stderr: str) -> Dict:
    """
    Parse Kraken2 stderr output to extract classification summary.
    
    Typical output:
        12345 sequences classified (85.23%)
        2156 sequences unclassified (14.77%)
    """
    stats = {
        'total': 0,
        'classified': 0,
        'unclassified': 0,
        'classified_percent': 0.0,
        'unclassified_percent': 0.0
    }
    
    for line in stderr.strip().split('\n'):
        line = line.strip()
        if 'sequences classified' in line:
            parts = line.split()
            stats['classified'] = int(parts[0])
            # Extract percentage from parentheses
            for p in parts:
                if '%' in p:
                    stats['classified_percent'] = float(p.strip('()%'))
        elif 'sequences unclassified' in line or 'unclassified' in line:
            parts = line.split()
            stats['unclassified'] = int(parts[0])
            for p in parts:
                if '%' in p:
                    stats['unclassified_percent'] = float(p.strip('()%'))
        elif 'processed' in line.lower() or 'sequences' in line.lower():
            # Try to get total from first number
            parts = line.split()
            try:
                stats['total'] = int(parts[0])
            except (ValueError, IndexError):
                pass
    
    if stats['total'] == 0:
        stats['total'] = stats['classified'] + stats['unclassified']
    
    return stats


def run_bracken(kraken_report: str, db_path: str, output_dir: str,
                read_length: Optional[int] = None, level: str = 'S',
                threshold: int = 10) -> Optional[str]:
    """
    Run Bracken abundance re-estimation on Kraken2 output.

    Args:
        kraken_report: Path to Kraken2 report file
        db_path: Path to Kraken2 database (containing Bracken DB files)
        output_dir: Directory for output files
        read_length: Bracken read length (-r). This should reflect the read
                     length distribution of the reads being classified, not the
                     database maximum. If None, the largest length the database
                     was built for is used as a fallback. If given, it must match
                     one of the lengths the database was built for; otherwise this
                     is an error (no silent substitution), because -r changes the
                     abundance estimates and must be reproducible.
        level: Taxonomic level for estimation:
               D=Domain, P=Phylum, C=Class, O=Order, F=Family, G=Genus, S=Species
        threshold: Minimum number of reads required for a taxon

    Returns:
        Path to Bracken output file, or None if Bracken unavailable/failed
    """
    if not check_bracken_available():
        logging.warning("Bracken not found on PATH — skipping abundance re-estimation")
        return None
    
    # Check for Bracken database files
    bracken_db_files = [f for f in os.listdir(db_path) 
                        if f.startswith('database') and f.endswith('.kmer_distrib')]
    
    if not bracken_db_files:
        logging.warning(f"No Bracken database files found in {db_path} — skipping")
        return None
    
    # Parse the read lengths the database was built for
    available_lengths = []
    for f in bracken_db_files:
        try:
            # Format: database{LENGTH}mers.kmer_distrib
            length_str = f.replace('database', '').replace('mers.kmer_distrib', '')
            available_lengths.append(int(length_str))
        except ValueError:
            continue
    
    if not available_lengths:
        logging.warning("Could not parse Bracken database read lengths — skipping")
        return None

    available_sorted = sorted(available_lengths)

    if read_length is not None:
        # Honor the explicit choice, but require an exact match to a length the
        # database supports. Do NOT silently substitute — a wrong -r changes the
        # estimates, so surface the mismatch loudly instead.
        if read_length not in available_lengths:
            logging.error(
                f"Requested Bracken read length {read_length} is not available in "
                f"this database. Available lengths: {available_sorted}. "
                f"Rebuild the Bracken DB for {read_length} or choose an available length."
            )
            raise ValueError(
                f"Bracken read length {read_length} not available (have {available_sorted})"
            )
        best_length = read_length
        logging.info(f"Using requested Bracken read length: {best_length} "
                     f"(available: {available_sorted})")
    else:
        # Fallback when no length is specified: use the largest available.
        best_length = max(available_lengths)
        logging.info(f"No Bracken read length specified; using database maximum "
                     f"{best_length} (available: {available_sorted})")
    
    bracken_output = os.path.join(output_dir, f'bracken_{level}.txt')
    bracken_report = os.path.join(output_dir, f'bracken_{level}_report.txt')
    
    bracken_output = os.path.join(output_dir, f'bracken_{level}.txt')
    bracken_report = os.path.join(output_dir, f'bracken_{level}_report.txt')
    
    cmd = [
        'bracken',
        '-d', db_path,
        '-i', kraken_report,
        '-o', bracken_output,
        '-w', bracken_report,
        '-r', str(best_length),
        '-l', level,
        '-t', str(threshold)
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        logging.info(f"Bracken complete: {bracken_output}")
        return bracken_output
        
    except subprocess.CalledProcessError as e:
        logging.warning(f"Bracken failed (non-fatal): {e.stderr}")
        return None


def parse_kraken_report(report_file: str, min_percent: float = 0.5,
                        max_entries: int = 20) -> list:
    """
    Parse a Kraken2 standard report into a structured list.
    
    Report format (tab-separated):
        percent  reads_clade  reads_direct  rank  taxid  name
    
    Args:
        report_file: Path to Kraken2 report file
        min_percent: Minimum percentage to include in results
        max_entries: Maximum number of entries to return
        
    Returns:
        List of dicts with keys: percent, reads_clade, reads_direct, rank, taxid, name
    """
    results = []
    
    try:
        with open(report_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                
                percent = float(parts[0].strip())
                reads_clade = int(parts[1].strip())
                reads_direct = int(parts[2].strip())
                rank = parts[3].strip()
                taxid = parts[4].strip()
                name = parts[5].strip()
                
                # Skip root and unclassified for the detailed breakdown
                if rank in ('R', 'R1') or name == 'root':
                    continue
                
                if percent >= min_percent:
                    results.append({
                        'percent': percent,
                        'reads_clade': reads_clade,
                        'reads_direct': reads_direct,
                        'rank': rank,
                        'taxid': taxid,
                        'name': name
                    })
    except Exception as e:
        logging.error(f"Error parsing Kraken2 report: {e}")
        return []
    
    # Sort by percentage descending
    results.sort(key=lambda x: x['percent'], reverse=True)
    
    return results[:max_entries]


def parse_bracken_output(bracken_file: str, min_percent: float = 0.5,
                          max_entries: int = 20) -> list:
    """
    Parse Bracken abundance output file.
    
    Bracken output is tab-separated with headers:
        name  taxonomy_id  taxonomy_lvl  kraken_assigned  added  new_est  fraction
    
    Args:
        bracken_file: Path to Bracken output file
        min_percent: Minimum fraction (as percent) to include
        max_entries: Maximum entries to return
        
    Returns:
        List of dicts with keys: name, taxid, level, reads, fraction_percent
    """
    results = []
    
    try:
        with open(bracken_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                fraction = float(row.get('fraction_total_reads', 0)) * 100
                if fraction >= min_percent:
                    results.append({
                        'name': row.get('name', 'Unknown'),
                        'taxid': row.get('taxonomy_id', ''),
                        'level': row.get('taxonomy_lvl', ''),
                        'reads': int(row.get('new_est_reads', 0)),
                        'fraction_percent': round(fraction, 2)
                    })
    except Exception as e:
        logging.error(f"Error parsing Bracken output: {e}")
        return []
    
    results.sort(key=lambda x: x['fraction_percent'], reverse=True)
    return results[:max_entries]


def run_taxonomy_analysis(fastq_file: str, db_path: str, output_dir: str,
                          threads: int = 4, confidence: float = 0.05,
                          bracken_level: str = 'S',
                          bracken_threshold: int = 10,
                          bracken_read_length: Optional[int] = None) -> Dict:
    """
    Run full taxonomic analysis pipeline: Kraken2 + optional Bracken.
    
    Args:
        fastq_file: Input FASTQ with unclassified reads
        db_path: Path to Kraken2/Bracken database
        output_dir: Output directory for taxonomy results
        threads: Number of threads
        confidence: Kraken2 confidence threshold
        bracken_level: Taxonomic level for Bracken (S=Species, G=Genus, etc.)
        bracken_threshold: Minimum reads for Bracken taxon
        bracken_read_length: Bracken read length (-r) to match the data's read
                             length distribution. If None, the database maximum
                             is used as a fallback.
        
    Returns:
        Dict with all taxonomy results for report generation:
        {
            'kraken_stats': {...},
            'kraken_report': [...],     # Parsed Kraken report entries
            'bracken_results': [...],   # Parsed Bracken results (or None)
            'bracken_level': str,       # Level used
            'domain_summary': {...},    # High-level domain breakdown
        }
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Run Kraken2
    kraken_result = run_kraken2(fastq_file, db_path, output_dir, threads, confidence)
    
    # Parse Kraken2 report
    kraken_entries = parse_kraken_report(kraken_result['report_file'])
    
    # Build domain-level summary from Kraken report
    domain_summary = _build_domain_summary(kraken_result['report_file'])
    
    # Run Bracken
    bracken_results = None
    bracken_file = run_bracken(
        kraken_result['report_file'], db_path, output_dir,
        read_length=bracken_read_length,
        level=bracken_level, threshold=bracken_threshold
    )
    
    if bracken_file and os.path.exists(bracken_file):
        bracken_results = parse_bracken_output(bracken_file)
    
    return {
        'kraken_stats': kraken_result['stats'],
        'kraken_report': kraken_entries,
        'bracken_results': bracken_results,
        'bracken_level': bracken_level,
        'domain_summary': domain_summary
    }


def _build_domain_summary(report_file: str) -> Dict:
    """
    Build a high-level domain breakdown from Kraken2 report.
    Extracts top-level domains: Bacteria, Viruses, Eukaryota, Archaea.
    
    Returns:
        Dict mapping domain names to read counts and percentages
    """
    domains = {}
    unclassified_reads = 0
    total_reads = 0
    
    try:
        with open(report_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                
                percent = float(parts[0].strip())
                reads_clade = int(parts[1].strip())
                rank = parts[3].strip()
                name = parts[5].strip()
                
                # Track unclassified
                if name == 'unclassified':
                    unclassified_reads = reads_clade
                    total_reads += reads_clade
                elif rank == 'R' and name == 'root':
                    total_reads += reads_clade
                
                # Domain level
                if rank == 'D':
                    domains[name] = {
                        'reads': reads_clade,
                        'percent': percent
                    }
    except Exception as e:
        logging.error(f"Error building domain summary: {e}")
    
    if total_reads > 0 and unclassified_reads > 0:
        domains['Unclassified'] = {
            'reads': unclassified_reads,
            'percent': round((unclassified_reads / total_reads) * 100, 2)
        }
    
    return domains
