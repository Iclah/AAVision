"""
Utility functions for AAV analysis pipeline.
Handles validation, logging, file operations, and tool checks.
"""

import os
import sys
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Optional, List, Dict


def setup_logging(verbose: bool = False, log_file: str = None) -> logging.Logger:
    """
    Configure logging for the pipeline.
    
    Args:
        verbose: If True, show DEBUG messages
        log_file: Optional path to log file
        
    Returns:
        Configured logger instance
    """
    level = logging.DEBUG if verbose else logging.INFO
    
    # Create formatters
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Setup root logger
    logger = logging.getLogger()
    logger.setLevel(level)
    
    # Remove existing handlers
    logger.handlers = []
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)  # Always save DEBUG to file
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logging.getLogger('aav_analyzer')


def check_required_tools() -> bool:
    """
    Check if all required bioinformatics tools are installed.
    
    Returns:
        True if all tools found, False otherwise
    """
    required_tools = ['minimap2', 'samtools', 'seqkit']
    missing_tools = []
    
    for tool in required_tools:
        if not shutil.which(tool):
            missing_tools.append(tool)
    
    if missing_tools:
        logging.error("Missing required tools: %s", ', '.join(missing_tools))
        logging.error("Install via: conda install -c bioconda %s", ' '.join(missing_tools))
        return False
    
    return True


def find_reference_file(base_name: str, search_dir: str = '.') -> Optional[str]:
    """
    Find reference file with .fasta or .fa extension.
    
    Args:
        base_name: Base filename without extension (e.g., 'vector')
        search_dir: Directory to search in
        
    Returns:
        Full path to reference file, or None if not found
    """
    for ext in ['.fasta', '.fa']:
        path = os.path.join(search_dir, f"{base_name}{ext}")
        if os.path.exists(path):
            return path
    return None


def validate_fasta(filepath: str) -> bool:
    """
    Basic validation that a file is in FASTA format.
    
    Args:
        filepath: Path to FASTA file
        
    Returns:
        True if valid FASTA format
    """
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
            if not first_line.startswith('>'):
                logging.error(f"Invalid FASTA format in {filepath}: missing '>' header")
                return False
        return True
    except Exception as e:
        logging.error(f"Error reading {filepath}: {e}")
        return False


def validate_fastq(filepath: str) -> bool:
    """
    Basic validation that a file is in FASTQ format.
    Handles both .fastq and .fastq.gz files.
    
    Args:
        filepath: Path to FASTQ file
        
    Returns:
        True if valid FASTQ format
    """
    import gzip
    
    try:
        # Determine if file is gzipped
        if filepath.endswith('.gz'):
            f = gzip.open(filepath, 'rt')
        else:
            f = open(filepath, 'r')
        
        with f:
            # Check first record
            header = f.readline().strip()
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            
            if not header.startswith('@'):
                logging.error(f"Invalid FASTQ format in {filepath}: missing '@' header")
                return False
            if plus != '+':
                logging.error(f"Invalid FASTQ format in {filepath}: missing '+' separator")
                return False
            if len(sequence) != len(quality):
                logging.error(f"Invalid FASTQ format in {filepath}: sequence/quality length mismatch")
                return False
                
        return True
    except Exception as e:
        logging.error(f"Error reading {filepath}: {e}")
        return False


def check_reference_files(ref_dir: str = '.') -> dict:
    """
    Check for all required and optional reference files.
    
    Args:
        ref_dir: Directory containing reference files
        
    Returns:
        Dictionary mapping reference names to file paths
    """
    required_refs = ['vector', 'backbone', 'helper', 'host']
    optional_refs = ['rep_cap', 'spike_in']
    
    refs = {}
    missing = []
    
    # Check required references
    for ref_name in required_refs:
        ref_path = find_reference_file(ref_name, ref_dir)
        if ref_path is None:
            missing.append(ref_name)
        else:
            if validate_fasta(ref_path):
                refs[ref_name] = ref_path
            else:
                missing.append(ref_name)
    
    if missing:
        logging.error("Missing or invalid required reference files: %s", ', '.join(missing))
        logging.error("Required: vector, backbone, helper, host (.fasta or .fa)")
        return None
    
    # Check optional references
    for ref_name in optional_refs:
        ref_path = find_reference_file(ref_name, ref_dir)
        if ref_path is not None and validate_fasta(ref_path):
            refs[ref_name] = ref_path
            logging.info(f"Optional reference detected: {ref_name}")
    
    return refs


def create_output_structure(base_dir: str) -> dict:
    """
    Create organized output directory structure.
    
    Args:
        base_dir: Base output directory
        
    Returns:
        Dictionary mapping directory types to paths
    """
    dirs = {
        'base': base_dir,
        'fastq': os.path.join(base_dir, 'fastq'),
        'bam': os.path.join(base_dir, 'bam'),
        'stats': os.path.join(base_dir, 'stats'),
        'logs': os.path.join(base_dir, 'logs')
    }
    
    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)
    
    return dirs


def run_command(cmd: List[str], description: str = None, check: bool = True) -> subprocess.CompletedProcess:
    """
    Run a shell command with error handling and logging.
    
    Args:
        cmd: Command as list of strings
        description: Human-readable description for logging
        check: Whether to raise exception on non-zero exit
        
    Returns:
        CompletedProcess instance
        
    Raises:
        subprocess.CalledProcessError if command fails and check=True
    """
    if description:
        logging.info(description)
    
    logging.debug("Running: %s", ' '.join(cmd))
    
    try:
        result = subprocess.run(
            cmd,
            check=check,
            capture_output=True,
            text=True
        )
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}")
        logging.error(f"Command: {' '.join(cmd)}")
        if e.stderr:
            logging.error(f"Error output: {e.stderr}")
        raise


def count_fastq_reads(fastq_path: str) -> int:
    """
    Count number of reads in a FASTQ file.
    Handles both .fastq and .fastq.gz files.
    
    Args:
        fastq_path: Path to FASTQ file
        
    Returns:
        Number of reads
    """
    if not os.path.exists(fastq_path) or os.path.getsize(fastq_path) == 0:
        return 0
    
    result = run_command(
        ['seqkit', 'seq', '-n', fastq_path],
        description=None,
        check=True
    )
    return len(result.stdout.strip().split('\n')) if result.stdout.strip() else 0


def get_file_size_mb(filepath: str) -> float:
    """
    Get file size in megabytes.
    
    Args:
        filepath: Path to file
        
    Returns:
        File size in MB
    """
    return os.path.getsize(filepath) / (1024 * 1024)


# ----------------------------------------------------------------------
# Multi-sample input handling
# ----------------------------------------------------------------------

def find_fastq_files(sample_dir: str) -> List[str]:
    """
    Find FASTQ files in a sample directory.
    
    Searches the sample directory itself and one level of subdirectories
    (e.g., sample_dir/ and sample_dir/fastq_pass/). Skips any subdirectory
    named 'references' to avoid mixing reference files with input data.
    Does not recurse deeper than one level — this keeps re-run output
    files (which live multiple levels deep under output/) from being
    picked up as input.
    
    Args:
        sample_dir: Path to the sample directory
        
    Returns:
        Sorted list of absolute paths to FASTQ files found
    """
    extensions = ('.fastq', '.fastq.gz')
    found = []
    
    if not os.path.isdir(sample_dir):
        return found
    
    # Level 0: files directly in sample_dir
    for entry in os.listdir(sample_dir):
        full_path = os.path.join(sample_dir, entry)
        if os.path.isfile(full_path) and entry.lower().endswith(extensions):
            found.append(os.path.abspath(full_path))
    
    # Level 1: files in immediate subdirectories (excluding 'references')
    for entry in os.listdir(sample_dir):
        sub_path = os.path.join(sample_dir, entry)
        if not os.path.isdir(sub_path):
            continue
        if entry == 'references':
            continue
        for sub_entry in os.listdir(sub_path):
            full_path = os.path.join(sub_path, sub_entry)
            if os.path.isfile(full_path) and sub_entry.lower().endswith(extensions):
                found.append(os.path.abspath(full_path))
    
    return sorted(found)


def merge_fastq_files(fastq_files: List[str], output_path: str) -> str:
    """
    Merge multiple FASTQ files into a single .fastq.gz file.
    
    Handles a mix of plain .fastq and gzipped .fastq.gz inputs. Uses
    `cat` for gzipped inputs (gzip streams concatenate losslessly) and
    pipes plain text through gzip for compression. The result is always
    a gzipped output regardless of input compression.
    
    If the output file already exists, this function does NOT overwrite
    it — call sites should check first and decide whether to skip.
    
    Args:
        fastq_files: List of input FASTQ file paths
        output_path: Path where the merged .fastq.gz will be written
        
    Returns:
        Path to the merged output file
        
    Raises:
        ValueError: if fastq_files is empty
        subprocess.CalledProcessError: if the underlying merge command fails
    """
    if not fastq_files:
        raise ValueError("Cannot merge: no FASTQ files provided")
    
    logging.info(f"Merging {len(fastq_files)} FASTQ files into {output_path}")
    
    # Build a shell pipeline that handles mixed compression cleanly:
    # - .gz files: cat them directly (gzip streams concatenate)
    # - plain .fastq files: pipe through gzip
    # We assemble per-file commands and concatenate their outputs to gzip.
    
    gz_files = [f for f in fastq_files if f.lower().endswith('.gz')]
    plain_files = [f for f in fastq_files if not f.lower().endswith('.gz')]
    
    # If everything is gzipped, simple cat is enough
    if plain_files == []:
        cmd = f"cat {' '.join(_shell_quote(f) for f in gz_files)} > {_shell_quote(output_path)}"
    elif gz_files == []:
        # Everything plain: cat then gzip
        cmd = (f"cat {' '.join(_shell_quote(f) for f in plain_files)} "
               f"| gzip > {_shell_quote(output_path)}")
    else:
        # Mixed: gzip the plain ones inline, then cat with the gz ones
        plain_part = (f"cat {' '.join(_shell_quote(f) for f in plain_files)} | gzip"
                      if plain_files else "")
        gz_part = f"cat {' '.join(_shell_quote(f) for f in gz_files)}" if gz_files else ""
        # Order doesn't matter for FASTQ semantically; put gz first then plain
        cmd = f"( {gz_part}; {plain_part} ) > {_shell_quote(output_path)}"
    
    logging.debug(f"Merge command: {cmd}")
    
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"FASTQ merge failed: {e.stderr}")
        # Clean up partial output if it exists
        if os.path.exists(output_path):
            try:
                os.remove(output_path)
            except OSError:
                pass
        raise
    
    size_mb = get_file_size_mb(output_path)
    logging.info(f"Merged FASTQ written: {output_path} ({size_mb:.1f} MB)")
    return output_path


def _shell_quote(s: str) -> str:
    """Quote a path for safe shell interpolation."""
    return "'" + s.replace("'", "'\\''") + "'"


def filter_fastq_by_length(input_fastq: str, output_fastq: str,
                           min_length: int) -> Dict[str, int]:
    """
    Filter a FASTQ file, keeping only reads >= min_length.

    Uses `seqkit seq -m <min_length>`. Counts reads before and after
    so the number dropped can be reported.

    Args:
        input_fastq: Path to input FASTQ (.fastq or .fastq.gz)
        output_fastq: Path to write the filtered FASTQ (.fastq.gz)
        min_length: Minimum read length to keep (reads shorter are dropped)

    Returns:
        Dict with keys: 'input_reads', 'output_reads', 'dropped_reads'

    Raises:
        subprocess.CalledProcessError: if seqkit fails
    """
    logging.info(f"Filtering reads shorter than {min_length} bp")

    input_count = count_fastq_reads(input_fastq)

    # seqkit seq -m filters by minimum length; pipe through gzip for output
    cmd = (f"seqkit seq -m {int(min_length)} {_shell_quote(input_fastq)} "
           f"| gzip > {_shell_quote(output_fastq)}")
    logging.debug(f"Length filter command: {cmd}")

    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Length filter failed: {e.stderr}")
        if os.path.exists(output_fastq):
            try:
                os.remove(output_fastq)
            except OSError:
                pass
        raise

    output_count = count_fastq_reads(output_fastq)
    dropped = input_count - output_count

    logging.info(f"Length filter: kept {output_count:,} / {input_count:,} reads "
                 f"({dropped:,} dropped, < {min_length} bp)")

    return {
        'input_reads': input_count,
        'output_reads': output_count,
        'dropped_reads': dropped
    }


def prepare_sample_fastq(sample_dir: str, sample_name: str,
                         min_length: int = 0,
                         work_dir: str = None) -> Optional[str]:
    """
    Prepare a single FASTQ file for analysis from a sample directory.

    Behavior:
    - If `merged.fastq.gz` already exists in the sample dir, use it
      (idempotent re-runs don't pay the merge cost twice).
    - If exactly one FASTQ file is found, use it directly (no merge).
    - If multiple FASTQ files are found, merge them into
      `<sample_dir>/merged.fastq.gz`. Original input files are not
      modified or deleted.
    - If a minimum length filter is requested (min_length > 0), the
      prepared FASTQ is then filtered into
      `<sample_dir>/merged.filtered.fastq.gz` and that path is returned.
      The unfiltered file is left untouched. A filter stats file
      (`<sample_dir>/length_filter_stats.json`) records how many reads
      were dropped.
    - If no FASTQ files are found, return None and log an error.

    Args:
        sample_dir: Path to the sample directory
        sample_name: Name of the sample (for logging)
        min_length: Minimum read length to keep (0 = no filtering)
        work_dir: Optional directory for writing merged/filtered intermediates
                  and the filter stats file. Defaults to sample_dir (standalone
                  behavior). Workflow managers (Nextflow) pass a writable task
                  directory here so the (possibly read-only / symlink-staged)
                  sample directory is never written to.

    Returns:
        Absolute path to the FASTQ file ready for analysis, or None if
        no FASTQ files were found.
    """
    import json as _json

    # Where intermediate/merged/filtered files get written
    out_dir = work_dir if work_dir else sample_dir
    if work_dir:
        os.makedirs(work_dir, exist_ok=True)

    # ---- Step 1: obtain a single prepared FASTQ (merge if needed) ----
    prepared = None

    # Check for pre-existing merged file first (in the output/work dir)
    pre_merged = os.path.join(out_dir, 'merged.fastq.gz')
    if os.path.exists(pre_merged) and os.path.getsize(pre_merged) > 0:
        logging.info(f"[{sample_name}] Using existing merged FASTQ: {pre_merged}")
        prepared = os.path.abspath(pre_merged)
    else:
        fastq_files = find_fastq_files(sample_dir)

        if not fastq_files:
            logging.error(f"[{sample_name}] No FASTQ files found in {sample_dir}")
            return None

        if len(fastq_files) == 1:
            logging.info(f"[{sample_name}] Single FASTQ file found, using directly: "
                         f"{fastq_files[0]}")
            prepared = fastq_files[0]
        else:
            logging.info(f"[{sample_name}] Found {len(fastq_files)} FASTQ files, merging")
            merged_path = os.path.join(out_dir, 'merged.fastq.gz')
            prepared = merge_fastq_files(fastq_files, merged_path)

    # ---- Step 2: apply length filter if requested ----
    if min_length and min_length > 0:
        filtered_path = os.path.join(out_dir, 'merged.filtered.fastq.gz')

        # Reuse an existing filtered file only if it's newer than the source
        if (os.path.exists(filtered_path) and os.path.getsize(filtered_path) > 0
                and os.path.getmtime(filtered_path) >= os.path.getmtime(prepared)):
            logging.info(f"[{sample_name}] Using existing filtered FASTQ: "
                         f"{filtered_path}")
            return os.path.abspath(filtered_path)

        stats = filter_fastq_by_length(prepared, filtered_path, min_length)

        # Persist filter stats for the report / audit trail
        stats_path = os.path.join(out_dir, 'length_filter_stats.json')
        stats_out = dict(stats)
        stats_out['min_length'] = int(min_length)
        try:
            with open(stats_path, 'w') as f:
                _json.dump(stats_out, f, indent=2)
        except OSError as e:
            logging.warning(f"[{sample_name}] Could not write length filter stats: {e}")

        return os.path.abspath(filtered_path)

    # No filtering requested
    return prepared


def discover_samples(run_dir: str) -> List[Dict[str, str]]:
    """
    Discover samples under a run directory.
    
    A "sample" is any direct subdirectory of `<run_dir>/samples/`. The
    sample's display name is the directory name itself.
    
    Args:
        run_dir: Path to the run directory containing a `samples/` subdir
        
    Returns:
        List of dicts with keys: 'name' (sample name), 'dir' (absolute
        path to the sample directory). Sorted by name.
        
    Raises:
        FileNotFoundError: if `<run_dir>/samples/` does not exist
    """
    samples_root = os.path.join(run_dir, 'samples')
    
    if not os.path.isdir(samples_root):
        raise FileNotFoundError(
            f"Expected samples directory not found: {samples_root}\n"
            f"The run directory must contain a 'samples/' subdirectory "
            f"with one folder per sample."
        )
    
    samples = []
    seen_names_lower = {}  # for case-collision detection
    
    for entry in sorted(os.listdir(samples_root)):
        sample_dir = os.path.join(samples_root, entry)
        if not os.path.isdir(sample_dir):
            continue
        
        # Detect case-only collisions (matters on case-insensitive filesystems
        # but worth flagging even on case-sensitive ones for portability)
        lower = entry.lower()
        if lower in seen_names_lower:
            raise ValueError(
                f"Sample name collision: '{entry}' and "
                f"'{seen_names_lower[lower]}' differ only in case. "
                f"Rename one of them."
            )
        seen_names_lower[lower] = entry
        
        samples.append({
            'name': entry,
            'dir': os.path.abspath(sample_dir)
        })
    
    return samples


def resolve_sample_references(sample_dir: str, run_dir: str,
                              sample_name: str) -> Optional[Dict[str, str]]:
    """
    Resolve reference files for a sample using per-file fallback.
    
    Resolution order for each reference file:
    1. <sample_dir>/references/<name>.{fasta,fa}  (sample-specific override)
    2. <run_dir>/references/<name>.{fasta,fa}      (run-level default)
    3. Not found → required refs cause an error; optional refs are skipped
    
    Logs the resolved source for each reference so users can audit which
    reference was used where.
    
    Args:
        sample_dir: Absolute path to the sample directory
        run_dir: Absolute path to the run directory
        sample_name: Name of the sample (for logging)
        
    Returns:
        Dict mapping reference name to absolute file path, or None if
        any required reference is missing.
    """
    required_refs = ['vector', 'backbone', 'helper', 'host']
    optional_refs = ['rep_cap', 'spike_in']
    
    sample_refs_dir = os.path.join(sample_dir, 'references')
    run_refs_dir = os.path.join(run_dir, 'references')
    
    resolved = {}
    missing_required = []
    
    def _find_in(name, ref_dir):
        if not os.path.isdir(ref_dir):
            return None
        path = find_reference_file(name, ref_dir)
        if path and validate_fasta(path):
            return os.path.abspath(path)
        return None
    
    # Required references
    for ref_name in required_refs:
        sample_path = _find_in(ref_name, sample_refs_dir)
        if sample_path:
            resolved[ref_name] = sample_path
            logging.info(f"[{sample_name}]   {ref_name}: {sample_path} "
                         f"(sample override)")
            continue
        
        run_path = _find_in(ref_name, run_refs_dir)
        if run_path:
            resolved[ref_name] = run_path
            logging.info(f"[{sample_name}]   {ref_name}: {run_path} "
                         f"(run-level)")
            continue
        
        missing_required.append(ref_name)
    
    if missing_required:
        logging.error(f"[{sample_name}] Missing required references: "
                      f"{', '.join(missing_required)}")
        logging.error(f"[{sample_name}] Searched in: {sample_refs_dir} "
                      f"and {run_refs_dir}")
        return None
    
    # Optional references
    for ref_name in optional_refs:
        sample_path = _find_in(ref_name, sample_refs_dir)
        if sample_path:
            resolved[ref_name] = sample_path
            logging.info(f"[{sample_name}]   {ref_name}: {sample_path} "
                         f"(sample override, optional)")
            continue
        
        run_path = _find_in(ref_name, run_refs_dir)
        if run_path:
            resolved[ref_name] = run_path
            logging.info(f"[{sample_name}]   {ref_name}: {run_path} "
                         f"(run-level, optional)")
            continue
        
        logging.info(f"[{sample_name}]   {ref_name}: not provided (skipped)")
    
    return resolved
