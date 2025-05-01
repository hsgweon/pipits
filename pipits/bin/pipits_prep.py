#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import shutil
import logging
import time
import gzip
import bz2
import statistics # Needed for mean
# *** ADDED traceback import for detailed error logging ***
import traceback

__version__    = "4.0" # Removed skipjoin option and placeholders

__author__     = "Hyun Soon Gweon"
__copyright__  = "Copyright 2015-2024, The PIPITS Project"
__credits__    = ["Hyun Soon Gweon", "Anna Oliver", "Joanne Taylor", "Tim Booth", "Melanie Gibbs", "Daniel S. Read", "Robert I. Griffiths", "Karsten Schonrogge"]
__license__    = "GPL"
__maintainer__ = "PIPITS Maintainers"

# Tool Paths (ensure these are correct for the execution environment)
VSEARCH                    = "vsearch"
FASTX_FASTQ_QUALITY_FILTER = "fastq_quality_filter"
FASTX_FASTQ_TO_FASTA       = "fastq_to_fasta"

class Colors:
    RESET = '\033[0m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    HEADER = '\033[95m'

# --- Logging Setup ---
def setup_logging(log_file_path, verbose):
    """Sets up logging to file and console."""
    log_format_file = '%(asctime)s - %(levelname)s - %(message)s'
    log_format_console = '%(levelname)s: %(message)s' # Keep console simple
    log_level = logging.DEBUG if verbose else logging.INFO

    logger = logging.getLogger()
    # Set the base level of the logger to DEBUG to capture everything
    logger.setLevel(logging.DEBUG)

    # Clear previous handlers to avoid duplicate logging
    if logger.hasHandlers():
        logger.handlers.clear()

    # File handler - always logs DEBUG and above
    file_handler = logging.FileHandler(log_file_path, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter(log_format_file)
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    # Console handler - logs at INFO or DEBUG level based on verbosity
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_formatter = logging.Formatter(log_format_console)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

# --- Command Execution ---
def run_cmd(command, verbose):
    """Runs a command using subprocess and logs output."""
    logging.debug(f"Running command: {command}")
    try:
        # Run command, wait for completion, capture output
        process = subprocess.run(command, shell=True, check=True,
                                 capture_output=True, text=True, encoding='utf-8', errors='ignore')

        # Log STDOUT/STDERR only if verbose is enabled
        # In non-verbose, stderr from tools might still indicate warnings/info
        stdout_lines = process.stdout.strip().split('\n') if process.stdout else []
        stderr_lines = process.stderr.strip().split('\n') if process.stderr else []

        if verbose:
            # Log all non-empty lines from stdout/stderr in verbose mode
            if stdout_lines:
                logging.debug("Command STDOUT:")
                for line in stdout_lines:
                    if line: logging.debug(line)
            if stderr_lines:
                logging.debug("Command STDERR:")
                for line in stderr_lines:
                     if line: logging.debug(line)
        else:
             # In non-verbose mode, still log non-empty STDERR lines at DEBUG level
             # as they might contain useful info/warnings from tools like VSEARCH
             if stderr_lines and any(line for line in stderr_lines):
                 logging.debug("Command STDERR (non-verbose):")
                 for line in stderr_lines:
                     if line: logging.debug(line)

        return process

    except subprocess.CalledProcessError as e:
        # Log detailed error information if command fails
        logging.error(f"Command failed with exit code {e.returncode}")
        logging.error(f"Failed command: {e.cmd}")
        # Capture and log stderr/stdout from the failed command
        stderr_output = e.stderr.strip() if e.stderr else "N/A"
        stdout_output = e.stdout.strip() if e.stdout else "N/A"
        logging.error(f"STDERR:\n---\n{stderr_output}\n---")
        if stdout_output != "N/A": # Only log stdout if it contains something
            logging.error(f"STDOUT:\n---\n{stdout_output}\n---")
        raise # Re-raise the exception to be caught by the main handler
    except FileNotFoundError:
         # Handle case where the command itself is not found
         cmd_name = command.split()[0]
         logging.error(f"Command not found: '{cmd_name}'. Please ensure '{cmd_name}' (e.g., {VSEARCH}, {FASTX_FASTQ_QUALITY_FILTER}) is installed and in your system's PATH.")
         raise # Re-raise
    except Exception as e:
        # Catch any other unexpected errors during command execution
        logging.error(f"An unexpected error occurred running command: {command}")
        logging.error(f"Error details: {e}")
        raise # Re-raise

# --- Sequence Utilities ---
def blocks(files, size=65536):
    """Helper function to read files in chunks (binary mode)."""
    while True:
        b = files.read(size)
        if not b: break
        yield b

def get_file_line_count(filename):
    """Counts lines in a potentially compressed file."""
    opener = open
    mode = "rb" # Read as bytes for efficient newline counting
    try:
        # Choose the correct opener based on file extension
        if filename.lower().endswith(".gz"):
            opener = gzip.open
        elif filename.lower().endswith(".bz2"):
            opener = bz2.open
        # Add other potential compressed formats if needed
        # else: assume plain text

        with opener(filename, mode) as f:
            # Count the occurrences of the newline byte
            count = sum(bl.count(b"\n") for bl in blocks(f))
        return count

    except FileNotFoundError:
        # Log as warning if file not found during counting, it might be expected sometimes
        logging.warning(f"File not found for line count: {filename}")
        return 0
    except Exception as e:
        # Log error for other issues during counting
        logging.error(f"Error counting lines in {filename}: {e}")
        return 0 # Return 0 on error

# --- Core Functions ---

def count_sequences(input_dir, sampleids_list, fastqs_f_list, fastqs_r_list, summary_file_handle):
    """
    Counts sequences in forward and reverse FASTQ files listed for each sample.
    Calculates and logs total, min, max, and average sequences per sample.
    Assumes input files are FASTQ (4 lines per sequence).
    """
    total_sequences = 0
    sequences_per_sample = [] # Store count for each sample (Fwd + Rev combined)

    logging.info("Counting initial sequences from provided list file...")

    # Basic check for list consistency
    if len(sampleids_list) != len(fastqs_f_list) or len(sampleids_list) != len(fastqs_r_list):
         logging.critical("Mismatch in lengths of sample ID and FASTQ file lists. Cannot proceed with counting.")
         summary_file_handle.write("Initial sequence count: ERROR - List length mismatch\n")
         return 0 # Indicate error / zero count

    for i in range(len(sampleids_list)):
        sample_id = sampleids_list[i]
        fwd_filename = fastqs_f_list[i]
        rev_filename = fastqs_r_list[i]

        # Construct full paths
        fwd_full_path = os.path.join(input_dir, fwd_filename)
        rev_full_path = os.path.join(input_dir, rev_filename)

        # Get line counts (handles file not found internally)
        fwd_lines = get_file_line_count(fwd_full_path)
        rev_lines = get_file_line_count(rev_full_path)

        # Calculate sequences assuming FASTQ format (4 lines/seq)
        fwd_seqs = 0
        if fwd_lines > 0:
            if fwd_lines % 4 == 0:
                fwd_seqs = fwd_lines // 4
            else:
                # Log warning if line count is invalid for FASTQ
                logging.warning(f"Forward file {fwd_filename} (sample {sample_id}) has line count ({fwd_lines}) not divisible by 4. Treating as 0 sequences.")

        rev_seqs = 0
        if rev_lines > 0:
            if rev_lines % 4 == 0:
                rev_seqs = rev_lines // 4
            else:
                logging.warning(f"Reverse file {rev_filename} (sample {sample_id}) has line count ({rev_lines}) not divisible by 4. Treating as 0 sequences.")

        # Calculate total for this sample and add to list for stats
        sample_total_seqs = fwd_seqs + rev_seqs
        logging.debug(f"Sample {sample_id}: Fwd={fwd_seqs}, Rev={rev_seqs}, Total={sample_total_seqs}")
        sequences_per_sample.append(sample_total_seqs)
        total_sequences += sample_total_seqs

    # Calculate Statistics only if we have counts
    if not sequences_per_sample:
        logging.error("Could not count sequences for any sample. Check input files listed in the list file exist and are valid FASTQ.")
        summary_file_handle.write("Initial sequence count: ERROR - No valid sequences counted\n")
        return 0

    min_seqs = min(sequences_per_sample)
    max_seqs = max(sequences_per_sample)
    avg_seqs = statistics.mean(sequences_per_sample) # Calculate mean

    # Log results clearly
    logging.info(f"{Colors.BLUE}... Total number of input reads (Fwd + Rev): {total_sequences}{Colors.RESET}")
    logging.info(f"{Colors.BLUE}... Minimum reads per sample: {min_seqs}{Colors.RESET}")
    logging.info(f"{Colors.BLUE}... Maximum reads per sample: {max_seqs}{Colors.RESET}")
    logging.info(f"{Colors.BLUE}... Average reads per sample: {avg_seqs:.2f}{Colors.RESET}")

    # Write results to summary file
    summary_file_handle.write(f"Total number of input reads (Fwd + Rev): {total_sequences}\n")
    summary_file_handle.write(f"Minimum reads per sample: {min_seqs}\n")
    summary_file_handle.write(f"Maximum reads per sample: {max_seqs}\n")
    summary_file_handle.write(f"Average reads per sample: {avg_seqs:.2f}\n")

    # Add a warning if the total count is zero after processing all samples
    if total_sequences == 0:
        logging.warning("Total initial sequence count is 0 after checking all listed files. Subsequent steps might produce empty output.")

    return total_sequences


def reindex_fastq(input_dir, output_dir, sampleids_list, filenames_list):
    """Reindexes FASTQ files, naming them based on sample IDs."""
    logging.info(f"Reindexing files into directory: {output_dir}")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Default opener and mode for plain text
    default_opener = open
    default_mode = 'rt' # Read text mode

    for i in range(len(filenames_list)):
        in_filename = filenames_list[i]
        in_filepath = os.path.join(input_dir, in_filename)
        sample_id = sampleids_list[i]
        # Output filename format: SampleID.fastq
        out_filename = f"{sample_id}.fastq"
        out_filepath = os.path.join(output_dir, out_filename)

        logging.debug(f"Reindexing {in_filename} to {out_filename}")

        # Determine opener and mode based on input file extension
        opener = default_opener
        file_open_mode = default_mode
        if in_filename.lower().endswith(".gz"):
            opener = gzip.open
            file_open_mode = 'rt' # gzip opens in text mode with 'rt'
        elif in_filename.lower().endswith(".bz2"):
            opener = bz2.open
            file_open_mode = 'rt' # bz2 opens in text mode with 'rt'
        elif not in_filename.lower().endswith((".fastq", ".fq")): # Log if not typical FASTQ extension
             logging.warning(f"Processing file with non-standard extension during reindexing: {in_filename}. Assuming plain text FASTQ.")

        try:
            # Use utf-8 encoding, ignore errors for robustness
            with opener(in_filepath, file_open_mode, encoding='utf-8', errors='ignore') as infile, \
                 open(out_filepath, "w", encoding='utf-8') as outfile: # Output is always plain text FASTQ

                line_number = 0
                sequence_number = 1
                for line in infile:
                    line_number += 1
                    # Process based on line number within FASTQ record (1-based index % 4)
                    line_index_in_record = line_number % 4
                    if line_index_in_record == 1: # Header line
                        # Basic check for '@' start
                        if not line.startswith('@'):
                             logging.warning(f"Line {line_number} in {in_filename} does not start with '@'. Input may not be FASTQ.")
                        # Write new header: @SampleID_SequenceNumber
                        outfile.write(f"@{sample_id}_{sequence_number}\n")
                        sequence_number += 1
                    elif line_index_in_record == 2: # Sequence line
                         outfile.write(line.rstrip('\r\n') + "\n") # Ensure single newline
                    elif line_index_in_record == 3: # Separator line
                        # Standardize separator to just '+'
                        outfile.write("+\n")
                    elif line_index_in_record == 0: # Quality score line (line 4, 8, etc.)
                         outfile.write(line.rstrip('\r\n') + "\n") # Ensure single newline

            # Post-processing checks
            if line_number == 0:
                 logging.warning(f"Input file {in_filepath} was empty. Output {out_filepath} is empty.")
            elif line_number % 4 != 0:
                 # This indicates a truncated or invalid FASTQ file
                 logging.warning(f"Input file {in_filepath} has {line_number} lines, not divisible by 4. Output {out_filepath} may be incomplete or corrupt.")

            # Double-check if output file is empty even if input wasn't (could happen on error)
            if os.path.exists(out_filepath) and os.stat(out_filepath).st_size == 0 and line_number > 0:
                 logging.warning(f"Output file {out_filepath} is empty after reindexing {in_filename}, although input seemed to have lines.")

        except FileNotFoundError:
             logging.error(f"Input file not found during reindexing: {in_filepath}. Creating empty output.")
             # Create empty file so downstream steps don't fail immediately on file not found
             open(out_filepath, 'a').close()
        except Exception as e:
             # Catch other errors like permission issues, corrupted compressed files, etc.
             logging.error(f"Error reindexing {in_filename}: {e}. Creating empty output.")
             # Attempt to remove potentially corrupt partial file before creating empty one
             try:
                 if os.path.exists(out_filepath): os.remove(out_filepath)
             except OSError: pass # Ignore error if removal fails
             open(out_filepath, 'a').close() # Create empty file

def join_reads(input_dir_f, input_dir_r, output_dir, sampleids_list,
               threads, verbose, summary_file_handle):
    """Joins paired-end reads using VSEARCH."""
    logging.info("Joining paired-end reads [VSEARCH]")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    total_joined_sequences = 0
    successful_joins = 0
    failed_joins = 0

    for sample_id in sampleids_list:
        fwd_read_path = os.path.join(input_dir_f, f"{sample_id}.fastq")
        rev_read_path = os.path.join(input_dir_r, f"{sample_id}.fastq")
        # Define the path for successfully joined reads
        final_joined_path = os.path.join(output_dir, f"{sample_id}.fastq")
        # Define paths for VSEARCH report files (useful for debugging)
        report_path = os.path.join(output_dir, f"{sample_id}.join.log")

        # Check if input files exist and are not empty before attempting join
        fwd_exists = os.path.exists(fwd_read_path) and os.stat(fwd_read_path).st_size > 0
        rev_exists = os.path.exists(rev_read_path) and os.stat(rev_read_path).st_size > 0

        if not fwd_exists or not rev_exists:
            if not fwd_exists: logging.warning(f"Forward read missing or empty: {fwd_read_path}")
            if not rev_exists: logging.warning(f"Reverse read missing or empty: {rev_read_path}")
            logging.warning(f"Skipping join for sample {sample_id}. Creating empty joined output file.")
            # Create empty output file for consistency downstream
            open(final_joined_path, 'a').close()
            failed_joins += 1
            continue # Move to the next sample

        # --- Run VSEARCH fastq_mergepairs ---
        try:
            # Construct command parts for clarity
            # Using typical PIPITS/PISPINO defaults for joining parameters
            cmd_parts = [
                VSEARCH,
                "--fastq_mergepairs", f'"{fwd_read_path}"', # Input forward reads
                "--reverse", f'"{rev_read_path}"',         # Input reverse reads
                "--fastqout", f'"{final_joined_path}"',    # Output for successfully merged reads
                # Optional: Output for reads that didn't merge (uncomment if needed)
                # "--fastqout_notmerged_fwd", f'"{os.path.join(output_dir, f"{sample_id}_notmerged_fwd.fastq")}"',
                # "--fastqout_notmerged_rev", f'"{os.path.join(output_dir, f"{sample_id}_notmerged_rev.fastq")}"',
                "--log", f'"{report_path}"',               # Log file for merge statistics
                "--threads", str(threads),                 # Number of threads
                # --- Joining Parameters ---
                "--fastq_allowmergestagger", # Allow merging reads that are staggered relative to each other
                "--fastq_maxdiffs", "500",   # Maximum number of mismatches allowed in the overlap region (high value)
                "--fastq_minovlen", "20",    # Minimum length of the overlap region required
                "--fastq_minmergelen", "100" # Minimum length of the final merged sequence
                # Add other VSEARCH options if needed
            ]
            cmd = " ".join(cmd_parts)
            # Run the command (will raise exception on failure)
            run_cmd(cmd, verbose)

            # Count sequences in the final joined file *if* it was created
            num_seqs_joined = 0
            if os.path.exists(final_joined_path):
                 line_count = get_file_line_count(final_joined_path)
                 if line_count > 0 and line_count % 4 == 0:
                     num_seqs_joined = line_count // 4
                     logging.debug(f"Sample {sample_id}: Found {num_seqs_joined} joined sequences.")
                 elif line_count > 0: # Invalid FASTQ line count
                      logging.warning(f"Sample {sample_id}: Joined file {final_joined_path} has unexpected line count ({line_count}).")
                 # If line_count is 0, num_seqs_joined remains 0
            else:
                 # Vsearch might not create the file if zero reads merge successfully
                 logging.debug(f"Joined output file not created for {sample_id} (zero reads merged?): {final_joined_path}")
                 # Ensure empty file exists if VSEARCH didn't create it
                 open(final_joined_path, 'a').close()

            total_joined_sequences += num_seqs_joined
            successful_joins += 1

        except Exception as e:
            # Catch errors from run_cmd (CalledProcessError, FileNotFoundError, etc.)
            logging.error(f"Error joining reads for sample {sample_id}: {e}")
            failed_joins += 1
            # Ensure an empty joined file exists even if the process failed mid-way
            if not os.path.exists(final_joined_path):
                open(final_joined_path, 'a').close()
            # Optionally try to clean up report file?
            # if os.path.exists(report_path): try: os.remove(report_path) except OSError: pass


    # Log summary results for joining step
    logging.info(f"Joining complete. Samples attempted: {len(sampleids_list)}, Successful: {successful_joins}, Failed/Skipped: {failed_joins}")
    if total_joined_sequences == 0 and successful_joins > 0:
        # This case indicates joins ran but yielded no merged sequences across all samples
        logging.warning("VSEARCH JOINING STEP YIELDED 0 SEQUENCES ACROSS ALL SAMPLES. Check input data quality, read orientation, and joining parameters (min overlap, max diffs, etc.). See VSEARCH logs in output directory for details.")
    elif successful_joins == 0:
        logging.error("VSEARCH JOINING FAILED OR WAS SKIPPED FOR ALL SAMPLES.")
    else:
        # Log the total count if > 0
        logging.info(f"{Colors.BLUE}... Total number of joined reads (VSEARCH): {total_joined_sequences}{Colors.RESET}")
        summary_file_handle.write(f"Number of joined reads (VSEARCH): {total_joined_sequences}\n")

    # Return the total count for potential use later (or just rely on logging)
    return total_joined_sequences


def quality_filter(input_dir, output_dir, sampleids_list,
                   base_phred_quality_score, FASTX_fastq_quality_filter_q,
                   FASTX_fastq_quality_filter_p, verbose, summary_file_handle):
    """Filters sequences based on quality using FASTX-Toolkit."""
    logging.info("Quality filtering [FASTX-Toolkit fastq_quality_filter]")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    total_filtered_sequences = 0
    successful_filters = 0
    failed_filters = 0

    for sample_id in sampleids_list:
        # Input is the result of the previous step (joining)
        input_file = os.path.join(input_dir, f"{sample_id}.fastq")
        # Output for filtered reads
        output_file = os.path.join(output_dir, f"{sample_id}.fastq")

        # Check if input file exists and is not empty
        if not os.path.exists(input_file) or os.stat(input_file).st_size == 0:
            logging.warning(f"Input file for quality filtering not found or empty: {input_file}. Skipping filter for {sample_id}, creating empty output.")
            open(output_file, 'a').close() # Create empty output
            failed_filters += 1
            continue # Skip to next sample

        try:
            # Construct command parts
            cmd_parts = [
                FASTX_FASTQ_QUALITY_FILTER,
                "-i", f'"{input_file}"',             # Input FASTQ
                "-o", f'"{output_file}"',            # Output FASTQ
                "-q", str(FASTX_fastq_quality_filter_q), # Minimum quality score
                "-p", str(FASTX_fastq_quality_filter_p), # Minimum percent of bases at that score
                "-Q", str(base_phred_quality_score)  # Input quality score encoding (e.g., 33 or 64)
                # Add -v for verbose output from the tool itself if needed for debugging
                # "-v"
            ]
            cmd = " ".join(cmd_parts)
            # Run the command (raises exception on failure)
            run_cmd(cmd, verbose)

            # Count sequences in the filtered output file
            num_seqs_filtered = 0
            if os.path.exists(output_file):
                line_count = get_file_line_count(output_file)
                if line_count > 0 and line_count % 4 == 0:
                    num_seqs_filtered = line_count // 4
                    logging.debug(f"Sample {sample_id}: Found {num_seqs_filtered} quality-filtered sequences.")
                elif line_count > 0: # Invalid FASTQ line count
                     logging.warning(f"Sample {sample_id}: Quality-filtered file {output_file} has unexpected line count ({line_count}).")
                # If line_count is 0, num_seqs_filtered remains 0
            else:
                # fastq_quality_filter might not create output if zero reads pass
                logging.debug(f"Quality-filtered output file not created for {sample_id} (zero reads passed?): {output_file}")
                # Ensure empty file exists
                open(output_file, 'a').close()

            total_filtered_sequences += num_seqs_filtered
            successful_filters += 1

        except Exception as e:
             # Catch errors from run_cmd or other issues
             logging.error(f"Error during quality filtering for sample {sample_id}: {e}")
             failed_filters += 1
             # Ensure empty output file exists on error
             if not os.path.exists(output_file): open(output_file, 'a').close()

    # Log summary for filtering step
    logging.info(f"Quality filtering complete. Samples attempted: {len(sampleids_list)}, Successful: {successful_filters}, Failed/Skipped: {failed_filters}")
    if total_filtered_sequences == 0 and successful_filters > 0:
        # Indicates filtering ran but removed all sequences
        logging.warning("QUALITY FILTERING YIELDED 0 SEQUENCES ACROSS ALL SAMPLES. Check quality parameters (q, p, Q) and input data quality.")
    elif successful_filters == 0:
        logging.error("QUALITY FILTERING FAILED OR WAS SKIPPED FOR ALL SAMPLES.")
    else:
        # Log total count if > 0
        logging.info(f"{Colors.BLUE}... Total number of quality-filtered reads: {total_filtered_sequences}{Colors.RESET}")
        summary_file_handle.write(f"Number of quality filtered reads: {total_filtered_sequences}\n")

    return total_filtered_sequences


def convert_to_fasta(input_dir, output_dir, sampleids_list,
                     base_phred_quality_score, FASTX_fastq_to_fasta_keep_n,
                     verbose, summary_file_handle):
    """Converts FASTQ to FASTA using FASTX-Toolkit, optionally keeping Ns."""
    logging.info("Converting FASTQ to FASTA [FASTX-Toolkit fastq_to_fasta]")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    total_converted_sequences = 0
    successful_conversions = 0
    failed_conversions = 0

    for sample_id in sampleids_list:
        # Input is the result of quality filtering
        input_file = os.path.join(input_dir, f"{sample_id}.fastq")
        # Output is FASTA format
        output_file = os.path.join(output_dir, f"{sample_id}.fasta")

        # Check if input FASTQ file exists and is not empty
        if not os.path.exists(input_file) or os.stat(input_file).st_size == 0:
            logging.warning(f"Input file for FASTA conversion not found or empty: {input_file}. Skipping conversion for {sample_id}, creating empty output.")
            open(output_file, 'a').close() # Create empty output
            failed_conversions += 1
            continue # Skip to next sample

        try:
            # Construct command parts
            cmd_parts = [
                FASTX_FASTQ_TO_FASTA,
                "-i", f'"{input_file}"',             # Input FASTQ
                "-o", f'"{output_file}"',            # Output FASTA
                "-Q", str(base_phred_quality_score)  # Input quality score encoding
            ]
            # Logic based on the --keep-Ns flag:
            # -n flag in fastq_to_fasta *removes* sequences containing Ns.
            # So, if user specified --keep-Ns (keep_n is True), we *don't* add the -n flag.
            # If user did *not* specify --keep-Ns (keep_n is False), we *do* add the -n flag.
            if not FASTX_fastq_to_fasta_keep_n:
                 cmd_parts.append("-n") # Add flag to discard sequences with Ns
                 logging.debug(f"Sample {sample_id}: Removing sequences with Ns during FASTA conversion.")
            else:
                 logging.debug(f"Sample {sample_id}: Keeping sequences with Ns during FASTA conversion.")

            cmd = " ".join(cmd_parts)
            # Run the command (raises exception on failure)
            run_cmd(cmd, verbose)

            # Count sequences in the converted FASTA file (2 lines per sequence)
            num_seqs_converted = 0
            if os.path.exists(output_file):
                line_count = get_file_line_count(output_file)
                if line_count > 0 and line_count % 2 == 0:
                    num_seqs_converted = line_count // 2
                    logging.debug(f"Sample {sample_id}: Found {num_seqs_converted} sequences after FASTA conversion.")
                elif line_count > 0: # Invalid FASTA line count (odd number)
                     logging.warning(f"Sample {sample_id}: Converted FASTA file {output_file} has unexpected line count ({line_count}).")
                # If line_count is 0, num_seqs_converted remains 0
            else:
                # fastq_to_fasta might not create output if zero reads pass/convert
                logging.debug(f"Converted FASTA file not created for {sample_id} (zero reads converted?): {output_file}")
                # Ensure empty file exists
                open(output_file, 'a').close()

            total_converted_sequences += num_seqs_converted
            successful_conversions += 1

        except Exception as e:
             # Catch errors from run_cmd or other issues
             logging.error(f"Error during FASTA conversion for sample {sample_id}: {e}")
             failed_conversions += 1
             # Ensure empty output file exists on error
             if not os.path.exists(output_file): open(output_file, 'a').close()

    # Log summary for conversion step
    logging.info(f"FASTA conversion complete. Samples attempted: {len(sampleids_list)}, Successful: {successful_conversions}, Failed/Skipped: {failed_conversions}")
    if total_converted_sequences == 0 and successful_conversions > 0:
        # Indicates conversion ran but yielded no sequences
        logging.warning("FASTA CONVERSION YIELDED 0 SEQUENCES ACROSS ALL SAMPLES. This might happen if quality filtering removed all reads, or if the '-n' flag removed all sequences with Ns.")
    elif successful_conversions == 0:
        logging.error("FASTA CONVERSION FAILED OR WAS SKIPPED FOR ALL SAMPLES.")
    else:
        # Log total count if > 0
        logging.info(f"{Colors.BLUE}... Total number of final prepped sequences: {total_converted_sequences}{Colors.RESET}")
        summary_file_handle.write(f"Number of final prepped sequences: {total_converted_sequences}\n")

    return total_converted_sequences


def merge_fasta(input_dir, output_dir, sampleids_list, summary_file_handle, final_output_name="prepped.fasta"):
    """Merges individual FASTA files into a single file."""
    final_output_path = os.path.join(output_dir, final_output_name)
    logging.info(f"Merging FASTA files from {input_dir} into {final_output_path}")

    total_merged_sequences = 0
    files_merged = 0
    files_missing_or_empty = 0

    try:
        # Open the final output file for writing
        with open(final_output_path, "w", encoding='utf-8') as outfile:
            # Iterate through each sample ID to find its corresponding FASTA file
            for sample_id in sampleids_list:
                input_file = os.path.join(input_dir, f"{sample_id}.fasta")

                # Check if the individual FASTA file exists and has content
                if os.path.exists(input_file) and os.stat(input_file).st_size > 0:
                    logging.debug(f"Merging {input_file}")
                    files_merged += 1
                    seq_count_in_file = 0
                    line_count_in_file = 0
                    # Open the individual file and append its content
                    try:
                        with open(input_file, "r", encoding='utf-8', errors='ignore') as infile:
                            for line in infile:
                                outfile.write(line)
                                line_count_in_file += 1
                        # Calculate sequence count for this file
                        if line_count_in_file > 0 and line_count_in_file % 2 == 0:
                             seq_count_in_file = line_count_in_file // 2
                        elif line_count_in_file > 0:
                             logging.warning(f"Input FASTA {input_file} had odd line count ({line_count_in_file}). Sequence count for this file might be inaccurate.")
                             seq_count_in_file = line_count_in_file // 2 # Estimate anyway

                        total_merged_sequences += seq_count_in_file
                        logging.debug(f"Merged {seq_count_in_file} sequences from {input_file}")

                    except IOError as e:
                         logging.error(f"Error reading input FASTA file during merge: {input_file} - {e}")
                         # Continue to next file even if one fails to read
                else:
                    # Log if a file is missing or empty
                    logging.warning(f"Input file for merging not found or empty, skipping: {input_file}")
                    files_missing_or_empty += 1

        # Log summary of merging process
        logging.info(f"Merging complete. Merged {files_merged} files. Skipped {files_missing_or_empty} missing/empty files.")
        if files_merged == 0:
             logging.error(f"No valid FASTA files found to merge in {input_dir}. Final output file {final_output_path} will be empty.")
             # Ensure file exists even if empty
             if not os.path.exists(final_output_path): open(final_output_path, 'a').close()
        else:
             logging.info(f"{Colors.BLUE}... Total sequences in merged file: {total_merged_sequences}{Colors.RESET}")
             logging.info(f"{Colors.BLUE}... Final output file: {final_output_path}{Colors.RESET}")
             summary_file_handle.write(f"Number of sequences in final merged file: {total_merged_sequences}\n")

        # Final check on the merged file itself
        if os.path.exists(final_output_path):
             final_file_line_count = get_file_line_count(final_output_path)
             if final_file_line_count % 2 != 0 and final_file_line_count > 0:
                 logging.warning(f"Final merged file {final_output_path} has an odd number of lines ({final_file_line_count}), indicating potential corruption.")
             elif final_file_line_count == 0 and total_merged_sequences > 0:
                  logging.error(f"Calculated {total_merged_sequences} merged sequences, but the final file {final_output_path} is empty.")
             elif final_file_line_count // 2 != total_merged_sequences:
                  logging.warning(f"Mismatch between calculated merged sequence count ({total_merged_sequences}) and count based on final file lines ({final_file_line_count // 2}).")
        else:
             logging.error(f"Final merged file {final_output_path} was not created.")


    except IOError as e:
        logging.error(f"Error writing final merged FASTA file {final_output_path}: {e}")
        # Ensure final output path exists even on error, possibly empty
        if not os.path.exists(final_output_path): open(final_output_path, 'a').close()
    except Exception as e:
        logging.error(f"An unexpected error occurred during FASTA merging: {e}")
        if not os.path.exists(final_output_path): open(final_output_path, 'a').close()


# *** WRAPPER FUNCTION for the main pipeline logic ***
def run_prep_pipeline(options):
    """Contains the main execution logic called by main()."""

    # --- Directory Setup ---
    if not os.path.exists(options.outputdir):
        os.makedirs(options.outputdir)
    tmpDir = os.path.join(options.outputdir, "tmp")
    if not options.retain and os.path.exists(tmpDir):
         logging.info(f"Removing existing temporary directory: {tmpDir}")
         shutil.rmtree(tmpDir)
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir)

    # --- Logging Setup (already done in main, but paths are defined here) ---
    summary_file_path = os.path.join(options.outputdir, "summary.log")

    # --- Start Timer & Log Basic Info ---
    start_time = time.time()
    logging.info(f"--- Starting pipits_prep v{__version__} ---")
    logging.info(f"Output directory: {os.path.abspath(options.outputdir)}")
    logging.info(f"Input data directory: {os.path.abspath(options.dataDir)}")
    logging.info(f"List file: {os.path.abspath(options.listfile)}")
    logging.info("Read joining method: VSEARCH")
    logging.info(f"Quality filter (FASTX): q={options.FASTX_fastq_quality_filter_q}, p={options.FASTX_fastq_quality_filter_p}, Q={options.base_phred_quality_score}")
    logging.info(f"Keep Ns during FASTA conversion: {options.FASTX_fastq_to_fasta_keep_n}")
    logging.info(f"Threads for joining: {options.threads}")
    logging.info(f"Retain intermediates: {options.retain}")
    logging.info(f"Verbose logging: {options.verbose}")
    logging.info("-" * 30)


    # --- Check Input Directory ---
    if not os.path.isdir(options.dataDir):
        logging.critical(f"Input data directory not found or not a directory: {options.dataDir}")
        return 1 # Exit error

    # --- Read List File ---
    sampleids = []
    fastqs_f = []
    fastqs_r = []
    logging.info(f"Reading sample list file: {options.listfile}")
    try:
        # Use utf-8 encoding, ignore errors for robustness
        with open(options.listfile, "r", encoding='utf-8', errors='ignore') as lf:
            line_num = 0
            for line in lf:
                line_num += 1
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith("#"): continue

                parts = line.split("\t")
                if len(parts) != 3:
                    logging.critical(f"Incorrect format in list file '{options.listfile}' on line {line_num}. Expected 3 tab-separated columns (SampleID, FwdFile, RevFile). Found {len(parts)} columns.")
                    return 1 # Exit error

                sample_id, fwd_file, rev_file = parts[0], parts[1], parts[2]

                # Validate Sample ID format - Crucial for PIPITS
                if not sample_id:
                     logging.critical(f"Empty SampleID found on line {line_num} in list file '{options.listfile}'.")
                     return 1
                if "_" in sample_id:
                     # Underscores in SampleID often break downstream parsing (e.g., uc2otutable)
                     logging.critical(f"SampleID '{sample_id}' on line {line_num} contains forbidden character '_'. Please rename.")
                     return 1
                if " " in sample_id:
                     logging.critical(f"SampleID '{sample_id}' on line {line_num} contains a forbidden space character.")
                     return 1

                # Check file existence *now* to provide early feedback
                fwd_path_check = os.path.join(options.dataDir, fwd_file)
                rev_path_check = os.path.join(options.dataDir, rev_file)

                if not os.path.exists(fwd_path_check):
                     # Log as warning, but continue reading list file. Reindexing step will handle missing file.
                     logging.warning(f"Forward read file '{fwd_file}' (sample '{sample_id}', line {line_num}) not found in '{options.dataDir}'.")
                if not os.path.exists(rev_path_check):
                     logging.warning(f"Reverse read file '{rev_file}' (sample '{sample_id}', line {line_num}) not found in '{options.dataDir}'.")

                # Add sample ID and filenames to lists
                sampleids.append(sample_id)
                fastqs_f.append(fwd_file)
                fastqs_r.append(rev_file)

        if not sampleids:
             # If the list file was empty or only contained comments/invalid lines
             logging.critical(f"No valid sample entries found in list file: {options.listfile}")
             return 1 # Exit error

        logging.info(f"Found {len(sampleids)} samples in list file.")

        # Check for duplicate SampleIDs - crucial for correct processing
        if len(sampleids) != len(set(sampleids)):
             duplicates = sorted(list(set([s for s in sampleids if sampleids.count(s) > 1])))
             logging.critical(f"Duplicate SampleIDs found in list file '{options.listfile}'. Duplicates: {', '.join(duplicates)}. Please ensure all SampleIDs are unique.")
             return 1 # Exit error
        logging.debug("No duplicate SampleIDs found in list file.")

    except FileNotFoundError:
        logging.critical(f"List file not found: {options.listfile}")
        return 1 # Exit error
    except Exception as e:
        # Catch other errors during list file reading/parsing
        logging.critical(f"Error reading list file {options.listfile}: {e}")
        if options.verbose: logging.debug(traceback.format_exc())
        return 1 # Exit error

    # --- Open Summary File ---
    summary_file_handle = None # Initialize
    try:
        summary_file_handle = open(summary_file_path, "w")
        summary_file_handle.write(f"# pipits_prep v{__version__} Summary\n")
        summary_file_handle.write(f"# Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        summary_file_handle.write(f"# Input directory: {os.path.abspath(options.dataDir)}\n")
        summary_file_handle.write(f"# Output directory: {os.path.abspath(options.outputdir)}\n")
        summary_file_handle.write(f"# List file: {os.path.abspath(options.listfile)}\n")
        summary_file_handle.write(f"# Samples processed: {len(sampleids)}\n")
        summary_file_handle.write("# --- Parameters ---\n")
        summary_file_handle.write(f"Quality filter min score (q): {options.FASTX_fastq_quality_filter_q}\n")
        summary_file_handle.write(f"Quality filter min percent (p): {options.FASTX_fastq_quality_filter_p}\n")
        summary_file_handle.write(f"Keep Ns during conversion: {options.FASTX_fastq_to_fasta_keep_n}\n")
        summary_file_handle.write(f"PHRED base: {options.base_phred_quality_score}\n")
        summary_file_handle.write(f"Threads (for joining): {options.threads}\n")
        summary_file_handle.write("# --- Processing Summary ---\n")

    except IOError as e:
         # If summary file cannot be opened, log critical error and exit
         logging.critical(f"Cannot open summary file for writing: {summary_file_path} - {e}")
         # No need to close handle here as it wasn't successfully opened
         return 1 # Exit error

    # === Start Pipeline Execution ===
    pipeline_failed = False
    final_seq_count = 0 # Initialize final count

    try: # Wrap the entire pipeline in a try block
        initial_seq_count = count_sequences(
            options.dataDir, sampleids, fastqs_f, fastqs_r, summary_file_handle
        )
        if initial_seq_count == 0:
             # Log warning but continue, subsequent steps should handle empty inputs
             logging.warning("Initial sequence count reported as 0. Processing will continue, but check input files and list.")

        # --- Step 1: Reindexing ---
        reindex_f_dir = os.path.join(tmpDir, "1_reindexed_F")
        reindex_r_dir = os.path.join(tmpDir, "2_reindexed_R")
        # Reindex forward reads
        reindex_fastq(options.dataDir, reindex_f_dir, sampleids, fastqs_f)
        # Reindex reverse reads
        reindex_fastq(options.dataDir, reindex_r_dir, sampleids, fastqs_r)
        # Note: reindex_fastq handles missing input files by creating empty outputs

        # --- Step 2: Joining ---
        joined_dir = os.path.join(tmpDir, "3_joined")
        seq_count_after_join = join_reads(
            reindex_f_dir, reindex_r_dir, joined_dir, sampleids,
            options.threads, options.verbose, summary_file_handle
        )
        # join_reads handles empty inputs and logs warnings/errors

        # --- Step 3: Quality Filtering ---
        qfiltered_dir = os.path.join(tmpDir, "4_quality_filtered")
        seq_count_after_qfilter = quality_filter(
            joined_dir, qfiltered_dir, sampleids,
            options.base_phred_quality_score, options.FASTX_fastq_quality_filter_q,
            options.FASTX_fastq_quality_filter_p, options.verbose, summary_file_handle
        )
        # quality_filter handles empty inputs and logs warnings/errors

        # --- Step 4: FASTA Conversion ---
        fasta_dir = os.path.join(tmpDir, "5_fasta_converted")
        final_seq_count = convert_to_fasta(
             qfiltered_dir, fasta_dir, sampleids,
             options.base_phred_quality_score, options.FASTX_fastq_to_fasta_keep_n,
             options.verbose, summary_file_handle
        )
        # convert_to_fasta handles empty inputs and logs warnings/errors

        # --- Step 5: Merging ---
        merge_fasta(
            fasta_dir, options.outputdir, sampleids,
            summary_file_handle=summary_file_handle, final_output_name="prepped.fasta"
        )
        # merge_fasta handles missing/empty inputs

    except (subprocess.CalledProcessError, FileNotFoundError, ValueError, Exception) as e:
        # Catch any critical error during the pipeline steps
        logging.critical(f"Pipeline execution failed: {type(e).__name__} - {e}")
        # Log traceback if verbose for debugging
        if options.verbose:
            logging.debug(traceback.format_exc())
        pipeline_failed = True # Mark pipeline as failed
    # === End Pipeline Execution ===

    # --- Clean Up Intermediate Files ---
    if not options.retain and os.path.exists(tmpDir):
        logging.info("Cleaning temporary directory...")
        try:
            shutil.rmtree(tmpDir)
            logging.info(f"{Colors.BLUE}... Temporary directory removed.{Colors.RESET}")
        except OSError as e:
            # Log warning if cleanup fails, but don't exit pipeline
            logging.warning(f"Could not remove temporary directory {tmpDir}: {e}")
    elif options.retain:
         logging.info(f"Intermediate files retained in: {tmpDir}")


    # --- Final Summary & Timing ---
    end_time = time.time()
    try: # Ensure start_time exists before calculating duration
        duration = end_time - start_time
        duration_str = f"{duration:.2f} seconds"
    except NameError:
         duration_str = "N/A (calculation failed)"

    final_status = "FAILED" if pipeline_failed else "completed successfully"
    final_result_file = os.path.join(options.outputdir, "prepped.fasta")

    logging.info("-" * 30)
    logging.info(f"pipits_prep {final_status}.")
    # Check final output file existence and content only if pipeline didn't fail
    if not pipeline_failed:
         if os.path.exists(final_result_file):
             # Verify final count matches calculated final_seq_count
             final_check_lines = get_file_line_count(final_result_file)
             final_check_seqs = 0
             if final_check_lines > 0 and final_check_lines % 2 == 0:
                 final_check_seqs = final_check_lines // 2
             elif final_check_lines > 0:
                 logging.warning(f"Final output file {final_result_file} has unexpected line count ({final_check_lines}).")

             logging.info(f"Final output file: {final_result_file} ({final_check_seqs} sequences)")
             # Add warning if final file is empty or count mismatch
             if final_check_seqs == 0 and final_seq_count > 0:
                  logging.warning(f"Pipeline seemed successful (calculated {final_seq_count} sequences), but the final prepped.fasta file contains 0 sequences.")
             elif final_check_seqs != final_seq_count:
                  logging.warning(f"Mismatch between final calculated sequence count ({final_seq_count}) and count based on final file ({final_check_seqs}).")

         else:
             # This case means pipeline completed but final file is missing
             logging.error(f"Pipeline completed but final output file not found: {final_result_file}")
             pipeline_failed = True # Mark as failed if final output is missing
             final_status = "FAILED (Missing Output)"

    logging.info(f"Total processing time: {duration_str}")
    logging.info("-" * 30)


    # --- Write Final Summary Info & Close Handle ---
    if summary_file_handle:
        try:
            summary_file_handle.write("# --- End Summary ---\n")
            summary_file_handle.write(f"End time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            summary_file_handle.write(f"Duration: {duration_str}\n")
            summary_file_handle.write(f"Status: {final_status.upper()}\n")
            summary_file_handle.close() # Close the file handle
        except Exception as e:
            logging.error(f"Failed to write final summary or close summary file: {e}")

    # Return appropriate exit code
    return 1 if pipeline_failed else 0


# *** Entry point function for setup.py ***
def main():
    """Parses arguments and runs the PREP pipeline."""
    parser = argparse.ArgumentParser(
        description="pipits_prep: Reindex, join (VSEARCH), quality filter, convert and merge Illumina FASTQ data.")
    # Define arguments using parser...
    parser.add_argument(
        "-i", action="store", dest="dataDir", metavar="<DIR>",
        help="[REQUIRED] Directory containing raw sequence files (FASTQ format, can be .gz or .bz2 compressed).", required=True)
    parser.add_argument(
        "-o", action="store", dest="outputdir", metavar="<DIR>", default="pipits_prep_output",
        help="Directory to output results [default: pipits_prep_output]", required=False)
    parser.add_argument(
        "-l", action="store", dest="listfile", metavar="<FILE>",
        help="[REQUIRED] Tab-separated file with three columns: SampleID, ForwardReadFilename, ReverseReadFilename. Only files listed here will be processed.", required=True)
    parser.add_argument(
        "--FASTX-q", action="store", dest="FASTX_fastq_quality_filter_q", metavar="<INT>", type=int, default=30,
        help="FASTX quality filter: Minimum quality score to keep [default: 30]", required=False)
    parser.add_argument(
        "--FASTX-p", action="store", dest="FASTX_fastq_quality_filter_p", metavar="<INT>", type=int, default=80,
        help="FASTX quality filter: Minimum percent of bases that must have '--FASTX-q' quality [default: 80]", required=False)
    parser.add_argument(
        "--keep-Ns", action="store_true", dest="FASTX_fastq_to_fasta_keep_n", default=False,
        help="Keep sequences containing 'N' nucleotides during FASTQ to FASTA conversion [default: False, sequences with N are removed by FASTX tools]", required=False)
    parser.add_argument(
        "-b", action="store", dest="base_phred_quality_score", metavar="<INT>", type=int, default=33,
        help="Base PHRED quality score encoding of input FASTQ files (e.g., 33 or 64) [default: 33]", required=False)
    parser.add_argument(
        "-r", "--retain", action="store_true", dest="retain", default=False,
        help="Retain intermediate files in the 'tmp' subdirectory.", required=False)
    parser.add_argument(
        "-v", "--verbose", action="store_true", dest="verbose", default=False,
        help="Verbose mode (log debug messages to console and file).", required=False)
    parser.add_argument(
        "-t", "--threads", action="store", dest="threads", metavar="<INT>", type=int, default=1,
        help="Number of threads for VSEARCH joining [default: 1]", required=False)
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    options = parser.parse_args()

    # --- Setup Logging ---
    log_file_path = os.path.join(options.outputdir, "pipits_prep.log")
    # Create output dir early if possible
    try:
        if not os.path.exists(options.outputdir): os.makedirs(options.outputdir)
    except OSError as e:
        # Use print for this early error as logging might not be set up
        print(f"CRITICAL: Could not create output directory {options.outputdir}: {e}", file=sys.stderr)
        sys.exit(1) # Exit early
    setup_logging(log_file_path, options.verbose)

    # --- Run Pipeline ---
    exit_code = run_prep_pipeline(options)
    sys.exit(exit_code) # Exit with the code returned by the pipeline function


# --- Main Execution Block ---
if __name__ == '__main__':
    main()
