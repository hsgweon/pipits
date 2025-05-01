#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import shutil
import multiprocessing
import logging
import time
import glob
from itertools import islice
# *** ADDED traceback import for detailed error logging ***
import traceback

# --- Configuration ---
__version__ = "4.0"
__author__ = "Hyun Soon Gweon"
__copyright__ = "Copyright 2015-2024, The PIPITS Project"
__license__ = "GPL"
__maintainer__ = "Hyun Soon Gweon"
__email__ = "h.s.gweon@reading.ac.uk"

VSEARCH = "vsearch"
ITSX = "ITSx"
SEQKIT = "seqkit"

class Colors:
    RESET = '\033[0m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'

# --- Custom Logger Formatter with Color ---
class ColorFormatter(logging.Formatter):
    LOG_COLORS = {
        logging.DEBUG: Colors.CYAN,
        logging.INFO: Colors.BLUE,
        logging.WARNING: Colors.YELLOW,
        logging.ERROR: Colors.RED,
        logging.CRITICAL: Colors.MAGENTA
    }

    def format(self, record):
        log_message = super().format(record)
        color = self.LOG_COLORS.get(record.levelno, Colors.RESET)
        handler = getattr(self, 'handler', None)
        is_console_tty = handler and hasattr(handler, 'stream') and handler.stream.isatty()

        if is_console_tty:
             return f"{color}{log_message}{Colors.RESET}"
        else:
             # Strip color codes if not a TTY
             log_message = log_message.replace(Colors.RED, '').replace(Colors.GREEN, '').replace(Colors.YELLOW, '') \
                                    .replace(Colors.BLUE, '').replace(Colors.MAGENTA, '').replace(Colors.CYAN, '') \
                                    .replace(Colors.RESET, '')
             return log_message


# --- Helper Functions ---

def setup_logging(log_file_path, verbose):
    """Sets up logging to file and colorful console output."""
    log_format_file = '%(asctime)s - %(levelname)s - %(message)s'
    log_format_console = '%(levelname)s: %(message)s'
    log_level = logging.DEBUG if verbose else logging.INFO

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG) # Capture all levels

    # Clear existing handlers to prevent duplicates if called multiple times
    if logger.hasHandlers():
        logger.handlers.clear()

    # File handler (always logs DEBUG level and up)
    file_handler = logging.FileHandler(log_file_path, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter(log_format_file)
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)

    # Console handler (logs INFO or DEBUG based on verbose flag)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level) # Set level based on verbosity
    color_formatter = ColorFormatter(log_format_console)
    color_formatter.handler = console_handler # Attach handler for TTY check
    console_handler.setFormatter(color_formatter)
    logger.addHandler(console_handler)


def run_cmd(command, log_stdout=True):
    """Runs a command using subprocess and logs."""
    logging.debug(f"Running command: {command}")
    try:
        is_shell = isinstance(command, str) # Check if command is a string
        process = subprocess.run(command, shell=is_shell, check=True,
                                 capture_output=True, text=True, encoding='utf-8', errors='ignore')
        # Log output only if DEBUG is enabled
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            stdout_output = process.stdout.strip() if process.stdout else None
            stderr_output = process.stderr.strip() if process.stderr else None
            if log_stdout and stdout_output:
                logging.debug(f"Command STDOUT:\n---\n{stdout_output}\n---")
            if stderr_output:
                # Log stderr even if not explicitly requested in debug mode
                logging.debug(f"Command STDERR:\n---\n{stderr_output}\n---")
        return process
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}")
        # Handle command format whether list or string
        failed_cmd_str = ' '.join(e.cmd) if isinstance(e.cmd, list) else e.cmd
        logging.error(f"Failed command: {failed_cmd_str}")
        stderr_output = e.stderr.strip() if e.stderr else "N/A"
        stdout_output = e.stdout.strip() if e.stdout else "N/A"
        logging.error(f"STDERR:\n---\n{stderr_output}\n---")
        # Log stdout on error only if it contains something
        if stdout_output != "N/A": logging.error(f"STDOUT:\n---\n{stdout_output}\n---")
        raise # Re-raise the exception to be caught by the main loop
    except FileNotFoundError:
         # Extract command name correctly
         cmd_name = command.split()[0] if isinstance(command, str) else command[0]
         logging.error(f"Command not found: '{cmd_name}'. Please ensure it's installed and in your PATH.")
         raise # Re-raise
    except Exception as e:
        # Catch other potential errors during execution
        failed_cmd_str = ' '.join(command) if isinstance(command, list) else command
        logging.error(f"An unexpected error occurred running command: {failed_cmd_str}")
        logging.error(f"Error details: {e}")
        raise # Re-raise


def blocks(files, size=65536):
    """Helper function to read files in chunks."""
    while True:
        b = files.read(size)
        if not b: break
        yield b

def get_fasta_count(fasta_file):
    """Counts sequences in a FASTA file."""
    count = 0
    try:
        # Use 'rb' to read bytes and count '>' bytes for efficiency
        with open(fasta_file, 'rb') as f:
            count = sum(bl.count(b'>') for bl in blocks(f))
    except FileNotFoundError:
        logging.warning(f"File not found for counting: {fasta_file}")
        return 0 # Return 0 if file doesn't exist
    except Exception as e:
         logging.error(f"Error counting FASTA sequences in {fasta_file}: {e}")
         return 0 # Return 0 on error
    return count

def parse_uc_file(uc_file):
    """Parses a VSEARCH .uc file, stripping ;size=N annotations."""
    rep_to_originals = {}
    original_to_rep = {}
    representative_set = set()

    logging.info(f"Parsing UC file: {uc_file}")
    try:
        with open(uc_file, 'r', encoding='utf-8', errors='ignore') as f:
            line_count = 0
            for line in f:
                line_count += 1
                if line_count % 500000 == 0: # Log progress periodically
                     logging.debug(f"... parsed {line_count} lines of UC file")

                if line.startswith('#'): # Skip comment lines
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 10: # Basic validation
                     logging.debug(f"Skipping malformed UC line {line_count}: {line.strip()}")
                     continue

                record_type = fields[0]
                # Strip ;size=N from query and target IDs
                query_id_full = fields[8]
                target_id_full = fields[9]
                query_id = query_id_full.split(';', 1)[0] # Get part before first ';'
                target_id = target_id_full.split(';', 1)[0]


                if record_type == 'S': # Seed (representative) line
                    rep_id = query_id # Representative is the query for S lines
                    if rep_id not in rep_to_originals:
                         rep_to_originals[rep_id] = []
                    rep_to_originals[rep_id].append(rep_id) # Add itself to its list
                    original_to_rep[rep_id] = rep_id     # Map itself
                    representative_set.add(rep_id)
                elif record_type == 'H': # Hit line
                    original_id = query_id # Original is the query for H lines
                    rep_id = target_id    # Representative is the target for H lines
                    if rep_id not in rep_to_originals:
                         # This *shouldn't* happen if the UC file is correctly formed
                         # (i.e., all H targets should have appeared as S queries)
                         logging.warning(f"Representative '{rep_id}' for hit '{original_id}' not seen in S lines (line {line_count}). Creating entry anyway.")
                         rep_to_originals[rep_id] = []
                         representative_set.add(rep_id) # Track this rep as well
                    rep_to_originals[rep_id].append(original_id) # Add original to rep's list
                    original_to_rep[original_id] = rep_id     # Map original to rep


    except FileNotFoundError:
         logging.error(f"UC file not found: {uc_file}")
         raise # Re-raise to be caught by main loop
    except Exception as e:
         logging.error(f"Error parsing UC file {uc_file}: {e}")
         raise # Re-raise

    logging.info(f"Finished parsing UC file. Found {len(representative_set)} unique representative IDs (clean).")
    logging.debug(f"Size of rep_to_originals map: {len(rep_to_originals)}")
    return rep_to_originals, original_to_rep

def simple_fasta_iterator(fasta_file):
    """Yields sequence ID (stripped of size) and sequence from a FASTA file."""
    try:
        with open(fasta_file, 'r', encoding='utf-8', errors='ignore') as f:
            seq_id = None
            sequence_parts = []
            for line in f:
                line = line.strip()
                if not line: # Skip empty lines
                    continue
                if line.startswith('>'):
                    # If we were processing a sequence, yield it first
                    if seq_id is not None:
                        yield seq_id, "".join(sequence_parts)
                    # Start new sequence
                    full_header = line[1:]
                    seq_id = full_header.split(';', 1)[0] # Strip > and ;size=N annotation
                    sequence_parts = []
                elif seq_id is not None: # Append sequence lines only if we have a current ID
                    sequence_parts.append(line)
            # Yield the last sequence in the file
            if seq_id is not None:
                yield seq_id, "".join(sequence_parts)
    except FileNotFoundError:
        logging.error(f"FASTA file not found for iteration: {fasta_file}")
        raise # Re-raise
    except Exception as e:
        logging.error(f"Error iterating through FASTA file {fasta_file}: {e}")
        raise # Re-raise


# *** WRAPPER FUNCTION for the main logic ***
def run_funits_pipeline(options):
    """Contains the main execution logic called by main()."""

    if options.threads < 1:
        logging.error("Number of parallel processes (-t) must be at least 1.")
        return 1 # Return error code
    if options.itsx_threads < 1:
        logging.error("Number of threads per ITSx process (--itsx-threads) must be at least 1.")
        return 1 # Return error code

    min_len_threshold = 100

    # --- Setup ---
    if not os.path.exists(options.outDir):
        os.makedirs(options.outDir)

    tmpDir = os.path.join(options.outDir, "intermediate")
    if not options.retain and os.path.exists(tmpDir):
         logging.info(f"Removing existing intermediate directory: {tmpDir}")
         shutil.rmtree(tmpDir)
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir)

    summary_file = os.path.join(options.outDir, "summary.log")
    version_file = os.path.join(options.outDir, "versions.log")

    summary_lines = []

    logging.info(f"--- Starting pipits_funits v{__version__} ---")
    logging.info(f"Output directory: {os.path.abspath(options.outDir)}")
    logging.info(f"Input file: {os.path.abspath(options.input)}")
    logging.info(f"ITS region: {options.ITSx_subregion}")
    logging.info(f"Parallel ITSx processes: {options.threads}")
    logging.info(f"Threads per ITSx process: {options.itsx_threads}")
    logging.info(f"Retain intermediates: {options.retain}")
    logging.info(f"Verbose logging: {options.verbose}")
    logging.info("-" * 30)


    # --- Log versions ---
    try:
        version_info = f"pipits_funits: {__version__}\n"
        for tool, tool_cmd in [("VSEARCH", f"{VSEARCH} --version"),
                               ("ITSx", f"{ITSX} -v"),
                               ("SEQKIT", f"{SEQKIT} version")]:
             try:
                 # Use subprocess.run for consistency and better error handling
                 process = subprocess.run(tool_cmd, shell=True, check=True, capture_output=True, text=True, timeout=10)
                 output = process.stdout.strip() or process.stderr.strip()
                 version_info += f"{tool}: {output.splitlines()[0]}\n" # Take first line
             except Exception as tool_e:
                 version_info += f"{tool}: Not found or error getting version ({type(tool_e).__name__})\n"
                 logging.warning(f"Could not get version for {tool}: {tool_e}")

        with open(version_file, "w") as vf:
             vf.write(version_info)
        logging.info(f"Tool versions logged to {version_file}")
    except Exception as e:
        logging.warning(f"Could not log tool versions: {e}")

    # --- Start Processing ---
    start_time = time.time()

    try: # Wrap main logic in try/except
        # --- Input Validation ---
        if not os.path.exists(options.input) or os.stat(options.input).st_size == 0:
            logging.critical(f"Input file not found or is empty: {options.input}")
            return 1 # Exit with error code

        logging.info("Checking input FASTA for spaces in headers...")
        # Header check simplified - assume simple_fasta_iterator handles ID extraction correctly
        logging.info("... header check skipped (IDs assumed handled by iterators).")


        initial_seq_count = get_fasta_count(options.input)
        logging.info(f"Initial sequence count: {initial_seq_count}")
        summary_lines.append(f"Number of input sequences: {initial_seq_count}")
        if initial_seq_count == 0:
            logging.critical("Input file contains 0 sequences based on '>' count.")
            return 1 # Exit with error code

        # --- 1. Dereplicate ---
        derep_fasta = os.path.join(tmpDir, "derep.fasta")
        derep_uc = os.path.join(tmpDir, "derep.uc")
        logging.info("Dereplicating sequences [VSEARCH]...")
        cmd_derep = f"{VSEARCH} --derep_fulllength \"{options.input}\" --output \"{derep_fasta}\" --uc \"{derep_uc}\" --fasta_width 0 --sizeout"
        run_cmd(cmd_derep) # Will raise error if fails

        derep_seq_count = get_fasta_count(derep_fasta)
        logging.info(f"Dereplicated sequence count: {derep_seq_count}")
        summary_lines.append(f"Number of dereplicated sequences: {derep_seq_count}")
        if derep_seq_count == 0:
             logging.error("Dereplication resulted in 0 sequences.")
             return 1 # Exit with error code


        # --- 2. Run ITSx ---
        itsx_processed_fasta = os.path.join(tmpDir, f"derep.{options.ITSx_subregion}.fasta")
        itsx_output_prefix_base = 'derep' # Base name for ITSx output files

        if options.threads > 1:
            # --- Parallel ITSx ---
            logging.info(f"Splitting dereplicated file for parallel ITSx (Processes: {options.threads}) [SEQKIT]...")
            split_dir = os.path.join(tmpDir, "split")
            if os.path.exists(split_dir): shutil.rmtree(split_dir) # Clean slate
            os.makedirs(split_dir)

            # Use seqkit to split the dereplicated FASTA into parts
            cmd_split = f"{SEQKIT} split \"{derep_fasta}\" -p {options.threads} -f --out-dir \"{split_dir}\" -w 0"
            run_cmd(cmd_split) # Will raise error if fails

            # Find the created split files
            split_files = sorted(glob.glob(os.path.join(split_dir, f"{itsx_output_prefix_base}.part_*.fasta")))
            if not split_files:
                 logging.critical("Failed to find split files after running seqkit.")
                 return 1 # Exit with error code
            if len(split_files) != options.threads:
                 logging.warning(f"Expected {options.threads} split files, but found {len(split_files)}. Proceeding with found files.")

            logging.info(f"Split into {len(split_files)} parts.")

            logging.info(f"Extracting {options.ITSx_subregion} using ITSx in parallel (Processes={options.threads}, Threads/Process={options.itsx_threads}) [ITSx]...")
            itsx_cmds = []
            for part_file in split_files:
                part_base_name = os.path.basename(part_file)
                # Define output prefix for this part *within* the split directory
                part_output_prefix = os.path.join(split_dir, os.path.splitext(part_base_name)[0])
                cmd = (f"{ITSX} -i \"{part_file}\" -o \"{part_output_prefix}\" --preserve T -t F --cpu {options.itsx_threads} "
                       f"--silent T --save_regions {options.ITSx_subregion}")
                itsx_cmds.append(cmd)

            # --- Execute ITSx commands in parallel ---
            processes = {}
            process_errors = {} # Store stderr for failed processes
            try:
                # Launch processes
                for i, cmd in enumerate(itsx_cmds):
                     logging.debug(f"Launching ITSx task {i+1}/{len(itsx_cmds)}: {cmd}")
                     # Use Popen for parallel execution
                     processes[i] = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8', errors='ignore')

                # Monitor processes
                completed_count = 0
                total_tasks = len(processes)
                while processes:
                    finished_indices = []
                    for i, process in processes.items():
                         if process.poll() is not None: # Process finished
                            completed_count += 1
                            logging.info(f"ITSx task {completed_count}/{total_tasks} completed [ITSx].")
                            if process.returncode != 0:
                                 logging.warning(f"ITSx subprocess failed (Index: {i}, Return Code: {process.returncode}) [ITSx]")
                                 stderr_output = process.stderr.read() # Capture stderr on failure
                                 process_errors[i] = stderr_output
                                 if stderr_output: logging.warning(f"ITSx STDERR (Index {i}):\n---\n{stderr_output}\n---")
                            else:
                                # Log stderr even on success if verbose
                                if options.verbose:
                                     stderr_output = process.stderr.read()
                                     if stderr_output: logging.debug(f"ITSx STDERR (Index {i}, Success):\n---\n{stderr_output}\n---")
                            finished_indices.append(i) # Mark for removal

                    # Remove finished processes from the monitoring dictionary
                    for i in finished_indices:
                         del processes[i]

                    # Small sleep to avoid busy-waiting
                    time.sleep(0.2)

                if process_errors:
                    logging.error("One or more parallel ITSx tasks failed. See warnings above.")
                    # Optionally, stop the pipeline here if any task fails severely
                    # return 1

            except Exception as e:
                 logging.critical(f"Error occurred during parallel ITSx execution management: {e}")
                 # Attempt to terminate any remaining processes
                 for i, process in processes.items():
                      logging.warning(f"Terminating potentially running ITSx process (Index {i}).")
                      try: process.terminate()
                      except Exception: pass # Ignore errors during termination
                 return 1 # Exit with error code

            # --- Concatenate results ---
            logging.info("Concatenating parallel ITSx results...")
            # Construct the expected output file pattern
            itsx_split_outputs = sorted(glob.glob(os.path.join(split_dir, f"{itsx_output_prefix_base}.part_*.{options.ITSx_subregion}.fasta")))

            if not itsx_split_outputs:
                 logging.warning(f"No ITSx output files found matching pattern: {split_dir}/{itsx_output_prefix_base}.part_*.{options.ITSx_subregion}.fasta")
                 # Create an empty file to avoid downstream errors
                 open(itsx_processed_fasta, 'w').close()
            else:
                 try:
                     with open(itsx_processed_fasta, 'wb') as outfile: # Open in binary write mode
                         for infile_path in itsx_split_outputs:
                             try:
                                 # Open each part in binary read mode
                                 with open(infile_path, 'rb') as infile:
                                     shutil.copyfileobj(infile, outfile) # Efficiently copy file object
                             except FileNotFoundError:
                                 # Log if a specific part is missing (might happen if a process failed)
                                 logging.warning(f"ITSx output file not found during concatenation: {infile_path}")
                             except Exception as copy_e:
                                  logging.error(f"Error copying file {infile_path}: {copy_e}")
                 except Exception as concat_e:
                      logging.error(f"Error concatenating ITSx results to {itsx_processed_fasta}: {concat_e}")
                      logging.error("Concatenation failed. Result file may be incomplete.")
                      # Consider exiting if concatenation fails critically
                      # return 1


        else: # --- Single thread ITSx process run ---
            logging.info(f"Extracting {options.ITSx_subregion} using ITSx (Processes=1, Threads/Process={options.itsx_threads}) [ITSx]...")
            # Define the output prefix directly in the intermediate directory
            itsx_output_prefix_single = os.path.join(tmpDir, itsx_output_prefix_base)
            cmd_itsx = (f"{ITSX} -i \"{derep_fasta}\" -o \"{itsx_output_prefix_single}\" --preserve T -t F --cpu {options.itsx_threads} "
                       f"--silent T --save_regions {options.ITSx_subregion}")
            run_cmd(cmd_itsx) # Will raise error if fails

            # Define the expected output file name based on the prefix and region
            expected_itsx_output_single = f"{itsx_output_prefix_single}.{options.ITSx_subregion}.fasta"

            # Move or copy the single output file to the standardized name
            if os.path.exists(expected_itsx_output_single):
                try:
                    shutil.move(expected_itsx_output_single, itsx_processed_fasta)
                    logging.debug(f"Moved ITSx output {expected_itsx_output_single} to {itsx_processed_fasta}")
                except Exception as move_e:
                    logging.warning(f"Could not move ITSx output file, attempting copy: {move_e}")
                    try:
                        shutil.copy2(expected_itsx_output_single, itsx_processed_fasta)
                        logging.debug(f"Copied ITSx output {expected_itsx_output_single} to {itsx_processed_fasta}")
                    except Exception as copy_e:
                        logging.error(f"Failed to copy ITSx output file {expected_itsx_output_single} to {itsx_processed_fasta}: {copy_e}")
                        # Create empty file if copy fails too
                        open(itsx_processed_fasta, 'w').close()
            else:
                # If the expected output doesn't exist, log warning and create empty file
                logging.warning(f"Expected ITSx output file not found after single-threaded run: {expected_itsx_output_single}")
                open(itsx_processed_fasta, 'w').close()


        itsx_seq_count = get_fasta_count(itsx_processed_fasta)
        logging.info(f"ITS sequences extracted (dereplicated): {itsx_seq_count}")
        summary_lines.append(f"Number of ITS sequences (dereplicated): {itsx_seq_count}")
        if itsx_seq_count == 0:
            logging.warning(f"No sequences identified as {options.ITSx_subregion} by ITSx.")
            # Don't exit here, allow length filtering to proceed (it will handle empty input)


        # --- 3. Filter by Length ---
        length_filtered_fasta = os.path.join(tmpDir, f"derep.{options.ITSx_subregion}.sizefiltered.fasta")
        logging.info(f"Filtering sequences shorter than {min_len_threshold}bp [VSEARCH]...")
        cmd_filter = (f"{VSEARCH} --fastx_filter \"{itsx_processed_fasta}\" --fastaout \"{length_filtered_fasta}\" "
                     f"--fastq_minlen {min_len_threshold} --fasta_width 0")

        # Only run filter if ITSx output exists and is not empty
        if os.path.exists(itsx_processed_fasta) and os.stat(itsx_processed_fasta).st_size > 0:
            run_cmd(cmd_filter) # Will raise error if fails
        else:
            logging.warning(f"Skipping length filtering as input file is missing or empty: {itsx_processed_fasta}")
            # Create empty output file for consistency
            open(length_filtered_fasta, 'w').close()


        filtered_seq_count = get_fasta_count(length_filtered_fasta)
        logging.info(f"Length-filtered sequences (dereplicated): {filtered_seq_count}")
        summary_lines.append(f"Number of length-filtered sequences (dereplicated): {filtered_seq_count}")
        if filtered_seq_count == 0 and itsx_seq_count > 0:
            logging.warning(f"Length filtering removed all {itsx_seq_count} ITSx sequences (threshold {min_len_threshold}bp).")
        # Don't exit if filtering yields 0, proceed to re-inflation


        # --- 4. Load Passed Representatives and UC Map ---
        logging.info("Loading dereplication map and filtered ITS sequences...")
        passed_reps = {} # Dictionary to store {rep_id: sequence} for passed reps
        try:
            # Only try loading if the filtered file exists and is not empty
            if os.path.exists(length_filtered_fasta) and os.stat(length_filtered_fasta).st_size > 0:
                 # Use the iterator that strips ;size=N if present (VSEARCH filter output shouldn't have it)
                 for rep_id, sequence in simple_fasta_iterator(length_filtered_fasta):
                    passed_reps[rep_id] = sequence # Store with clean ID
            else:
                 logging.warning(f"Length filtered file is empty or missing: {length_filtered_fasta}. No representatives passed filters.")
        except Exception as e:
             logging.error(f"Failed to load filtered sequences from {length_filtered_fasta}: {e}")
             return 1 # Exit with error code
        logging.info(f"Loaded {len(passed_reps)} representative sequences that passed filters.")

        try:
            # Use the modified parse_uc_file which strips ;size=N
            rep_to_originals, original_to_rep = parse_uc_file(derep_uc)
        except Exception as e:
             logging.error(f"Failed to parse UC file: {e}")
             # No need to raise here, already logged, just exit
             return 1 # Exit with error code

        # --- 5. Generate Final Output (Optimized Approach) ---
        final_output_fasta = os.path.join(options.outDir, "ITS.fasta")
        logging.info(f"Generating final output FASTA (Optimized Method): {final_output_fasta}")

        output_count = 0
        itsx_passed_output_count = 0
        original_kept_output_count = 0
        processed_originals = set() # Keep track of originals already written

        try:
            with open(final_output_fasta, 'w') as f_out:
                # --- Part A: Write sequences whose representatives passed filters ---
                logging.info("Writing sequences whose representatives passed filters...")
                reps_processed_count = 0
                # Iterate through the representatives that passed ITSx/length filters (keys are clean IDs)
                for rep_id, its_sequence in passed_reps.items():
                    reps_processed_count += 1
                    if reps_processed_count % 1000 == 0: # Log progress
                        logging.debug(f"... processing passed representative {reps_processed_count}/{len(passed_reps)}")

                    # Find all original sequences that map to this representative (using clean rep_id)
                    original_ids_for_rep = rep_to_originals.get(rep_id, [])

                    # If a passed rep has no originals (shouldn't happen with correct parsing), log and skip
                    if not original_ids_for_rep:
                         logging.warning(f"Representative '{rep_id}' passed filters but has no originals mapped in UC file? Skipping.")
                         continue

                    # Write the *extracted ITS sequence* for *each* original ID mapped to this rep
                    for original_id in original_ids_for_rep:
                        if original_id not in processed_originals: # Avoid duplicates if UC has strange entries
                            f_out.write(f">{original_id}\n")
                            f_out.write(f"{its_sequence}\n")
                            output_count += 1
                            itsx_passed_output_count += 1
                            processed_originals.add(original_id)

                logging.info(f"Finished writing {itsx_passed_output_count} sequences derived from passed representatives.")

                # --- Part B: Write original sequences whose representatives failed ---
                logging.info("Processing original sequences whose representatives failed filters...")
                failed_written_count = 0
                original_checked_count = 0
                # Iterate through the *original* input file sequences
                for original_id, original_sequence in simple_fasta_iterator(options.input):
                    original_checked_count += 1
                    if original_checked_count % 50000 == 0: # Log progress
                         logging.debug(f"... checked {original_checked_count}/{initial_seq_count} original sequences for failures")

                    # Skip if this original was already processed via its passed representative
                    if original_id in processed_originals:
                        continue

                    # Find its representative (original_to_rep maps clean original -> clean rep)
                    representative_id = original_to_rep.get(original_id)

                    # Determine if the representative FAILED the ITSx/length filters
                    # Conditions:
                    # 1. The original had a representative (representative_id is not None)
                    # 2. That representative is NOT in the set of passed representatives (passed_reps)
                    if representative_id and representative_id not in passed_reps:
                         # Write the *original* full sequence back
                         f_out.write(f">{original_id}\n")
                         f_out.write(f"{original_sequence}\n")
                         output_count += 1
                         original_kept_output_count += 1
                         failed_written_count += 1
                         processed_originals.add(original_id) # Mark as processed
                    elif not representative_id:
                         # This case means the original_id was likely filtered *before* dereplication
                         # or it was a singleton representative itself that failed ITSx/length.
                         # We should write the original back.
                         logging.debug(f"Original sequence '{original_id}' not found in UC map or its representative failed. Writing original.")
                         f_out.write(f">{original_id}\n")
                         f_out.write(f"{original_sequence}\n")
                         output_count += 1
                         original_kept_output_count += 1
                         failed_written_count += 1
                         processed_originals.add(original_id) # Mark as processed


                logging.info(f"Finished processing failed sequences. Wrote back {failed_written_count} original sequences.")

                # Final count sanity check
                if output_count != initial_seq_count:
                     logging.warning(f"Final output count ({output_count}) does not match initial count ({initial_seq_count}). Processed originals set size: {len(processed_originals)}. Check logs for details.")
                else:
                     logging.info("Final output count matches initial input count.")


        except Exception as e:
             logging.critical(f"Error generating final output: {e}")
             if options.verbose:
                  logging.debug(traceback.format_exc()) # Log full traceback if verbose
             return 1 # Exit with error code

        logging.info("Final output generation complete.")
        logging.info(f"Total sequences written to final output: {Colors.GREEN}{output_count}{Colors.RESET}")
        logging.info(f"Sequences derived from ITSx+filter pass: {Colors.GREEN}{itsx_passed_output_count}{Colors.RESET}")
        logging.info(f"Original sequences kept (no ITSx hit/filter fail): {Colors.YELLOW}{original_kept_output_count}{Colors.RESET}")

        summary_lines.append(f"Total number of sequences with ITS subregion (re-inflated): {itsx_passed_output_count}")
        summary_lines.append(f"Total number of sequences without conserved region / failing filters: {original_kept_output_count}")
        summary_lines.append(f"Total number of sequences in the final output (FUNITS): {output_count}")


        # --- Cleanup ---
        if not options.retain:
            logging.info("Cleaning temporary directory...")
            try:
                shutil.rmtree(tmpDir)
                logging.info("... done.")
            except OSError as e:
                logging.warning(f"Could not remove temporary directory {tmpDir}: {e}")
        else:
             logging.info(f"Intermediate files retained in: {tmpDir}")

        # --- Finish ---
        end_time = time.time()
        total_time = end_time - start_time
        logging.info(f"pipits_funits finished successfully in {Colors.GREEN}{total_time:.2f} seconds.{Colors.RESET}")
        logging.info(f"Final ITS sequences are in: {final_output_fasta}")
        logging.info(f"Next step: pipits_process [ Example: pipits_process -i \"{final_output_fasta}\" -o pipits_process ]")

        # --- Write Summary ---
        try:
             with open(summary_file, "w") as sf:
                  sf.write("--- PIPITS_FUNITS Summary ---\n")
                  for line in summary_lines:
                       sf.write(line + "\n")
                  sf.write(f"Processing completed in {total_time:.2f} seconds.\n")
        except Exception as e:
             logging.warning(f"Could not write summary file {summary_file}: {e}")

        return 0 # Success exit code

    # --- Catch errors from main logic ---
    except (subprocess.CalledProcessError, FileNotFoundError, IOError, ImportError, RuntimeError, ValueError) as e:
         logging.critical(f"A critical error occurred during processing: {type(e).__name__} - {e}")
         if options.verbose: logging.debug(traceback.format_exc())
         return 1 # Failure exit code
    except Exception as e:
         logging.critical(f"An unexpected critical error occurred: {type(e).__name__} - {e}")
         if options.verbose: logging.debug(traceback.format_exc())
         return 1 # Failure exit code


# *** Entry point function for setup.py ***
def main():
    """Parses arguments and runs the FUNITS pipeline."""
    parser = argparse.ArgumentParser(description=f"PIPITS_FUNITS v{__version__}: Extract ITS1 or ITS2 regions.")
    parser.add_argument("-i", dest="input", metavar="<FILE>", help="[REQUIRED] Input FASTA file (e.g., from pipits_prep).", required=True)
    parser.add_argument("-o", dest="outDir", metavar="<DIR>", help="[REQUIRED] Directory to output results.", required=True)
    parser.add_argument("-x", dest="ITSx_subregion", help="[REQUIRED] Subregion (ITS1 or ITS2).", required=True, choices=["ITS1", "ITS2"])
    parser.add_argument("-t", dest="threads", metavar="<INT>", help="Number of parallel ITSx *processes* to run [default: 1].", default=1, type=int, required=False)
    parser.add_argument("--itsx-threads", dest="itsx_threads", metavar="<INT>", help="Number of threads *per* ITSx process [default: 1]. Total threads used = -t * --itsx-threads.", default=1, type=int, required=False)
    parser.add_argument("-r", dest="retain", action="store_true", help="Retain intermediate files.", required=False)
    parser.add_argument("-v", dest="verbose", action="store_true", help="Verbose mode (more detailed logging).", required=False)
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    options = parser.parse_args()

    # Setup logging here, after args are parsed, so verbosity is known
    log_file = os.path.join(options.outDir, "pipits_funits.log")
    # Create output dir early if possible to ensure log file can be written
    try:
        if not os.path.exists(options.outDir): os.makedirs(options.outDir)
    except OSError as e:
        print(f"CRITICAL: Could not create output directory {options.outDir}: {e}", file=sys.stderr)
        sys.exit(1) # Exit early if output dir fails
    setup_logging(log_file, options.verbose)

    # Run the main pipeline logic and get exit code
    exit_code = run_funits_pipeline(options)
    sys.exit(exit_code) # Exit with the code from the pipeline function


# --- Main Execution Block ---
if __name__ == '__main__':
    main()
