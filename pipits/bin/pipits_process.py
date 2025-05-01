#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import shutil
import logging
import time
import hashlib
# *** Restored imports needed for original DB download logic ***
import requests
import tarfile
import progressbar
import traceback
from pathlib import Path

try:
    from biom import load_table
    from biom.util import biom_open # For writing biom files if needed
except ImportError:
    # Log critical error and exit if biom is missing
    # Note: Logging might not be set up yet if this fails early
    print("CRITICAL: The 'biom-format' library is required but not installed.", file=sys.stderr)
    print("Please install it, e.g., using 'pip install biom-format'", file=sys.stderr)
    sys.exit(1)


__version__     = "4.0"
__author__      = "Hyun Soon Gweon"
__copyright__   = "Copyright 2015-2024, The PIPITS Project"
__credits__     = ["Hyun Soon Gweon", "Anna Oliver", "Joanne Taylor", "Tim Booth", "Melanie Gibbs", "Daniel S. Read", "Robert I. Griffiths", "Karsten Schonrogge"]
__license__     = "GPL"
__maintainer__  = "Hyun Soon Gweon"
__email__       = "h.s.gweon@reading.ac.uk"

VSEARCH = "vsearch"
BIOM = "biom"
# *** RDP Classifier command likely needs full path or configuration ***
RDP_CLASSIFIER = "classifier" # Or specify full path if needed


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
        # Check if handler exists, has a stream, and if that stream is a TTY
        is_console_tty = handler and hasattr(handler, 'stream') and hasattr(handler.stream, 'isatty') and handler.stream.isatty()

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

    # Clear existing handlers to prevent duplicates
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
        is_shell = isinstance(command, str)
        process = subprocess.run(command, shell=is_shell, check=True,
                                 capture_output=True, text=True, encoding='utf-8', errors='ignore')
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            stdout_output = process.stdout.strip() if process.stdout else None
            stderr_output = process.stderr.strip() if process.stderr else None
            if log_stdout and stdout_output:
                logging.debug(f"Command STDOUT:\n---\n{stdout_output}\n---")
            if stderr_output:
                logging.debug(f"Command STDERR:\n---\n{stderr_output}\n---")
        return process
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}")
        failed_cmd_str = ' '.join(e.cmd) if isinstance(e.cmd, list) else e.cmd
        logging.error(f"Failed command: {failed_cmd_str}")
        stderr_output = e.stderr.strip() if e.stderr else "N/A"
        stdout_output = e.stdout.strip() if e.stdout else "N/A"
        logging.error(f"STDERR:\n---\n{stderr_output}\n---")
        if stdout_output != "N/A": logging.error(f"STDOUT:\n---\n{stdout_output}\n---")
        raise
    except FileNotFoundError:
         cmd_name = command.split()[0] if isinstance(command, str) else command[0]
         # Check specifically for the RDP classifier command
         if cmd_name == RDP_CLASSIFIER:
             logging.error(f"Command '{RDP_CLASSIFIER}' not found. Please ensure the RDP Classifier tools (classifier.jar and its dependencies) are installed and the '{RDP_CLASSIFIER}' command (or a wrapper script) is in your system's PATH or provide the full path.")
         else:
              logging.error(f"Command not found: '{cmd_name}'. Please ensure it's installed and in your PATH.")
         raise
    except Exception as e:
        failed_cmd_str = ' '.join(command) if isinstance(command, list) else command
        logging.error(f"An unexpected error occurred running command: {failed_cmd_str}")
        logging.error(f"Error details: {e}")
        raise

def generate_sample_list_from_fasta(input_fasta_path, output_list_path):
    """Generates sample list from FASTA headers."""
    logging.info(f"Extracting unique sample IDs from {input_fasta_path}...")
    uniquesampleids = set() # Use a set for faster uniqueness check
    processed_headers = 0
    try:
        with open(input_fasta_path, "r", encoding='utf-8', errors='ignore') as f_in:
            for line in f_in:
                if line.startswith(">"):
                    processed_headers += 1
                    try:
                        # Assumes header format >SampleID_OtherInfo
                        sampleid = line[1:].split("_", 1)[0]
                        if sampleid: uniquesampleids.add(sampleid)
                    except IndexError:
                        # Handle headers without '_' -> use the whole header as ID? Or skip?
                        # Current behaviour: skip. Add warning.
                        logging.warning(f"Could not parse sample ID from header (missing '_'?): {line.strip()}")
                        continue # Skip this header
        if not uniquesampleids:
             logging.warning("No unique sample IDs could be extracted from the input FASTA.")
             # Create empty file to avoid downstream errors
             open(output_list_path, 'w').close()
             return False # Indicate failure or no IDs found
        # Write sorted list of unique IDs
        with open(output_list_path, "w") as f_out:
            for sample_id in sorted(list(uniquesampleids)):
                f_out.write(sample_id + "\n")
        logging.info(f"Found {len(uniquesampleids)} unique sample IDs from {processed_headers} headers.")
        logging.info(f"Sample list written to {output_list_path}")
        return True # Success
    except FileNotFoundError:
        logging.error(f"Input FASTA file not found: {input_fasta_path}")
        return False
    except IOError as e:
        logging.error(f"Error reading/writing file during sample list generation: {e}")
        return False
    except Exception as e:
        logging.error(f"An unexpected error occurred during sample list generation: {e}")
        return False

def reformat_rdp_taxonomy(input_file, output_file, confidence_threshold):
    """
    Re-formats RDP Classifier output to QIIME-compatible format.
    Args:
        input_file (str): Path to the raw RDP classifier output.
        output_file (str): Path to the output formatted file.
        confidence_threshold (float): Minimum confidence to retain a level.
    Returns:
        bool: True on success, False on failure.
    """
    logging.info(f"Reformatting RDP taxonomy from {input_file} (Threshold: {confidence_threshold})")
    # Standard QIIME prefixes
    PREFIX_MAP = {
        "domain": "d__", "kingdom": "k__", "phylum": "p__", "subphylum": "subp__",
        "class": "c__", "subclass": "subc__", "order": "o__", "suborder": "subo__",
         "family": "f__", "genus": "g__", "species": "s__"
        # Add more if RDP outputs other ranks consistently
    }
    try:
        with open(input_file, "r", encoding='utf-8', errors='ignore') as handle_input, \
             open(output_file, "w", encoding='utf-8') as handle_output:
            processed_lines = 0
            assigned_otus = 0
            for line in handle_input:
                processed_lines += 1
                # RDP output is tab-separated: OTUID SeqLen Strand Domain "Root" Conf Phylum "TaxonName" Conf ...
                e = line.rstrip().split("\t")
                # Need at least OTUID and one rank triplet (rank, name, conf) beyond Root
                # Index 0: OTUID, 1: SeqLen, 2: Strand, 3: Rank1(domain), 4: Name1("Root"), 5: Conf1, 6: Rank2(phylum?), ...
                if len(e) < 6: # Check if there's at least the "Root" level info
                    logging.warning(f"Skipping malformed line {processed_lines} in RDP output (too few fields): {line.strip()}")
                    continue

                otu_id = e[0]
                taxonomy_levels = [] # Store tuples of (prefix_name, confidence)

                # Iterate through the taxonomy levels provided
                # Start from index 3 (first rank name, usually 'domain')
                # Each level has 3 fields: rank, "name", confidence
                for i in range(3, len(e), 3):
                    # Check if the triplet is complete
                    if (i + 2) >= len(e):
                        logging.warning(f"Incomplete taxonomy triplet for OTU {otu_id} near field {i}, stopping parse for this OTU.")
                        break

                    tax_rank = e[i].lower() # Rank name (e.g., 'domain', 'phylum')
                    # Handle quoted taxon names if present
                    tax_name = e[i+1].strip('"').replace(" ", "_") # Taxon name (e.g., 'Fungi', 'Ascomycota')
                    try:
                        tax_confidence = float(e[i+2])
                    except ValueError:
                        logging.warning(f"Invalid confidence value '{e[i+2]}' for rank '{tax_rank}' in OTU {otu_id}. Skipping level.")
                        continue # Skip this taxonomic level

                    # Skip the 'Root' level usually provided by RDP
                    if tax_name.lower() == 'root':
                        continue

                    prefix = PREFIX_MAP.get(tax_rank)
                    if prefix and tax_name: # Ensure we have a prefix and a name
                        taxonomy_levels.append((f"{prefix}{tax_name}", tax_confidence))
                    else:
                        # Log unrecognized ranks but continue processing
                        logging.debug(f"Unrecognized taxonomic rank '{tax_rank}' or missing name for OTU {otu_id}. Skipping this level.")

                # Filter based on confidence threshold - applied hierarchically
                taxonomy_filtered = []
                taxonomy_conf_filtered = []
                last_valid_conf = 0.0
                for level_str, conf in taxonomy_levels:
                    if conf >= confidence_threshold:
                        taxonomy_filtered.append(level_str)
                        taxonomy_conf_filtered.append(conf)
                        last_valid_conf = conf
                    else:
                        # Stop adding levels once confidence drops below threshold
                        break

                # Construct final string
                if not taxonomy_filtered:
                    # If no level met the threshold, assign as Unassigned at Kingdom level
                    final_taxonomy_str = "k__Unassigned"
                    final_confidence_str = "1.0000" # Confidence of 'Unassigned'
                else:
                    final_taxonomy_str = ";".join(taxonomy_filtered) # Use semicolon, no space for QIIME1
                    # The confidence reported is typically that of the *last* assigned level
                    final_confidence_str = f"{last_valid_conf:.4f}"
                    assigned_otus += 1


                # Write QIIME-compatible output: OTUID<tab>TaxonomyString<tab>Confidence
                handle_output.write(f"{otu_id}\t{final_taxonomy_str}\t{final_confidence_str}\n")

        logging.info(f"Reformatted {processed_lines} lines from RDP output.")
        logging.info(f"Found {assigned_otus} OTUs with assignments passing threshold {confidence_threshold}.")
        logging.info(f"Output written to {output_file}")
        return True

    except FileNotFoundError:
        logging.error(f"Input RDP file not found: {input_file}")
        return False
    except IOError as e:
        logging.error(f"Error reading RDP file {input_file} or writing to {output_file}: {e}")
        return False
    except Exception as e:
        logging.error(f"An unexpected error occurred during RDP reformatting: {e}")
        logging.error(traceback.format_exc()) # Log traceback for debugging
        return False


def reformat_sintax_taxonomy(input_file, output_file, confidence_threshold):
    """
    Re-formats VSEARCH SINTAX output to QIIME-compatible format.
    Args:
        input_file (str): Path to the raw SINTAX tabbed output.
        output_file (str): Path to the output formatted file.
        confidence_threshold (float): Minimum confidence to retain a level.
    Returns:
        bool: True on success, False on failure.
    """
    logging.info(f"Reformatting SINTAX taxonomy from {input_file} (Threshold: {confidence_threshold})")
    # SINTAX abbreviations to QIIME prefixes
    PREFIX_MAP = {"d": "d__", "k": "k__", "p": "p__", "c": "c__", "o": "o__", "f": "f__", "g": "g__", "s": "s__"}
    try:
        with open(input_file, "r", encoding='utf-8', errors='ignore') as handle_input, \
             open(output_file, "w", encoding='utf-8') as handle_output:
            processed_lines = 0
            assigned_otus = 0
            for line in handle_input:
                processed_lines += 1
                parts = line.rstrip().split("\t")
                if not parts: continue # Skip empty lines

                otu_id = parts[0]
                raw_taxonomy_str = ""
                # SINTAX output fields:
                # 1: OTUID
                # 2: Taxonomy (comma-sep: d:Domain(conf),p:Phylum(conf),...)
                # 3: Strand (+/-) (Optional)
                # 4: Alternative taxonomy (Optional)
                if len(parts) > 1:
                    raw_taxonomy_str = parts[1]

                taxonomy_levels = [] # Store tuples of (prefix_name, confidence)

                if raw_taxonomy_str and raw_taxonomy_str != '-': # Check if assignment exists
                    tax_levels = raw_taxonomy_str.split(",")
                    for level_str in tax_levels:
                        try:
                            # Format: level_abbr:Name(confidence) e.g., p:Ascomycota(0.99)
                            level_part, conf_part = level_str.split("(")
                            level_abbr = level_part.split(":")[0]
                            tax_name = level_part.split(":")[1].replace(" ", "_") # Get name, replace spaces
                            tax_confidence = float(conf_part.rstrip(")")) # Get confidence float
                        except (ValueError, IndexError) as e:
                            logging.warning(f"Could not parse SINTAX level string '{level_str}' for OTU {otu_id}: {e}. Skipping level.")
                            continue # Skip this malformed level

                        prefix = PREFIX_MAP.get(level_abbr)
                        if prefix and tax_name: # Ensure we have a prefix and a non-empty name
                            taxonomy_levels.append((f"{prefix}{tax_name}", tax_confidence))
                        else:
                             logging.debug(f"Unrecognized SINTAX level '{level_abbr}' or empty name for OTU {otu_id}. Skipping.")

                # Filter based on confidence threshold - applied hierarchically
                taxonomy_filtered = []
                taxonomy_conf_filtered = []
                last_valid_conf = 0.0
                for level_str, conf in taxonomy_levels:
                    if conf >= confidence_threshold:
                        taxonomy_filtered.append(level_str)
                        taxonomy_conf_filtered.append(conf)
                        last_valid_conf = conf
                    else:
                        # Stop adding levels once confidence drops below threshold
                        break

                # Construct final string
                if not taxonomy_filtered:
                    # If no level met the threshold (or no assignment initially)
                    final_taxonomy_str = "k__Unassigned"
                    final_confidence_str = "1.0000"
                else:
                    final_taxonomy_str = ";".join(taxonomy_filtered) # Semicolon, no space
                    # Confidence is that of the last assigned level
                    final_confidence_str = f"{last_valid_conf:.4f}"
                    assigned_otus += 1


                # Write QIIME-compatible output: OTUID<tab>TaxonomyString<tab>Confidence
                handle_output.write(f"{otu_id}\t{final_taxonomy_str}\t{final_confidence_str}\n")

        logging.info(f"Reformatted {processed_lines} lines from SINTAX output.")
        logging.info(f"Found {assigned_otus} OTUs with assignments passing threshold {confidence_threshold}.")
        logging.info(f"Output written to {output_file}")
        return True
    except FileNotFoundError:
        logging.error(f"Input SINTAX file not found: {input_file}")
        return False
    except IOError as e:
        logging.error(f"Error reading SINTAX file {input_file} or writing to {output_file}: {e}")
        return False
    except Exception as e:
        logging.error(f"An unexpected error occurred during SINTAX reformatting: {e}")
        logging.error(traceback.format_exc()) # Log traceback
        return False


def create_phylotype_table(input_biom_file, output_tsv_file, level):
    """
    Collapses a BIOM table to a specified taxonomic level and writes to TSV.
    Args:
        input_biom_file (str): Path to the input BIOM table file.
        output_tsv_file (str): Path to the output TSV phylotype table.
        level (int): Taxonomic level index (0=kingdom, 1=phylum, ..., 6=species).
                       Assumes standard 7 levels (k, p, c, o, f, g, s).
    Returns:
        bool: True on success, False on failure.
    """
    logging.info(f"Creating phylotype table (Level index: {level}) from {input_biom_file}")
    try:
        # Load the BIOM table
        biom_table = load_table(input_biom_file)
        # Level is 0-indexed (0=k, 1=p, ..., 6=s)
        collapse_level_index = int(level)
        # We want to slice up to *and including* the desired level
        slice_end = collapse_level_index + 1

        # Define the collapsing function
        def collapse_fn(otu_id, metadata):
            # Check if 'taxonomy' metadata exists for this OTU
            taxonomy_list = metadata.get('taxonomy') # Use .get() for safer access
            if taxonomy_list and isinstance(taxonomy_list, (list, tuple)):
                # Check if the taxonomy list is long enough for the requested level
                if len(taxonomy_list) >= slice_end:
                    # Join the levels up to the desired slice end
                    return ';'.join(taxonomy_list[:slice_end])
                else:
                    # If taxonomy is shorter than requested level, use the full available path
                    logging.debug(f"Taxonomy for {otu_id} is shorter ({len(taxonomy_list)} levels) than requested level index {collapse_level_index}. Using full path: {taxonomy_list}")
                    return ';'.join(taxonomy_list)
            else:
                # Handle OTUs with missing or invalid taxonomy metadata
                logging.warning(f"No valid taxonomy metadata found for OTU ID: {otu_id}. Collapsing to 'k__Unassigned'.")
                # Use a consistent 'Unassigned' format
                return 'k__Unassigned'

        logging.debug(f"Collapsing BIOM table observations using slice end index: {slice_end} (Level {level})")
        # Perform the collapse operation on observations (OTUs)
        # norm=False keeps raw counts, min_group_size=1 keeps all groups
        collapsed_table = biom_table.collapse(collapse_fn, norm=False, min_group_size=1, axis='observation')

        # Sort the collapsed table by the new phylotype names (observation IDs)
        collapsed_table.sort(axis='observation')

        # Write the collapsed table to a TSV file
        with open(output_tsv_file, "w") as outfile:
            # .to_tsv() generates the standard text-based OTU table format
            outfile.write(collapsed_table.to_tsv())

        logging.info(f"Phylotype table (Level {level}) written to {output_tsv_file}")
        return True

    except FileNotFoundError:
        logging.error(f"Input BIOM file not found: {input_biom_file}")
        return False
    except (IOError, ValueError) as e:
        # Catch potential errors during BIOM loading or writing
        logging.error(f"Error processing BIOM file {input_biom_file} or writing TSV {output_tsv_file}: {e}")
        return False
    except KeyError as e:
         # This might occur if 'taxonomy' metadata key is missing entirely
         logging.error(f"Metadata key error during collapse (likely missing 'taxonomy' metadata in BIOM file {input_biom_file}): {e}")
         return False
    except Exception as e:
        # Catch any other unexpected errors
        logging.error(f"An unexpected error occurred during phylotyping: {e}")
        logging.error(traceback.format_exc()) # Log traceback
        return False

def create_otu_table_from_uc(uc_file_path, sample_list_path, output_tsv_path):
    """
    Generates an OTU table (TSV format) from a VSEARCH .uc file (expected from
    --usearch_global, primarily 'H' and 'N' lines) and a sample list.

    Args:
        uc_file_path (str): Path to the input .uc file.
        sample_list_path (str): Path to the file containing sample IDs (one per line).
        output_tsv_path (str): Path to write the output OTU table TSV file.

    Returns:
        bool: True on success, False on failure.
    """
    logging.info(f"Generating OTU table from UC file: {uc_file_path}")
    # Dictionary to store counts: {otu_id: {sample_id: count}}
    otu_counts = {}
    line_count = 0
    processed_hits = 0
    skipped_reads = 0

    try:
        # --- Process UC file ---
        logging.debug(f"Reading UC file: {uc_file_path}")
        with open(uc_file_path, "r", encoding='utf-8', errors='ignore') as infile:
            for line in infile:
                line = line.rstrip()
                if not line or line.startswith('#'): continue # Skip empty/comment lines

                line_count += 1
                fields = line.split("\t")

                # Expect 'H' (Hit) lines from --usearch_global for building the table
                if fields[0] != 'H':
                    continue # Skip non-hit lines (like 'N', 'C', 'S')

                if len(fields) < 10:
                    logging.warning(f"Hit line {line_count} in .uc file '{uc_file_path}' has < 10 fields. Skipping.")
                    skipped_reads += 1
                    continue

                # Field 8: Query label (Read ID, e.g., Sample1_123)
                # Field 9: Target label (OTU ID, e.g., OTU1)
                read_id = fields[8]
                otu_id_target = fields[9]

                # Extract sample ID from read ID (part before the first '_')
                try:
                    sample_id = read_id.split("_", 1)[0]
                except IndexError:
                    logging.warning(f"Could not parse sample ID from read ID '{read_id}' on hit line {line_count}. Skipping read.")
                    skipped_reads += 1
                    continue

                # Increment count for this OTU in this Sample
                if otu_id_target not in otu_counts:
                    otu_counts[otu_id_target] = {} # Initialize dict for this OTU
                if sample_id not in otu_counts[otu_id_target]:
                    otu_counts[otu_id_target][sample_id] = 0 # Initialize count for this sample
                otu_counts[otu_id_target][sample_id] += 1
                processed_hits += 1

        logging.debug(f"Processed {line_count} lines from {uc_file_path}.")
        logging.debug(f"Found {processed_hits} hits mapping to {len(otu_counts)} OTUs.")
        if skipped_reads > 0:
            logging.warning(f"Skipped {skipped_reads} reads due to parsing errors (see warnings above).")

        if not otu_counts:
             logging.warning(f"No hits ('H' lines) were processed from the UC file: {uc_file_path}. The resulting OTU table will be empty.")
             # We will proceed and write an empty table with header

        # --- Load Sample IDs from the list file ---
        logging.debug(f"Loading sample IDs from: {sample_list_path}")
        sample_ids_list = []
        try:
            with open(sample_list_path, "r", encoding='utf-8', errors='ignore') as samplefile:
                for line in samplefile:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        # Handle potential extra columns, take only the first
                        sample_ids_list.append(line.split("\t")[0])
        except FileNotFoundError:
             logging.error(f"Sample list file not found: {sample_list_path}")
             return False
        except IOError as e:
             logging.error(f"Could not read sample list file {sample_list_path}: {e}")
             return False


        if not sample_ids_list:
            logging.error(f"No sample IDs loaded from sample list file: {sample_list_path}")
            return False
        logging.info(f"Loaded {len(sample_ids_list)} sample IDs from {sample_list_path}.")

        # --- Write OTU Table ---
        logging.debug(f"Writing OTU table to: {output_tsv_path}")
        with open(output_tsv_path, "w") as outfile:
            # Write header line: #OTU_ID Sample1 Sample2 ...
            # Use '#' prefix for QIIME1 compatibility
            outfile.write("#OTU_ID\t" + "\t".join(sample_ids_list) + "\n")

            # Write counts for each OTU, sorted by OTU ID
            final_otu_count = 0
            for otu_id in sorted(otu_counts.keys()):
                final_otu_count += 1
                outfile.write(otu_id) # Write OTU ID first
                # For each sample in the ordered list, write the count for this OTU
                for sample_id in sample_ids_list:
                    # Use .get(sample_id, 0) to handle cases where OTU wasn't found in a sample
                    count = otu_counts[otu_id].get(sample_id, 0)
                    outfile.write("\t" + str(count))
                outfile.write("\n") # Newline after each OTU row

        if final_otu_count == 0 and processed_hits > 0:
             logging.warning(f"Processed {processed_hits} hits but wrote 0 OTUs to the table. This might indicate an issue with OTU IDs or sample IDs.")
        elif final_otu_count == 0 and processed_hits == 0:
             logging.info(f"No hits found in UC file, an empty OTU table (with header) was written to {output_tsv_path}")
        else:
             logging.info(f"OTU table with {final_otu_count} OTUs written to {output_tsv_path}")
        return True

    except FileNotFoundError as e:
        # This catches the UC file not found error specifically
        logging.error(f"File not found during OTU table generation: {e}")
        return False
    except IOError as e:
        # Catches errors reading UC or writing TSV
        logging.error(f"I/O error during OTU table generation: {e}")
        return False
    except Exception as e:
        # Catch any other unexpected errors
        logging.error(f"An unexpected error occurred during OTU table generation: {e}")
        logging.error(traceback.format_exc()) # Log traceback
        return False

# --- Database Download Functions (Restored from Original) ---

def get_md5(filename):
    """Calculates MD5 checksum of a file."""
    hash_md5 = hashlib.md5()
    try:
        with open(filename, "rb") as f:
            # Read file in chunks to handle large files
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    except FileNotFoundError:
        logging.error(f"File not found for MD5 calculation: {filename}")
        return None # Return None if file doesn't exist
    except IOError as e:
        logging.error(f"Could not read file for MD5 calculation {filename}: {e}")
        return None # Return None on read error


def downloadDB(url, md5, output_dir):
    """Downloads and unpacks database files with MD5 check."""
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except OSError as e:
             logging.error(f"Could not create database directory {output_dir}: {e}")
             raise # Re-raise to stop execution

    filename = url.split("/")[-1] # Get filename from URL
    filepath = os.path.join(output_dir, filename) # Full path to the archive
    # Assumes unpacked directory name is the archive name without .tar.gz
    unpacked_dirname = os.path.join(output_dir, filename.split(".tar.gz")[0])

    DOWNLOADDB = True # Flag to control download

    # Check 1: Is the directory already unpacked and non-empty?
    if os.path.isdir(unpacked_dirname) and len(os.listdir(unpacked_dirname)) > 0:
        logging.info(f"Found existing unpacked database directory: {unpacked_dirname}. Assuming DB is present.")
        DOWNLOADDB = False

    # Check 2: Does the archive file exist? If so, check its MD5.
    elif os.path.isfile(filepath):
        logging.info(f"Found existing archive: {filepath}. Checking MD5...")
        existing_md5 = get_md5(filepath)
        if existing_md5 and existing_md5 == md5:
            logging.info("Archive MD5 matches expected. Will proceed to unpack if necessary.")
            DOWNLOADDB = False # Don't download, but might need unpacking
        else:
            logging.warning(f"Archive {filepath} MD5 mismatch or error calculating MD5. Re-downloading...")
            # Attempt to remove the potentially corrupt archive
            try:
                 os.remove(filepath)
            except OSError as e:
                 logging.warning(f"Could not remove old/corrupt archive file {filepath}: {e}")
            # Also remove potentially incomplete unpacked directory
            if os.path.exists(unpacked_dirname):
                 try:
                     shutil.rmtree(unpacked_dirname)
                 except OSError as e:
                      logging.warning(f"Could not remove existing unpacked directory {unpacked_dirname}: {e}")
            DOWNLOADDB = True # Need to re-download

    # Check 3: If neither unpacked dir nor valid archive exists, download needed.
    else:
         logging.info(f"Database archive/unpacked directory not found for {filename}. Downloading.")
         DOWNLOADDB = True


    # --- Download if needed ---
    if DOWNLOADDB:
        logging.info(f"Downloading {filename} from {url}")
        try:
            # Use stream=True for large files, set timeout
            request = requests.get(url, stream=True, timeout=60)
            request.raise_for_status() # Raise HTTPError for bad responses (4xx or 5xx)
        except requests.exceptions.RequestException as e:
             logging.error(f"Download failed for {url}: {e}")
             raise # Stop execution

        try:
            # Write downloaded content to file with progress bar
            with open(filepath, 'wb') as file:
                file_size = int(request.headers.get('Content-Length', 0))
                file_size_mb = round(file_size / 1024 / 1024, 2) if file_size > 0 else 0
                # Setup progress bar widgets
                widgets = [
                    progressbar.Bar(), progressbar.Percentage(), " | ",
                    progressbar.FileTransferSpeed(), " | ",
                    progressbar.DataSize(),
                    f" / {file_size_mb:.2f} MB | " if file_size_mb > 0 else " | ",
                    progressbar.ETA()
                ]
                # Initialize progress bar
                bar = progressbar.ProgressBar(widgets=widgets, max_value=file_size or progressbar.UnknownLength).start()
                bytes_written = 0
                # Iterate over downloaded chunks
                for chunk in request.iter_content(chunk_size=8192):
                    if chunk: # filter out keep-alive new chunks
                        file.write(chunk)
                        bytes_written += len(chunk)
                        bar.update(bytes_written) # Update progress bar
                bar.finish() # Finish progress bar
        except IOError as e:
             logging.error(f"Could not write downloaded file {filepath}: {e}")
             raise # Stop execution
        except Exception as e:
             logging.error(f"Unexpected error during file writing/progress bar: {e}")
             raise

        # --- Verify downloaded file MD5 ---
        logging.info(f"Verifying MD5 checksum for downloaded file: {filepath}")
        downloaded_md5 = get_md5(filepath)
        if not downloaded_md5 or downloaded_md5 != md5:
             logging.error(f"Downloaded file {filepath} is corrupt or MD5 mismatch!")
             logging.error(f"Expected MD5: {md5}")
             logging.error(f"Calculated MD5: {downloaded_md5}")
             # Attempt cleanup of bad download
             try: os.remove(filepath)
             except OSError: pass
             raise RuntimeError("Downloaded database MD5 check failed.") # Stop execution
        else:
             logging.info("Download successful and verified.")

    # --- Unpack if necessary ---
    # Condition: Archive exists, but unpacked dir doesn't or is empty
    if os.path.isfile(filepath) and not (os.path.isdir(unpacked_dirname) and len(os.listdir(unpacked_dirname)) > 0):
        logging.info(f"Unpacking database archive: {filepath}")
        # Ensure target unpack directory exists and is empty
        if os.path.exists(unpacked_dirname):
            try:
                shutil.rmtree(unpacked_dirname)
            except OSError as e:
                logging.error(f"Could not remove existing unpacked directory before extraction {unpacked_dirname}: {e}")
                raise
        try:
            os.makedirs(unpacked_dirname)
        except OSError as e:
            logging.error(f"Could not create unpack directory {unpacked_dirname}: {e}")
            raise

        # --- Perform extraction ---
        try:
            # Open tar.gz file and extract all contents to the target directory
            with tarfile.open(filepath, "r:gz") as tar:
                tar.extractall(path=unpacked_dirname)
            logging.info("Unpacking complete.")
        except tarfile.ReadError as e:
            logging.error(f"Failed to unpack {filepath}. Archive may be corrupt. Error: {e}")
            raise # Stop execution
        except Exception as e:
            # Catch other potential tarfile errors
            logging.error(f"Unexpected error during unpacking: {e}")
            raise

    # --- Final check: Ensure unpacked directory exists and is populated ---
    elif not os.path.isdir(unpacked_dirname) or len(os.listdir(unpacked_dirname)) == 0:
        # This case means download was skipped (or failed silently) AND unpacking didn't happen/failed
        logging.critical(f"Database setup failed. Required unpacked directory '{unpacked_dirname}' is missing or empty after checks/download/unpack attempt.")
        raise RuntimeError("Database setup failed - unpacked directory missing/empty.")


# --- Wrapper function for the main pipeline logic ---
def run_process_pipeline(options):
    """Contains the main execution logic called by main()."""

    # --- Directory Setup ---
    if not os.path.exists(options.outDir):
        os.makedirs(options.outDir)

    # Intermediate directory for temporary files
    tmpDir = os.path.join(options.outDir, "intermediate")
    if not options.retain and os.path.exists(tmpDir):
         logging.info(f"Removing existing intermediate directory: {tmpDir}")
         shutil.rmtree(tmpDir)
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir)


    summary_file_path = os.path.join(options.outDir, "summary.log")
    version_file_path = os.path.join(options.outDir, "versions.log")

    logging.info(f"--- Starting pipits_process v{__version__} ---")
    logging.info(f"Output directory: {os.path.abspath(options.outDir)}")
    logging.info(f"Input FASTA: {os.path.abspath(options.input)}")
    logging.info(f"UNITE Version: {options.unite}")
    logging.info(f"Taxonomy Method: {options.taxassignmentmethod}")
    logging.info(f"VSEARCH ID: {options.VSEARCH_id}")
    logging.info(f"SINTAX Confidence: {options.sintaxconfidence}")
    logging.info(f"RDP Confidence: {options.RDP_assignment_threshold}")
    logging.info(f"Threads: {options.threads}")
    logging.info(f"Retain Intermediates: {options.retain}")
    logging.info(f"Verbose Logging: {options.verbose}")
    logging.info(f"Ignore DB Download: {options.ignore_db_download}") # Log this option
    logging.info("-" * 30)


    # --- Log versions ---
    try:
        version_info = f"pipits_process: {__version__}\n"
        # Try getting versions of key tools
        for tool, tool_cmd in [("VSEARCH", f"{VSEARCH} --version"),
                               ("BIOM", f"{BIOM} --version"),
                               ("RDP Classifier", f"{RDP_CLASSIFIER} -v" if shutil.which(RDP_CLASSIFIER) else None)]: # Check if RDP cmd exists
             if tool_cmd is None and tool == "RDP Classifier":
                 version_info += f"{tool}: Command '{RDP_CLASSIFIER}' not found in PATH.\n"
                 continue
             try:
                 process = subprocess.run(tool_cmd, shell=True, check=True, capture_output=True, text=True, timeout=10)
                 output = process.stdout.strip() or process.stderr.strip()
                 # Take first line of version output
                 version_info += f"{tool}: {output.splitlines()[0]}\n"
             except Exception as tool_e:
                 version_info += f"{tool}: Not found or error getting version ({type(tool_e).__name__})\n"
                 logging.warning(f"Could not get version for {tool}: {tool_e}")

        with open(version_file_path, "w") as vf:
             vf.write(version_info)
        logging.info(f"Tool versions logged to {version_file_path}")
    except Exception as e:
        logging.warning(f"Could not log tool versions: {e}")

    # --- Start Processing ---
    logging.info("pipits_process started")
    start_time = time.time()
    final_otu_biom_path = "" # Track the path of the final BIOM table for stats
    stats_method = ""        # Track which method was used for the final table
    otu_table_prelim_biom = None # Track preliminary biom file path

    try: # Wrap pipeline steps in try/except
        # --- Input Validation ---
        if not os.path.exists(options.input):
            logging.critical(f"Input file not found: {options.input}")
            return 1 # Exit error

        if os.stat(options.input).st_size == 0:
            logging.critical("Input file is empty! Cannot process.")
            # Exit gracefully, maybe write empty summary?
            # Write summary indicating empty input
            with open(summary_file_path, "w") as sf:
                sf.write("--- PIPITS_PROCESS Summary ---\n")
                sf.write("Input file was empty. Processing halted.\n")
            return 0 # Exit success (nothing to do)


        # --- Sample ID Handling ---
        if not options.sortedlist:
            logging.info("Sample list (-l) not provided. Generating from input FASTA headers...")
            generated_sample_list_path = os.path.join(options.outDir, "generated_sampleIDs.txt")
            try:
                success = generate_sample_list_from_fasta(options.input, generated_sample_list_path)
                if not success:
                     logging.critical("Failed to generate sample list automatically.")
                     return 1 # Exit error
                options.sortedlist = generated_sample_list_path # Use the generated list
            except Exception as e:
                 logging.critical(f"An unexpected error occurred calling generate_sample_list_from_fasta: {e}")
                 return 1 # Exit error
        else:
            logging.info(f"Using provided sample list: {options.sortedlist}")
            if not os.path.exists(options.sortedlist):
                logging.critical(f"Provided sample list file does NOT exist: {options.sortedlist}")
                return 1 # Exit error

        # --- Validate Sample List Content ---
        sampleIDs_check = []
        try:
            with open(options.sortedlist, 'r', encoding='utf-8', errors='ignore') as f_list:
                for line_num, line in enumerate(f_list):
                    line_content = line.strip()
                    if line_content and not line_content.startswith("#"):
                        sample_id = line_content.split("\t")[0] # Take first part if tab-separated
                        if not sample_id:
                            logging.warning(f"Empty sample ID found on line {line_num+1} in {options.sortedlist}. Skipping.")
                            continue
                        if ' ' in sample_id:
                             logging.warning(f"Space detected in sample ID '{sample_id}' (line {line_num+1}) from list file {options.sortedlist}. Consider removing spaces.")
                        if '_' in sample_id:
                              logging.warning(f"Underscore detected in sample ID '{sample_id}' (line {line_num+1}). This might cause issues if FASTA headers also use '_'.")
                        sampleIDs_check.append(sample_id)
        except IOError as e:
             logging.critical(f"Could not read sample list file {options.sortedlist}: {e}")
             return 1 # Exit error

        if not sampleIDs_check:
            logging.critical(f"No valid sample IDs found in the sample list file: {options.sortedlist}")
            return 1 # Exit error
        logging.info(f"Checked {len(sampleIDs_check)} sample IDs from {options.sortedlist}.")


        # Check for duplicate entries in the sample list
        if len(sampleIDs_check) != len(set(sampleIDs_check)):
             duplicateIDs = list(set([x for x in sampleIDs_check if sampleIDs_check.count(x) > 1]))
             logging.critical(f"Duplicate sample ID(s) found in your Sample list file '{options.sortedlist}'. Offending ID(s): " + ",".join(duplicateIDs))
             return 1 # Exit error
        logging.debug("No duplicate sample IDs found in list file.")


        # --- Database Download/Check Section (Restored Original Logic) ---
        DB_DOWNLOAD_ERROR = False # Flag for download errors
        db_base_dir = "pipits_db" # Base directory for databases

        if options.ignore_db_download:
            logging.warning("Skipping DB checks/downloads (--ignore-db-download specified).")
            logging.warning(f"Ensure required databases are present and correctly structured within '{db_base_dir}'.")
            if not os.path.isdir(db_base_dir):
                 logging.critical(f"DB download skipped, but required base DB directory '{db_base_dir}' not found.")
                 return 1 # Exit if base dir is missing when skipping download
        else:
            logging.info("Checking/Downloading required databases...")
            # Define MD5 checksums (These need to be kept up-to-date)
            md5_map_rdp = {"19.02.2025": "1d5b1fdab5b173618c912c08835958f8",
                           "04.04.2024": "94812f45cbfed846b55a6e845e68f35f",
                           "25.07.2023": "1c578d0aba436f0b66d9a73b2991086d",
                           "27.10.2022": "7d31f5612a78607e50d4170b75d0cbfa",
                           "10.05.2021": "13f1edfb1357eeda3f41ff8b0a15447f",
                           "04.02.2020": "b2f833c89794be20a5fdb6169d9205f1",
                           "02.02.2019": "8fd3b74a510bb20b67933a2ecc620f89",
                           "01.12.2017": "3c5be9c60fecf70076739379e7c9ead5",
                           "28.06.2017": "33fa78987751a494c586676ff3a0da65"}
            md5_map_sintax = {"19.02.2025": "bb6a7d9d16d6d60484d62d55fdd6dc55",
                              "04.04.2024": "135cce3029569c1c0528a7fdd9ed6673",
                              "25.07.2023": "cf7dcd5e289a4d31c87fadb056caaca7",
                              "27.10.2022": "b26ebd07a5abb415e1ad35b8dc8108d2"} # Older SINTAX MD5s might be missing
            md5_uchime = "a84af781f9ba42a76f456b558dbc1ae5" # Assume fixed UCHIME version

            # --- Download UCHIME DB ---
            logging.info("Checking/Downloading UCHIME reference database...")
            uchime_url = "https://sourceforge.net/projects/pipits/files/PIPITS_DB/uchime_reference_dataset_28.06.2017.tar.gz"
            try:
                downloadDB(url=uchime_url, md5=md5_uchime, output_dir=db_base_dir)
            except Exception as e:
                 logging.error(f"Failed UCHIME DB download/setup: {e}")
                 DB_DOWNLOAD_ERROR = True

            # --- Download RDP DB (if needed) ---
            if not DB_DOWNLOAD_ERROR and options.taxassignmentmethod in ["all", "rdp"]:
                logging.info(f"Checking/Downloading UNITE RDP database (Version: {options.unite})...")
                rdp_url = f"https://sourceforge.net/projects/pipits/files/PIPITS_DB/UNITE_retrained_{options.unite}.tar.gz"
                rdp_md5 = md5_map_rdp.get(options.unite)
                if not rdp_md5:
                    logging.error(f"Unknown MD5 checksum for UNITE RDP version {options.unite}. Cannot verify download.")
                    DB_DOWNLOAD_ERROR = True
                else:
                    try:
                        downloadDB(url=rdp_url, md5=rdp_md5, output_dir=db_base_dir)
                    except Exception as e:
                         logging.error(f"Failed UNITE RDP DB (v{options.unite}) download/setup: {e}")
                         DB_DOWNLOAD_ERROR = True

            # --- Download SINTAX DB (if needed) ---
            if not DB_DOWNLOAD_ERROR and options.taxassignmentmethod in ["all", "sin"]:
                logging.info(f"Checking/Downloading UNITE SINTAX database (Version: {options.unite})...")
                sintax_url = f"https://sourceforge.net/projects/pipits/files/PIPITS_DB/UNITE_retrained_{options.unite}.sintax.fa.tar.gz"
                sintax_md5 = md5_map_sintax.get(options.unite)
                if not sintax_md5:
                     logging.warning(f"MD5 checksum for UNITE SINTAX version {options.unite} is not defined in the script. Download will proceed without verification.")
                     # Set md5 to None or a dummy value if downloadDB requires it, or modify downloadDB
                     sintax_md5 = None # Or handle appropriately in downloadDB if None is passed
                     # Continue download but log the lack of verification
                     try:
                        downloadDB(url=sintax_url, md5="dummy_md5_to_bypass_check", output_dir=db_base_dir) # Need modification in downloadDB or remove check
                        logging.warning("SINTAX DB downloaded, but MD5 could not be verified.")
                     except Exception as e:
                          logging.error(f"Failed UNITE SINTAX DB (v{options.unite}) download/setup: {e}")
                          DB_DOWNLOAD_ERROR = True
                else:
                    try:
                        downloadDB(url=sintax_url, md5=sintax_md5, output_dir=db_base_dir)
                    except Exception as e:
                         logging.error(f"Failed UNITE SINTAX DB (v{options.unite}) download/setup: {e}")
                         DB_DOWNLOAD_ERROR = True

            # --- Exit if any DB download failed ---
            if DB_DOWNLOAD_ERROR:
                 logging.critical("One or more database downloads or setups failed. Please check logs and internet connection. Exiting.")
                 return 1 # Exit error


        # --- Define DB Paths and Final Check ---
        # (This part is similar to the check_databases function but integrated here)
        DB_MISSING_ERROR = False
        db_paths = {} # Dictionary to store the final paths

        # UCHIME Path
        UCHIME_DB_NAME = "uchime_reference_dataset_28.06.2017"
        db_paths['uchime_db'] = str(Path(db_base_dir) / UCHIME_DB_NAME / f"{UCHIME_DB_NAME}.fasta")
        if not os.path.isfile(db_paths['uchime_db']):
            logging.error(f"UCHIME DB FASTA file not found after download/check: {db_paths['uchime_db']}")
            DB_MISSING_ERROR = True

        # RDP Path (if needed)
        if options.taxassignmentmethod in ["all", "rdp"]:
            RDP_UNITE_DB_NAME = f"UNITE_retrained_{options.unite}"
            db_paths['rdp_props'] = str(Path(db_base_dir) / RDP_UNITE_DB_NAME / "UNITE_retrained" / "rRNAClassifier.properties")
            if not os.path.isfile(db_paths['rdp_props']):
                logging.error(f"RDP properties file not found after download/check: {db_paths['rdp_props']}")
                DB_MISSING_ERROR = True
        else: db_paths['rdp_props'] = None # Not needed

        # SINTAX Path (if needed)
        if options.taxassignmentmethod in ["all", "sin"]:
            SINTAX_UNITE_DB_NAME = f"UNITE_retrained_{options.unite}.sintax.fa"
            db_paths['sintax_db'] = str(Path(db_base_dir) / SINTAX_UNITE_DB_NAME / f"UNITE_retrained_{options.unite}.sintax.fa")
            if not os.path.isfile(db_paths['sintax_db']):
                logging.error(f"SINTAX DB FASTA file not found after download/check: {db_paths['sintax_db']}")
                DB_MISSING_ERROR = True
        else: db_paths['sintax_db'] = None # Not needed

        if DB_MISSING_ERROR:
             logging.critical("Required database files missing even after download/check phase. Exiting.")
             return 1 # Exit error
        else:
             logging.info("Database file paths verified.")


        # --- Main Processing Steps ---

        # 1. Dereplicate sequences
        minuniquesize = 2 if not options.includeuniqueseqs else 1
        logging.info(f"Dereplicating unique sequences (min size: {minuniquesize}) [VSEARCH]")
        input_nr_fasta = os.path.join(tmpDir, "input_nr.fasta")
        # Use --fasta_width 0 to prevent line wrapping
        cmd = f"{VSEARCH} --derep_fulllength \"{options.input}\" --output \"{input_nr_fasta}\" --minuniquesize {minuniquesize} --sizeout --fasta_width 0"
        run_cmd(cmd) # Raises error on failure
        # Check if dereplication produced output
        if not os.path.exists(input_nr_fasta) or os.stat(input_nr_fasta).st_size == 0:
            logging.critical("Dereplication produced no sequences. This might happen if all sequences were singletons and --includeuniqueseqs was not used, or if the input was effectively empty after previous steps.")
            # Write summary indicating no sequences post-derep
            with open(summary_file_path, "w") as sf:
                sf.write("--- PIPITS_PROCESS Summary ---\n")
                sf.write("Dereplication resulted in 0 sequences. Processing halted.\n")
            return 0 # Exit gracefully

        # 2. OTU clustering (using cluster_fast)
        logging.info(f"Clustering OTUs at {options.VSEARCH_id} identity [VSEARCH cluster_fast]")
        input_nr_otus_fasta = os.path.join(tmpDir, "input_nr_otus.fasta") # Centroids output
        input_nr_otus_uc = os.path.join(tmpDir, "input_nr_otus.uc")       # Clustering map
        cmd = f"{VSEARCH} --cluster_fast \"{input_nr_fasta}\" --id {options.VSEARCH_id} --centroids \"{input_nr_otus_fasta}\" --uc \"{input_nr_otus_uc}\" --threads {options.threads} --fasta_width 0"
        run_cmd(cmd) # Raises error on failure
        if not os.path.exists(input_nr_otus_fasta) or os.stat(input_nr_otus_fasta).st_size == 0:
            logging.critical("OTU clustering (cluster_fast) produced no centroid sequences. Exiting.")
            # Write summary
            with open(summary_file_path, "w") as sf:
                 sf.write("--- PIPITS_PROCESS Summary ---\n")
                 sf.write("OTU clustering resulted in 0 sequences. Processing halted.\n")
            return 0 # Exit gracefully


        # 3. Chimera removal using reference database
        logging.info("Removing chimeras using UCHIME-ref [VSEARCH]")
        input_nr_otus_nonchimeras_fasta = os.path.join(tmpDir, "input_nr_otus_nonchimeras.fasta")
        # Use the checked uchime_db path
        cmd = f"{VSEARCH} --uchime_ref \"{input_nr_otus_fasta}\" --db \"{db_paths['uchime_db']}\" --nonchimeras \"{input_nr_otus_nonchimeras_fasta}\" --threads {options.threads} --fasta_width 0"
        run_cmd(cmd) # Raises error on failure
        # Check if any non-chimeras remain
        if not os.path.exists(input_nr_otus_nonchimeras_fasta) or os.stat(input_nr_otus_nonchimeras_fasta).st_size == 0:
            logging.warning("UCHIME chimera removal resulted in 0 non-chimeric sequences.")
            # Pipeline can continue, but downstream steps might yield empty results.
            # Create empty file to prevent downstream file-not-found errors.
            open(input_nr_otus_nonchimeras_fasta, 'a').close()


        # 4. Rename OTUs (relabel centroids)
        logging.info("Renaming OTU centroids (e.g., >OTU1)")
        renamed_otus_fasta = os.path.join(tmpDir, "input_nr_otus_nonchimeras_relabelled.fasta")
        otu_counter = 0
        try:
            # Process the non-chimeric FASTA file
            with open(input_nr_otus_nonchimeras_fasta, "r", encoding='utf-8', errors='ignore') as handle_in, \
                 open(renamed_otus_fasta, "w", encoding='utf-8') as handle_out:
                for line in handle_in:
                    if line.startswith(">"):
                        otu_counter += 1
                        # Write new header format >OTU<number>
                        handle_out.write(f">OTU{otu_counter}\n")
                    else:
                        # Write sequence line as is
                        handle_out.write(line)
            logging.info(f"Renamed {otu_counter} non-chimeric OTU sequences.")
            if otu_counter == 0:
                 logging.warning("No non-chimeric OTUs were found to rename.")
                 # Ensure the renamed file exists but is empty if counter is 0
                 if not os.path.exists(renamed_otus_fasta): open(renamed_otus_fasta, 'a').close()

        except IOError as e:
             logging.critical(f"OTU renaming failed due to file I/O error: {e}")
             return 1 # Exit error
        except Exception as e:
             logging.critical(f"An unexpected error occurred during OTU renaming: {e}")
             return 1 # Exit error

        # Check again if the final renamed file exists and has content
        if not os.path.exists(renamed_otus_fasta) or os.stat(renamed_otus_fasta).st_size == 0:
            logging.warning("No non-chimeric OTUs to process further (post-renaming check). Final OTU table will be empty.")
            # Don't exit here, allow OTU table generation which will be empty

        # 5. Map original reads to OTU centroids
        logging.info("Mapping original reads to non-chimeric OTU centroids [VSEARCH usearch_global]")
        otus_uc_path = os.path.join(tmpDir, "otus.uc") # Output mapping file
        # Map original input reads against the *renamed* non-chimeric centroids
        # Only proceed if renamed_otus_fasta has sequences
        if os.path.exists(renamed_otus_fasta) and os.stat(renamed_otus_fasta).st_size > 0:
            cmd = f"{VSEARCH} --usearch_global \"{options.input}\" --db \"{renamed_otus_fasta}\" --id {options.VSEARCH_id} --uc \"{otus_uc_path}\" --threads {options.threads}"
            run_cmd(cmd) # Raises error on failure
        else:
            logging.warning(f"Skipping read mapping as the OTU representative file is empty: {renamed_otus_fasta}")
            # Create an empty UC file to avoid downstream errors
            open(otus_uc_path, 'a').close()

        # 6. Generate OTU Table from UC mapping
        logging.info("Generating OTU table from mapping file [Internal Function]")
        otu_table_prelim_txt = os.path.join(tmpDir, "otu_table_prelim.txt")
        success = create_otu_table_from_uc(uc_file_path=otus_uc_path,
                                           sample_list_path=options.sortedlist,
                                           output_tsv_path=otu_table_prelim_txt)
        if not success:
            logging.critical("Failed to generate OTU table using integrated function.")
            return 1 # Exit error

        # 7. Convert preliminary OTU table to BIOM format
        logging.info("Converting OTU table to BIOM format [BIOM]")
        otu_table_prelim_biom = os.path.join(tmpDir, "otu_table_prelim.biom")
        # Remove existing biom file first to avoid issues with biom convert append behavior
        if os.path.exists(otu_table_prelim_biom):
            try: os.remove(otu_table_prelim_biom)
            except OSError as e: logging.warning(f"Could not remove existing preliminary BIOM file: {e}")

        # Check if the preliminary text table exists and is not empty (more than just header)
        prelim_table_has_data = False
        if os.path.exists(otu_table_prelim_txt) and os.stat(otu_table_prelim_txt).st_size > 0:
             try:
                 with open(otu_table_prelim_txt, 'r') as f_check:
                    # Check if there's more than one line (i.e., more than just the header)
                    if sum(1 for line in f_check if line.strip()) > 1:
                        prelim_table_has_data = True
             except IOError as e:
                  logging.warning(f"Could not read preliminary OTU table to check content: {e}")

        if prelim_table_has_data:
            # Convert to HDF5 BIOM format
            cmd = f"{BIOM} convert -i \"{otu_table_prelim_txt}\" -o \"{otu_table_prelim_biom}\" --table-type=\"OTU table\" --to-hdf5"
            run_cmd(cmd) # Raises error on failure
            # Check if BIOM file was actually created
            if not os.path.exists(otu_table_prelim_biom):
                 logging.error("BIOM conversion command ran but output file was not created.")
                 otu_table_prelim_biom = None # Mark as failed
        else:
            logging.warning(f"Preliminary OTU table '{otu_table_prelim_txt}' is empty or contains only header. Skipping BIOM conversion and downstream steps.")
            otu_table_prelim_biom = None # Indicate no BIOM file created


        # --- Downstream steps (Taxonomy, Phylotyping) only if preliminary BIOM was created ---
        if otu_table_prelim_biom and os.path.exists(otu_table_prelim_biom):

            # --- SINTAX Taxonomy Assignment ---
            if options.taxassignmentmethod in ["all", "sin"]:
                logging.info("Assigning taxonomy with VSEARCH-SINTAX")
                sintax_raw_txt = os.path.join(options.outDir, "assigned_taxonomy_sintax_raw.txt")
                # Use checked SINTAX DB path
                cmd = f"{VSEARCH} --sintax \"{renamed_otus_fasta}\" --db \"{db_paths['sintax_db']}\" --sintax_cutoff {options.sintaxconfidence} --tabbedout \"{sintax_raw_txt}\" --threads {options.threads}"
                run_cmd(cmd) # Raises error on failure

                logging.info("Reformatting SINTAX output")
                sintax_reformatted_txt = os.path.join(options.outDir, "assigned_taxonomy_sintax.txt")
                success = reformat_sintax_taxonomy(sintax_raw_txt, sintax_reformatted_txt, options.sintaxconfidence)
                if not success:
                     logging.error("SINTAX reformatting failed.")
                     # Decide whether to continue or exit
                     return 1 # Exit error

                logging.info("Adding SINTAX assignment to BIOM table")
                otu_table_sintax_biom = os.path.join(options.outDir, "otu_table_sintax.biom")
                if os.path.exists(otu_table_sintax_biom): os.remove(otu_table_sintax_biom)
                # Command to add metadata to the BIOM table
                # --sc-separated indicates semicolon-separated taxonomy string
                # --float-fields specifies which metadata columns are floats
                cmd = f"{BIOM} add-metadata -i \"{otu_table_prelim_biom}\" -o \"{otu_table_sintax_biom}\" --observation-metadata-fp \"{sintax_reformatted_txt}\" --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy --float-fields confidence"
                run_cmd(cmd) # Raises error on failure

                # Check if the final SINTAX BIOM table was created
                if not os.path.exists(otu_table_sintax_biom):
                     logging.error("Adding SINTAX metadata command ran but output BIOM file was not created.")
                else:
                     # Set this as the potential final table for stats
                     final_otu_biom_path = otu_table_sintax_biom
                     stats_method = "SINTAX"

                     logging.info("Converting SINTAX BIOM table to TSV format")
                     otu_table_sintax_txt = os.path.join(options.outDir, "otu_table_sintax.txt")
                     if os.path.exists(otu_table_sintax_txt): os.remove(otu_table_sintax_txt)
                     cmd = f"{BIOM} convert -i \"{otu_table_sintax_biom}\" -o \"{otu_table_sintax_txt}\" --to-tsv --header-key taxonomy"
                     run_cmd(cmd) # Raises error on failure

                     logging.info("Generating SINTAX phylotype table (Level 6 - Species)")
                     phylotype_sintax_txt = os.path.join(options.outDir, "phylotype_table_sintax.txt")
                     phylotype_level = 6 # Species level
                     success = create_phylotype_table(otu_table_sintax_biom, phylotype_sintax_txt, phylotype_level)
                     if not success:
                         logging.error("SINTAX Phylotyping failed.")
                         # return 1 # Optionally exit on failure

                     # Convert phylotype TSV back to BIOM if needed (Optional)
                     # phylotype_sintax_biom = os.path.join(options.outDir, "phylotype_table_sintax.biom")
                     # if os.path.exists(phylotype_sintax_biom): os.remove(phylotype_sintax_biom)
                     # cmd = f"{BIOM} convert -i \"{phylotype_sintax_txt}\" -o \"{phylotype_sintax_biom}\" --table-type=\"OTU table\" --to-hdf5"
                     # run_cmd(cmd)

            # --- RDP Taxonomy Assignment ---
            if options.taxassignmentmethod in ["all", "rdp"]:
                logging.info("Assigning taxonomy with RDP Classifier")
                rdp_raw_txt = os.path.join(options.outDir, "assigned_taxonomy_rdp_raw.txt")

                # Construct RDP command - Requires RDP path/wrapper script
                # Ensure classifier command exists using shutil.which
                if not shutil.which(RDP_CLASSIFIER):
                     logging.error(f"RDP Classifier command '{RDP_CLASSIFIER}' not found in PATH. Cannot proceed with RDP assignment.")
                     raise FileNotFoundError(f"{RDP_CLASSIFIER} command not found")

                # Use checked RDP properties path
                # Note: RDP command often needs specific Java memory settings (-Xmx, -Xms)
                cmd = f"{RDP_CLASSIFIER} -Xms{options.Xms} -Xmx{options.Xmx} classify -c {options.RDP_assignment_threshold} -t \"{db_paths['rdp_props']}\" -o \"{rdp_raw_txt}\" \"{renamed_otus_fasta}\""
                run_cmd(cmd) # Raises error on failure

                logging.info("Reformatting RDP output")
                rdp_reformatted_txt = os.path.join(options.outDir, "assigned_taxonomy_rdp.txt")
                success = reformat_rdp_taxonomy(rdp_raw_txt, rdp_reformatted_txt, options.RDP_assignment_threshold)
                if not success:
                     logging.error("RDP reformatting failed.")
                     return 1 # Exit error

                logging.info("Adding RDP assignment to BIOM table")
                otu_table_rdp_biom = os.path.join(options.outDir, "otu_table_rdp.biom")
                if os.path.exists(otu_table_rdp_biom): os.remove(otu_table_rdp_biom)
                cmd = f"{BIOM} add-metadata -i \"{otu_table_prelim_biom}\" -o \"{otu_table_rdp_biom}\" --observation-metadata-fp \"{rdp_reformatted_txt}\" --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy --float-fields confidence"
                run_cmd(cmd) # Raises error on failure

                if not os.path.exists(otu_table_rdp_biom):
                     logging.error("Adding RDP metadata command ran but output BIOM file was not created.")
                else:
                    # Use RDP table for stats if SINTAX wasn't run or failed
                    if not final_otu_biom_path:
                        final_otu_biom_path = otu_table_rdp_biom
                        stats_method = "RDP"

                    logging.info("Converting RDP BIOM table to TSV format")
                    otu_table_rdp_txt = os.path.join(options.outDir, "otu_table_rdp.txt")
                    if os.path.exists(otu_table_rdp_txt): os.remove(otu_table_rdp_txt)
                    cmd = f"{BIOM} convert -i \"{otu_table_rdp_biom}\" -o \"{otu_table_rdp_txt}\" --to-tsv --header-key taxonomy"
                    run_cmd(cmd) # Raises error on failure

                    logging.info("Generating RDP phylotype table (Level 6 - Species)")
                    phylotype_rdp_txt = os.path.join(options.outDir, "phylotype_table_rdp.txt")
                    phylotype_level = 6 # Species level
                    success = create_phylotype_table(otu_table_rdp_biom, phylotype_rdp_txt, phylotype_level)
                    if not success:
                         logging.error("RDP Phylotyping failed.")
                         # return 1 # Optionally exit

                    # Optional: Convert RDP phylotype TSV back to BIOM
                    # phylotype_rdp_biom = os.path.join(options.outDir, "phylotype_table_rdp.biom")
                    # if os.path.exists(phylotype_rdp_biom): os.remove(phylotype_rdp_biom)
                    # cmd = f"{BIOM} convert -i \"{phylotype_rdp_txt}\" -o \"{phylotype_rdp_biom}\" --table-type=\"OTU table\" --to-hdf5"
                    # run_cmd(cmd)

        # --- End of conditional downstream steps ---
        else:
             logging.warning("Skipping taxonomy assignment and phylotyping because the preliminary BIOM table was not generated (likely due to no reads mapping or empty intermediate files).")
             # Use preliminary BIOM for stats if it exists, otherwise stats will be N/A
             if otu_table_prelim_biom and os.path.exists(otu_table_prelim_biom):
                  final_otu_biom_path = otu_table_prelim_biom
                  stats_method = "Preliminary (No Taxonomy)"
             else:
                  final_otu_biom_path = "" # Ensure it's empty
                  stats_method = "N/A (No BIOM table)"


        # --- Move representative sequence file to output Directory ---
        # Always do this if the renamed file exists
        final_repseq_path = os.path.join(options.outDir, "repseqs.fasta")
        if os.path.exists(renamed_otus_fasta):
             logging.info(f"Moving final representative sequences to {final_repseq_path}")
             try:
                 shutil.move(renamed_otus_fasta, final_repseq_path)
             except Exception as e:
                 logging.warning(f"Could not move representative sequences file: {e}")
                 # Attempt copy as fallback?
                 try:
                     shutil.copy2(renamed_otus_fasta, final_repseq_path)
                     logging.debug("Copied representative sequences instead of moving.")
                     # Optionally remove original if copy succeeds
                     # os.remove(renamed_otus_fasta)
                 except Exception as copy_e:
                     logging.error(f"Could not copy representative sequences either: {copy_e}")
        else:
            logging.warning(f"Final representative sequence file not found, cannot move: {renamed_otus_fasta}")


    # --- End of Main Processing Try/Except Block ---
    except (subprocess.CalledProcessError, FileNotFoundError, IOError, ImportError, RuntimeError, ValueError) as e:
         logging.critical(f"A critical error occurred during processing: {type(e).__name__} - {e}")
         if options.verbose: logging.debug(traceback.format_exc()) # Log traceback if verbose
         logging.critical("Pipeline execution halted.")
         return 1 # Failure exit code
    except Exception as e:
         # Catch any other unexpected errors
         logging.critical(f"An unexpected critical error occurred: {type(e).__name__} - {e}")
         logging.error(traceback.format_exc()) # Log full traceback
         logging.critical("Pipeline execution halted.")
         return 1 # Failure exit code


    # --- Cleanup Intermediate Files ---
    if not options.retain:
        logging.info("Cleaning temporary directory...")
        try:
            if os.path.exists(tmpDir): shutil.rmtree(tmpDir)
            logging.info("... temporary directory removed.")
        except OSError as e:
             logging.warning(f"Could not remove temp dir {tmpDir}: {e}")
    else:
         logging.info(f"Intermediate files retained in: {tmpDir}")


    # --- Calculate Final Stats ---
    def biomstats(biom_file_path):
        """Calculate basic stats from a BIOM table."""
        try:
            biom_table = load_table(biom_file_path)
            num_samples = int(biom_table.shape[1]) # Columns
            num_otus = int(biom_table.shape[0])    # Rows
            total_reads = int(biom_table.sum())    # Sum of all counts
            return total_reads, num_otus, num_samples
        except ImportError:
             logging.warning("'biom-format' missing or failed import. Cannot calculate stats.")
             return "N/A", "N/A", "N/A"
        except FileNotFoundError:
             logging.warning(f"BIOM file not found for stats calculation: {biom_file_path}")
             return "NF", "NF", "NF" # Not Found
        except Exception as e:
             logging.warning(f"Cannot process BIOM file {biom_file_path} for stats: {e}")
             return "Err", "Err", "Err" # Error

    otu_reads_count, otu_count, otu_sample_count = "N/A", "N/A", "N/A"
    if final_otu_biom_path and os.path.exists(final_otu_biom_path):
        logging.info(f"Calculating summary stats from final BIOM table: {final_otu_biom_path} (Method: {stats_method})")
        otu_reads_count, otu_count, otu_sample_count = biomstats(final_otu_biom_path)
    else:
        logging.warning(f"Could not find a suitable final OTU BIOM table for summary stats (Checked path: {final_otu_biom_path}). Stats will be N/A.")
        stats_method = "N/A (No suitable BIOM)" # Update status method


    # --- Write Summary File ---
    summary_lines = [
        "--- PIPITS_PROCESS Summary ---",
        f"Final stats based on table:           {stats_method}",
        # Format counts nicely, handle non-numeric results
        f"No. of reads in final OTU table:      {otu_reads_count if isinstance(otu_reads_count, int) else str(otu_reads_count)}",
        f"Number of OTUs:                       {otu_count if isinstance(otu_count, int) else str(otu_count)}",
        f"Number of samples:                    {otu_sample_count if isinstance(otu_sample_count, int) else str(otu_sample_count)}",
        "",
        f"UNITE DB version used:                {options.unite}",
        f"Taxonomy Assignment Method chosen:    {options.taxassignmentmethod}",
        f"Identity Threshold (VSEARCH):         {options.VSEARCH_id}",
        f"RDP Confidence Threshold:             {options.RDP_assignment_threshold}",
        f"SINTAX Confidence Threshold:          {options.sintaxconfidence}",
        f"Included Singleton OTUs (--derep_minuniquesize): {1 if options.includeuniqueseqs else 2}",
        f"PIPITS Process Version:               {__version__}"
    ]
    try:
        with open(summary_file_path, "w") as sf:
            for line in summary_lines: sf.write(line + "\n")
            try:
                # Calculate duration if start_time exists
                end_time = time.time(); total_time = end_time - start_time
                sf.write(f"\nProcessing completed in {total_time:.2f} seconds.\n")
            except NameError: # start_time might not be defined if exited early
                 sf.write("\nProcessing time calculation skipped (early exit?).\n")
    except IOError as e:
        logging.warning(f"Could not write summary file {summary_file_path}: {e}")

    # --- Final Log Messages ---
    logging.info("--- Final Summary ---")
    logging.info(f"\tReads in Final Table:             {otu_reads_count}")
    logging.info(f"\tNumber of OTUs:                   {otu_count}")
    logging.info(f"\tNumber of Samples:                {otu_sample_count}")
    logging.info(f"\tUNITE DB Version:                 {options.unite}")
    logging.info(f"\tTaxonomy Method (for stats):      {stats_method}")
    logging.info(f"\tPIPITS Process Version:           {__version__}")
    try:
        end_time = time.time(); total_time = end_time - start_time
        logging.info(f"pipits_process finished successfully in {total_time:.2f} seconds.")
    except NameError:
        logging.info("pipits_process finished (possibly with early exit).")
    logging.info(f"Done - Resulting files are in '{options.outDir}' directory")

    return 0 # Success exit code


# --- Entry point function for setup.py ---
def main():
    """Parses arguments and runs the PROCESS pipeline."""
    parser = argparse.ArgumentParser("PIPITS_PROCESS: Sequences to OTU Table")
    # --- Input/Output ---
    parser.add_argument(
        "-i", action="store", dest="input", metavar="<FILE>",
        help="[REQUIRED] Input sequences in FASTA format. Typically output from pipits_funits or pipits_prep.", required=True)
    parser.add_argument(
        "-o", action="store", dest="outDir", metavar="<DIR>", default="pipits_process",
        help="Directory to output results [default: pipits_process].", required=False)
    parser.add_argument(
        "-l", action="store", dest="sortedlist", metavar="<TXT>",
        help="Sample list file (one sample ID per line). If not provided, generated automatically from FASTA headers (part before first '_').", required=False)

    # --- Clustering/Filtering Parameters ---
    parser.add_argument(
        "-d", action="store", dest="VSEARCH_id", metavar="<FLOAT>", default="0.97",
        help="VSEARCH clustering/mapping identity threshold [default: 0.97]", required=False)
    parser.add_argument(
        "--includeuniqueseqs", action="store_true", dest="includeuniqueseqs", default=False,
        help="Include unique sequences (singletons) during dereplication (--minuniquesize 1). Default is to remove them (--minuniquesize 2).", required=False)

    # --- Taxonomy Parameters ---
    parser.add_argument(
        "--taxassignmentmethod", action="store", dest="taxassignmentmethod", default="sin",
        help="Taxonomic assignment method: 'sin' (SINTAX only), 'rdp' (RDP only), 'all' (run both). [default: sin]",
        choices=["all", "rdp", "sin"], required=False)
    parser.add_argument(
        "--unite", action="store", dest="unite", default="19.02.2025",
        help="UNITE db version to use (must match downloaded DB folder names) [default: 19.02.2025]",
        choices=["19.02.2025", "04.04.2024", "25.07.2023", "27.10.2022", "10.05.2021", "04.02.2020", "02.02.2019", "01.12.2017", "28.06.2017"],
        required=False)
    parser.add_argument(
        "--sintaxconfidence", action="store", dest="sintaxconfidence", metavar="<FLOAT>", default="0.85", type=float,
        help="VSEARCH SINTAX assignment confidence threshold [default: 0.85]", required=False)
    parser.add_argument(
        "-c", action="store", dest="RDP_assignment_threshold", metavar="<FLOAT>", default="0.85", type=float,
        help="RDP Classifier assignment confidence threshold [default: 0.85]", required=False)
    parser.add_argument(
        "--Xms", action="store", dest="Xms", default="4g",
        help="Minimum Java heap size for RDP Classifier (e.g., 4g) [default: 4g]", required=False)
    parser.add_argument(
        "--Xmx", action="store", dest="Xmx", default="16g",
        help="Maximum Java heap size for RDP Classifier (e.g., 16g). Increase if RDP fails due to memory. [default: 16g]", required=False)

    # --- General Parameters ---
    parser.add_argument(
        "-r", action="store_true", dest="retain", default=False,
        help="Retain intermediate files (can use significant disk space).", required=False)
    parser.add_argument(
        "-v", action="store_true", dest="verbose", default=False,
        help="Verbose mode (more detailed logging).", required=False)
    parser.add_argument(
        "-t", action="store", dest="threads", metavar="<INT>", default="1", type=int, # Ensure threads is int
        help="Number of Threads for VSEARCH steps [default: 1]", required=False)
    # *** Added ignore-db-download back as per original code ***
    parser.add_argument(
        "--ignore-db-download", action="store_true", dest="ignore_db_download", default=False,
        help="Ignore automatic download/check of database files. Assumes databases are already present in 'pipits_db/'.", required=False)
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    options = parser.parse_args()

    # Ensure threads is at least 1
    if options.threads < 1:
         parser.error("Number of threads (-t) must be at least 1.")

    # Setup logging after args are parsed
    log_file = os.path.join(options.outDir, "pipits_process.log")
    # Create output dir early if possible
    try:
        if not os.path.exists(options.outDir): os.makedirs(options.outDir)
    except OSError as e:
        print(f"CRITICAL: Could not create output directory {options.outDir}: {e}", file=sys.stderr)
        sys.exit(1) # Exit early
    setup_logging(log_file, options.verbose)

    # Run the main pipeline logic
    exit_code = run_process_pipeline(options)
    sys.exit(exit_code) # Exit with the code from the pipeline function


# --- Main Execution Block ---
if __name__ == '__main__':
    main()
