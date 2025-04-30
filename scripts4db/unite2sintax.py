#!/usr/bin/env python3

import argparse
import sys
from Bio import SeqIO

def parse_unite_header(header):
    try:
        parts = header.split("|")
        if len(parts) < 5:
            print(f"Warning: Skipping record, unexpected header format: {header}", file=sys.stderr)
            return None

        tax_string = parts[4]
        sh_identifier = parts[2]

        # Use dict comprehension and .get() for concise rank extraction
        found_ranks = {item.split('__')[0]: item.split('__')[1] for item in tax_string.split(';') if '__' in item}

        kingdom = found_ranks.get('k', 'Unknown')
        phylum = found_ranks.get('p', 'Unknown')
        tax_class = found_ranks.get('c', 'Unknown')
        order = found_ranks.get('o', 'Unknown')
        family = found_ranks.get('f', 'Unknown')
        genus = found_ranks.get('g', 'Unknown')
        species = found_ranks.get('s', 'Unknown')

        return kingdom, phylum, tax_class, order, family, genus, species, sh_identifier

    except IndexError:
        print(f"Warning: Skipping record due to IndexError in header: {header}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Warning: Skipping record due to error ({type(e).__name__}) in header: {header}", file=sys.stderr)
        return None

def process_fasta(input_handle, output_handle):
    counter = 1
    for record in SeqIO.parse(input_handle, "fasta"):
        parsed_data = parse_unite_header(record.description.strip())

        if parsed_data:
            kingdom, phylum, tax_class, order, family, genus, species, sh_identifier = parsed_data
            # Construct the new header using an f-string
            new_header = f">{counter};tax=k:{kingdom},p:{phylum},c:{tax_class},o:{order},f:{family},g:{genus},s:{species}_{sh_identifier}"
            output_handle.write(new_header + "\n")
            output_handle.write(str(record.seq) + "\n")
            counter += 1

def main():
    parser = argparse.ArgumentParser(description="Parses UNITE FASTA file headers for SINTAX.")
    parser.add_argument("-i", "--input", dest="input_fasta", metavar="INPUT_FASTA",
                        help="[REQUIRED] Input UNITE FASTA file.", required=True)
    parser.add_argument("-o", "--output", dest="output_fasta", metavar="OUTPUT_FASTA",
                        help="[REQUIRED] Output reformatted FASTA file.", required=True)
    options = parser.parse_args()

    try:
        # Use 'with' for safe file handling
        with open(options.input_fasta, "r") as in_fas, \
             open(options.output_fasta, "w") as out_refseq:
            process_fasta(in_fas, out_refseq)
        print(f"Processing complete. Output: {options.output_fasta}")

    except FileNotFoundError as e:
        print(f"Error: Input file not found - {e}", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error: File I/O error - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
