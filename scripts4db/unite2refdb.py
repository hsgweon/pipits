#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
parser = argparse.ArgumentParser(description="Parses UNITE FASTA file for PIPITS reference database.")
parser.add_argument("-i", "--input",
                    dest="input_fasta",
                    metavar="INPUT_FASTA",
                    help="[REQUIRED] Input FASTA file from UNITE.",
                    required=True)
parser.add_argument("-r", "--refseq",
                    dest="out_refseq",
                    metavar="OUTPUT_REFSEQ",
                    help="[REQUIRED] Output Reference FASTA file.",
                    required=True)
parser.add_argument("-t", "--taxonomy",
                    dest="out_reftax",
                    metavar="OUTPUT_TAXONOMY",
                    help="[REQUIRED] Output Reference TAXONOMY file.",
                    required=True)
options = parser.parse_args()

import sys
from Bio import SeqIO
from unidecode import unidecode

def remove_non_ascii(text):
    """
    Transliterates Unicode text into ASCII using unidecode.
    Handles both str and bytes input.
    """
    if isinstance(text, bytes):
        text = text.decode('utf-8', errors='replace')
    return unidecode(text)

if __name__ == '__main__':

    counter = 0

    try:
        with open(options.input_fasta, "r") as infile, \
             open(options.out_refseq, "w") as out_refseq, \
             open(options.out_reftax, "w") as out_reftax:

            print(f"Processing FASTA file: {options.input_fasta}")
            print(f"Writing reference sequences to: {options.out_refseq}")
            print(f"Writing taxonomy to: {options.out_reftax}")

            for record in SeqIO.parse(infile, "fasta"):
                counter += 1

                try:
                    parts = record.description.split("|")
                    if len(parts) < 5:
                        print(f"Warning: Skipping record {counter}. Unexpected header format: {record.description}", file=sys.stderr)
                        continue

                    SH = parts[2]
                    full_lineage_info = parts[4]

                    header_part_for_lineage = full_lineage_info + "_" + SH
                    cleaned_lineage_info = remove_non_ascii(header_part_for_lineage)

                    fl = cleaned_lineage_info.split(";")

                    if len(fl) < 7:
                         print(f"Warning: Skipping record {counter}. Incomplete lineage ({len(fl)} levels found): {cleaned_lineage_info}", file=sys.stderr)
                         continue

                    try:
                        Kingdom_raw = fl[0].split("__")[1] if "__" in fl[0] else fl[0]
                        Phylum_raw  = fl[1].split("__")[1] if "__" in fl[1] else fl[1]
                        Class_raw   = fl[2].split("__")[1] if "__" in fl[2] else fl[2]
                        Order_raw   = fl[3].split("__")[1] if "__" in fl[3] else fl[3]
                        Family_raw  = fl[4].split("__")[1] if "__" in fl[4] else fl[4]
                        Genus_raw   = fl[5].split("__")[1] if "__" in fl[5] else fl[5]
                        Species_raw = fl[6].split("__")[1] if "__" in fl[6] else fl[6]

                    except IndexError:
                         print(f"Warning: Skipping record {counter}. Error splitting rank from identifier (e.g., missing '__'): {cleaned_lineage_info}", file=sys.stderr)
                         continue

                    Kingdom = Kingdom_raw
                    Phylum  = f"{Kingdom}|{Phylum_raw}"
                    Class   = f"{Phylum}|{Class_raw}"
                    Order   = f"{Class}|{Order_raw}"
                    Family  = f"{Order}|{Family_raw}"
                    Genus   = f"{Family}|{Genus_raw}"
                    Species = f"{Genus}|{Species_raw}"

                    reformatted_lineage = f"Root;{Kingdom};{Phylum};{Class};{Order};{Family};{Genus};{Species}"

                    out_refseq.write(f">{counter}\t{reformatted_lineage}\n")
                    out_refseq.write(f"{str(record.seq)}\n")

                    out_reftax.write(f"{counter}\t{reformatted_lineage}\n")

                except Exception as e:
                    print(f"Error processing record {counter} ({record.id}): {e}", file=sys.stderr)
                    print(f"Problematic Header: {record.description}", file=sys.stderr)
                    continue

            print(f"\nFinished processing. Total sequences processed: {counter}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {options.input_fasta}", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error opening or writing file: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)
