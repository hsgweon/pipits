#!/usr/bin/env python

import sgtk_SeqIO as SeqIO

############################################################
# Argument Options
# Header needs to have phylum -> species all separated by '_'

import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description = "[] indicates optional input (order unimportant)\n\nDereplicates, extracts fungal ITS, re-orientates where necessary, re-inflates sequences. Output should be ready for pipits_process.",
                                 formatter_class = RawTextHelpFormatter)

parser.add_argument("-1",
                    action = "store",
                    dest = "input_1",
                    metavar = "input_1",
                    help = "[REQUIRED] FASTA 1",
                    required = True)
parser.add_argument("-2",
                    action = "store",
                    dest = "input_2",
                    metavar = "input_2",
                    help = "[REQUIRED] FASTA 2",
                    required = True)
parser.add_argument("-o",
                    action = "store",
                    dest = "output",
                    metavar = "output",
                    help = "[REQUIRED] FASTA",
                    required = False)
options = parser.parse_args()

############################################################

outfile = open(options.output, "w")

record_1 = SeqIO.FastaParser(options.input_1)
record_2 = SeqIO.FastaParser(options.input_2)

#combined_keys = record_1.keys() + record_2.keys()

#print(combined_keys)
exit(0)

for i in sorted(list(set(combined_keys))):
    
    try:
        outfile.write(">" + record_1[i].description + "\n")
    except KeyError:
        outfile.write(">" + record_2[i].description + "\n")

    try:
        seq1 = record_1[i].seq
    except KeyError:
        seq1 = ""
    
    try:
        seq2 = record_2[i].seq
    except KeyError:
        seq2 = ""

    outfile.write(seq1 + seq2 + "\n")
    
exit(0)
