#!/usr/bin/python

############################################################
# Argument Options

import argparse
parser = argparse.ArgumentParser("Creates a sample list file from a fasta sequences.")
parser.add_argument("-i",
                    action = "store", 
                    dest = "infile", 
                    metavar = "infile",
                    help = "[REQUIRED]", 
                    required = True)
parser.add_argument("-o",
                    action = "store",
                    dest = "outfile",
                    metavar = "outfile",
                    help = "[REQUIRED]",
                    required = True)
options = parser.parse_args()

############################################################

infile = open(options.infile, "r")
outfile = open(options.outfile, "w")

uniquesampleids = []
with open(options.infile, "r") as f:
    for l in f:
        if l.startswith(">"):
            sampleid = l[1:].split("_")[0]
            if not sampleid in uniquesampleids:
                uniquesampleids.append(sampleid)

for i in uniquesampleids:
    outfile.write(i + "\n")

outfile.close()
