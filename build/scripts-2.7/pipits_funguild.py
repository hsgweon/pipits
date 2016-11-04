#!/usr/bin/python

############################################################
# Argument Options
# Header needs to have phylum -> species all separated by '_'

import argparse
parser = argparse.ArgumentParser("Converts OTU tables into FUNGuild formatted OTU table.")
parser.add_argument("-i",
                    action = "store", 
                    dest = "infile", 
                    metavar = "infile",
                    help = "[REQUIRED] Input OTU table generated from pipits_process.", 
                    required = True)
parser.add_argument("-o",
                    action = "store",
                    dest = "outfile",
                    metavar = "outfile",
                    help = "[REQUIRED] Output FUNGuild formatted OTU table.",
                    required = True)
options = parser.parse_args()

############################################################

infile = open(options.infile, "r")
outfile = open(options.outfile, "w")

for line in infile:
    if line.startswith("# Constructed from biom file"):
        continue

    if line.startswith("#OTU"):
        outfile.write("OTU ID\t" + "\t".join(line.split("\t")[1:]))
        continue

    taxonomy = line.rstrip().split("\t")[-1]
    taxonomy = taxonomy.replace("; ", ";") # Also remove space!

    outfile.write("\t".join(line.split("\t")[:-1]) + "\t" + taxonomy + "\n")

infile.close()
outfile.close()

