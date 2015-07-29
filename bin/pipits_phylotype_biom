#!/usr/bin/env python

####################
# Argument Options #
####################

import argparse
parser = argparse.ArgumentParser(description = "OTU BIOM to PHYLO TXT")
parser.add_argument("-i",
                    action = "store",
                    dest = "biominputfile",
                    metavar = "biomfile",
                    help = "[REQUIRED] BIOM file",
                    required = True)
parser.add_argument("-o",
                    action = "store",
                    dest = "outfile",
                    metavar = "FILE",
                    help = "[REQUIRED] PHYLOTYPE txt file",
                    required = True)
options = parser.parse_args()


#############################
# Import json formatted OTU #
#############################

import json, sys
jsondata = open(options.biominputfile)
biom = json.load(jsondata)

sampleSize = int(biom["shape"][1])
otus = int(biom["shape"][0])

taxonomies = []


#for i in range(len(biom["rows"])):
#	taxonomies.append(biom["rows"][i]["metadata"]["taxonomy"])
#print(taxonomies)


for i in range(len(biom["rows"])):
	taxonomies.append("; ".join(biom["rows"][i]["metadata"]["taxonomy"]))

sampleids = []
for i in range(len(biom["columns"])):
	sampleids.append(biom["columns"][i]["id"])

import numpy as np

# BIOM table into matrix
matrix = np.zeros(shape=(otus, sampleSize))
for i in biom["data"]:
	matrix[i[0], i[1]] = i[2]
totalCount = matrix.sum()

phylotypes = {}
for otu in range(otus):
	if taxonomies[otu] not in phylotypes:
		phylotypes[taxonomies[otu]] = matrix[otu]
	else:
		new_matrix = phylotypes[taxonomies[otu]] + matrix[otu]
		phylotypes[taxonomies[otu]] = new_matrix

outfile = open(options.outfile, "w")
outfile.write("#phylotype" + "\t" + "\t".join(sampleids) + "\t" + "taxonomy" + "\n")

count = 1
for k, a in phylotypes.items():
	outfile.write("phylotype_" + str(count) + "\t" + "\t".join(map(str,a.tolist())) + "\t" + k + "\n")
	count += 1

exit(0)
