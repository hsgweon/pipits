#!/usr/bin/env python

import sys

############################################################
# Argument Options
# Header needs to have phylum -> species all separated by '_'

import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description = "[] indicates optional input (order unimportant)\n\nDereplicates FASTA", formatter_class = RawTextHelpFormatter)
parser.add_argument("-i",
                    action = "store",
                    dest = "input",
                    metavar = "input",
                    help = "[REQUIRED] FASTA",
                    required = True)
parser.add_argument("-o",
                    action = "store",
                    dest = "output",
                    metavar = "output",
                    help = "[REQUIRED] FASTA",
                    required = True)
parser.add_argument("--cluster",
                    action = "store",
                    dest = "cluster",
                    metavar = "cluster",
                    help = "[REQUIRED] JSON",
                    required = True)
options = parser.parse_args()

############################################################

from itertools import groupby

isheader = lambda x : x.startswith('>')

d = {}
with open(options.input, "r") as infile:
    header = None
    for h, group in groupby(infile, isheader):
        if h:
            header = group.next().strip()[1:]
        else:
            seq = ''.join(line.strip() for line in group)
            if seq not in d:
                d[seq] = [header]
            else:
                d[seq].append(header)

d4json = {}
with open(options.output, 'w') as outfile:
    
    for key, value in d.iteritems():
        outfile.write(">" + value[0] + ";size=" + str(len(value)) + "\n" + key + "\n")
        d4json[value[0]] = [i for i in value[1:]]

import json
with open(options.cluster, "w") as cluster:
    json.dump(d4json, cluster)
