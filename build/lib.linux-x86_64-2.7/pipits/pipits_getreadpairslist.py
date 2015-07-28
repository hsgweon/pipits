#!/usr/bin/env python

import sys, os, subprocess, ConfigParser, shutil
import tcolours as tc
from time import strftime

__author__ = "Hyun Soon Gweon"
__copyright__ = "Copyright 2015, The PIPITS Project"
__credits__ = ["Hyun Soon Gweon", "Anna Oliver", "Joanne Taylor", "Tim Booth", "Melanie Gibbs", "Daniel S. Read", "Robert I. Griffiths", "Karsten Schonrogge"]
__license__ = "GPL"
__maintainer__ = "Hyun Soon Gweon"
__email__ = "hyugwe@ceh.ac.uk"


def run(options):

    print(tc.OKGREEN + "Generating read-pair list file from the input directory...")

    fastqs = []
    for file in os.listdir(options.dataDir):
        if file.endswith(".fastq.gz") or file.endswith(".fastq.bz2") or file.endswith(".fastq"):
            fastqs.append(file)

    if len(fastqs) % 2 != 0:
        logger.error("There are missing pair(s) in the Illumina sequences. Check your files and labelling")
        exit(1)

    fastqs_l = []
    fastqs_f = []
    fastqs_r = []

    coin = True
    for fastq in sorted(fastqs):
        if coin == True:
            fastqs_f.append(fastq)
        else:
            fastqs_r.append(fastq)
        coin = not coin

    for i in range(len(fastqs_f)):
        if fastqs_f[i].split("_")[0] != fastqs_r[i].split("_")[0]:
            logger.error("Problem with labelling the files.")
            exit(1)
        fastqs_l.append(fastqs_f[i].split("_")[0])

    if not options.output:
        outfile_name = options.dataDir + "_readpairslist.txt"
    else:
        outfile_name = options.output

    outfile_fastqslist = open(outfile_name, "w")
    outfile_fastqslist.write("# Lines beginning with \"#\" is ignored. \n")
    outfile_fastqslist.write("# SampleID\tFilename for forward-reads\tFilename for reverse-reads\n")
    count = 1
    for i in range(len(fastqs_f)):
        label = ""
        if options.label_add_c:
            label = fastqs_l[i] + options.label_add_c
        elif options.label_add_reindex_c:
            label = options.label_add_reindex_c + str(count).zfill(3)
        else:
            label = fastqs_l[i]
        count += 1
        outfile_fastqslist.write(label + "\t" + fastqs_f[i] + "\t" + fastqs_r[i] + "\n")
    outfile_fastqslist.close()

    print(tc.OKYELLOW + "Done. \"" + outfile_name + "\" created.")
