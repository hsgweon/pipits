#!/usr/bin/env python

import dependencies as pd
import subprocess, textwrap, os

__author__ = "Hyun Soon Gweon"
__copyright__ = "Copyright 2015, The PIPITS Project"
__credits__ = ["Hyun Soon Gweon", "Anna Oliver", "Joanne Taylor", "Tim Booth", "Melanie Gibbs", "Daniel S. Read", "Robert I. Griffiths", "Karsten Schonrogge"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Hyun Soon Gweon"
__email__ = "hyugwe@ceh.ac.uk"


def run(options):

    # http://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/ 
    cmd = " ".join(["java", 
                    "-jar", pd.RDP_CLASSIFIER_JAR, 
                    "train", 
                    "-o", options.outDir,
                    "-s", options.fasta, 
                    "-t", options.taxonomy])
    p = subprocess.Popen(cmd, shell=True)
    print("Please be patient. This process can take awhile.")
    p.wait()
    if p.returncode != 0:
        exit(1)

    outfile = open(options.outDir + "/rRNAClassifier.properties", "w")
    properties = "bergeyTree=bergeyTrainingTree.xml\nprobabilityList=genus_wordConditionalProbList.txt\nprobabilityIndex=wordConditionalProbIndexArr.txt\nwordPrior=logWordPrior.txt\nclassifierVersion=RDP Naive Bayesian rRNA Classifier Version\n"
    outfile.write(properties)

    print("Done, now edit the \"pipits.config\" file to let PIPITS know the location of the output directory (\"" + options.outDir + "\").")
    print("")
    print("Suggestion:")
    print("\tUNITE_TRAINED_SET_DIR = " + os.getcwd() + "/unite_retrained")
    print("")
    exit(0)


