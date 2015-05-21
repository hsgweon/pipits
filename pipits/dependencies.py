#!/usr/bin/env python

import subprocess, os, ConfigParser

__author__ = "Hyun Soon Gweon"
__copyright__ = "Copyright 2015, The PIPITS Project"
__credits__ = ["Hyun Soon Gweon", "Anna Oliver", "Joanne Taylor", "Tim Booth", "Melanie Gibbs", "Daniel S. Read", "Robert I. Griffiths", "Karsten Schonrogge"]
__license__ = "GPL"
__maintainer__ = "Hyun Soon Gweon"
__email__ = "hyugwe@ceh.ac.uk"

def init(config_file):

    config = ConfigParser.ConfigParser()
    try:
        with open(config_file):
            config.read(config_file)
    except IOError:
        print("Error: \"" + config_file + "\" not found.")
        exit(1)

    global FASTX_FASTQ_QUALITY_FILTER
    try:
        FASTX_FASTQ_QUALITY_FILTER = config.get("DEPENDENCIES", "FASTX_FASTQ_QUALITY_FILTER")
    except ConfigParser.NoOptionError:
        FASTX_FASTQ_QUALITY_FILTER = "fastq_quality_filter"

    global FASTX_FASTQ_TO_FASTA
    try:
        FASTX_FASTQ_TO_FASTA = config.get("DEPENDENCIES", "FASTX_FASTQ_TO_FASTA") 
    except ConfigParser.NoOptionError:
        FASTX_FASTQ_TO_FASTA = "fastq_to_fasta"

    global PEAR
    try:
        PEAR = config.get("DEPENDENCIES", "PEAR")
    except ConfigParser.NoOptionError:
        PEAR = "pear"

    global ITSx
    try:
        ITSx = config.get("DEPENDENCIES", "ITSx")
    except ConfigParser.NoOptionError:
        ITSx = "ITSx"

    global VSEARCH
    try:
        VSEARCH = config.get("DEPENDENCIES", "VSEARCH")
    except ConfigParser.NoOptionError:
        VSEARCH = "vsearch"

    global BIOM
    try:
        BIOM = config.get("DEPENDENCIES", "BIOM")
    except ConfigParser.NoOptionError:
        BIOM = "biom"

    global RDP_CLASSIFIER_JAR
    RDP_CLASSIFIER_JAR = config.get("DEPENDENCIES", "RDP_CLASSIFIER_JAR")
    
    global UNITE_REFERENCE_DATA_CHIMERA
    UNITE_REFERENCE_DATA_CHIMERA = config.get("DB", "UNITE_REFERENCE_DATA_CHIMERA")
    
    global UNITE_RETRAINED_DIR
    UNITE_RETRAINED_DIR = config.get("DB", "UNITE_RETRAINED_DIR")
    if UNITE_RETRAINED_DIR.endswith("/"): 
        UNITE_RETRAINED_DIR = UNITE_RETRAINED_DIR[:-1]

'''
    # Some crude checks
    for i in [UNITE_REFERENCE_DATA_CHIMERA]:
        p = subprocess.Popen("ls " + i, shell=True, stdout=subprocess.PIPE)
        p.wait()
        if p.returncode!= 0:
            print("Cannot find " + i + ". Make sure it is correcly set up and specified in pipits.config")
            exit(1)
'''

# This is for well-behaving third party tools i.e. ones that give proper STDOUT and STDERR
def run_cmd(command, logger, verbose):
    logger.debug(command)
    FNULL = open(os.devnull, 'w')
    if verbose:
        p = subprocess.Popen(command, shell=True)
    else:
        p = subprocess.Popen(command, shell=True, stdout=FNULL)
    p.wait()
    FNULL.close()
    if p.returncode != 0:
        logger.error("None zero returncode: " + command)
        exit(1)


# Run ITSx. Chop reads into regions. Re-orientate where needed
# ITSx always prints something to STDERR and outputs nothing to STDOUT, so need to supress stdout in non-verbose mode
# Returncode is always 0 no matter what... so way to tell whether it quits with an error or not other than by capturing STDERR with a phrase "FATAL ERROR" - not implemented 
def run_cmd_ITSx(command, logger, verbose):
    logger.debug(command)
    FNULL = open(os.devnull, 'w')
    if verbose:
        p = subprocess.Popen(command, shell=True)
    else:
        p = subprocess.Popen(command, shell=True, stderr=FNULL)
    p.wait()
    FNULL.close()


# VSEARCH outputs copyright info and licence to STDOUT; and the running outputs to STDERRDATA
def run_cmd_VSEARCH(command, logger, verbose):
    logger.debug(command)
    FNULL = open(os.devnull, 'w')
    if verbose:
        p = subprocess.Popen(command, shell=True)
    else:
        p = subprocess.Popen(command, shell=True, stdout=FNULL, stderr=FNULL)
    p.wait()
    FNULL.close()
    if p.returncode != 0:
        logger.error("None zero returncode: " + command)
        exit(1)

