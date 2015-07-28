#!/usr/bin/env python

import sys, os, argparse, subprocess, shutil, textwrap
import sgtk_SeqIO as SeqIO

import runcmd as rc
import tcolours as tc
import dependencies as pd

__author__ = "Hyun Soon Gweon"
__copyright__ = "Copyright 2015, The PIPITS Project"
__credits__ = ["Hyun Soon Gweon", "Anna Oliver", "Joanne Taylor", "Tim Booth", "Melanie Gibbs", "Daniel S. Read", "Robert I. Griffiths", "Karsten Schonrogge"]
__license__ = "GPL"
__maintainer__ = "Hyun Soon Gweon"
__email__ = "hyugwe@ceh.ac.uk"


def run(options):

    if not os.path.exists(options.outDir):
        os.mkdir(options.outDir)
    else:
        shutil.rmtree(options.outDir)
        os.mkdir(options.outDir)
    tmpDir = options.outDir + "/intermediate"
    if not os.path.exists(tmpDir):
        os.mkdir(tmpDir)


    # Logging
    import logging
    logger = logging.getLogger("pipits_funits")
    logger.setLevel(logging.DEBUG)

    streamLoggerFormatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", tc.HEADER + "%Y-%m-%d %H:%M:%S" + tc.ENDC)

    streamLogger = logging.StreamHandler()
    if options.verbose:
        streamLogger.setLevel(logging.DEBUG)
    else:
        streamLogger.setLevel(logging.INFO)
    streamLogger.setFormatter(streamLoggerFormatter)
    logger.addHandler(streamLogger)

    fileLoggerFormatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s", tc.HEADER + "%Y-%m-%d %H:%M:%S" + tc.ENDC)
    fileLogger = logging.FileHandler(options.outDir + "/log.txt", "w")
    fileLogger.setLevel(logging.DEBUG)
    fileLogger.setFormatter(fileLoggerFormatter)
    logger.addHandler(fileLogger)


    # Summary file
    summary_file = open(options.outDir + "/summary_pipits_funits.txt", "w")


    # Start! 
    logger.info(tc.OKBLUE + "PIPITS FUNITS started" + tc.ENDC)


    # Scripts
    EXE_DIR = os.path.dirname(os.path.realpath(__file__))
    PIPITS_SCRIPTS_DIR = EXE_DIR

    # Check integrity of the input file
    logger.info("Checking input FASTA for illegal characters")
    record = SeqIO.FastaParser(options.input)
    for i in record.keys():
        description = record[i].description
        if description.find(" ") != -1:
            logger.error("Error: \"  \" found in the headers. Please remove \" \" from headers in your FASTA file before proceeding to the next stage.")


    # For summary 1:
    logger.info("Counting input sequences")
    numberofsequences = 0
    cmd = " ".join(["grep \"^>\"", options.input, "|", "wc -l"])
    p = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE)
    numberofsequences += int(p.communicate()[0])
    p.wait()
    logger.info("\t" + tc.RED + "Number of input sequences: " + str(numberofsequences) + tc.ENDC)
    summary_file.write("Number of input sequences: " + str(numberofsequences) + "\n")


    # Dereplicate
    logger.info("Dereplicating sequences for efficiency")
    cmd = " ".join(["python", PIPITS_SCRIPTS_DIR + "/dereplicate_fasta.py", "-i", options.input, "-o", tmpDir + "/derep.fasta", "--cluster", tmpDir + "/derep.json"])
    rc.run_cmd(cmd, logger, options.verbose)


    # For summary 2:
    logger.debug("Counting dereplicated sequences")
    numberofsequences = 0
    cmd = " ".join(["grep \"^>\"", tmpDir + "/derep.fasta", "|", "wc -l"])
    p = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE)
    numberofsequences += int(p.communicate()[0])
    p.wait()
    logger.debug("\t" + tc.RED + "Number of dereplicated sequences: " + str(numberofsequences) + tc.ENDC)


    # Run ITSx. Chop reads into regions. Re-orientate where needed
    # ITSx always prints something to STDERR and outputs nothing to STDOUT, so need to supress stdout in non-verbose mode
    # Returncode is always 0 no matter what...
    # No way to tell whether it quits with an error or not other than by capturing STDERR with a phrase "FATAL ERROR" - not implemented
    logger.info("Extracting " + options.ITSx_subregion + " from sequences [ITSx]")
    cmd = " ".join([pd.ITSx, "-i", tmpDir + "/derep.fasta", "-o", tmpDir + "/derep", "--preserve", "T", "-t", "F", "--cpu", options.threads, "--save_regions", options.ITSx_subregion])
    rc.run_cmd_ITSx(cmd, logger, options.verbose)


    # Removing short sequences (<100bp)
    logger.info("Removing sequences below < 100bp")
    cmd = " ".join(["python", PIPITS_SCRIPTS_DIR + "/fasta_filter_by_length.py", 
                    "-i", tmpDir + "/derep." + options.ITSx_subregion + ".fasta", 
                    "-o", tmpDir + "/derep." + options.ITSx_subregion + ".sizefiltered.fasta", 
                    "-l 100"])
    rc.run_cmd(cmd, logger, options.verbose)


    # Re-inflate
    logger.info("Re-inflating sequences")
    cmd = " ".join(["python", PIPITS_SCRIPTS_DIR + "/inflate_fasta.py", 
                    "-i", tmpDir + "/derep." + options.ITSx_subregion + ".sizefiltered.fasta", 
                    "-o", options.outDir + "/ITS.fasta", "--cluster", tmpDir + "/derep.json"])
    rc.run_cmd(cmd, logger, options.verbose)


    # Count number of ITS
    logger.info("Counting sequences after re-inflation")
    numberofsequences = 0

    cmd = " ".join(["grep \"^>\"", options.outDir + "/ITS.fasta", "|", "wc -l"])
    p = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    numberofsequences = int(p.communicate()[0])
    p.wait()

    if numberofsequences == 0:
        logger.info(tc.RED + "\tNumber of sequences with ITS subregion: " + str(numberofsequences) + tc.ENDC)
        logger.info(tc.RED + "Have you chosen the right subregion? Exiting as no sequences to process." + tc.ENDC)
        summary_file.write("Number of sequences with ITS subregion: " + str(numberofsequences) + "\n")
        exit(1)
    else:
        logger.info(tc.RED + "\tNumber of sequences with ITS subregion: " + str(numberofsequences) + tc.ENDC)
        summary_file.write("Number of sequences with ITS subregion: " + str(numberofsequences) + "\n")


    '''
    # Concatenating ITS1 and ITS2
    logger.info("Concatenating ITS1 and ITS2 ...")
    cmd = " ".join(["python", PIPITS_SCRIPTS_DIR + "/concatenate_fasta.py", 
                    "-1", options.outDir + "/ITS1.fasta" , 
                    "-2", options.outDir + "/ITS2.fasta",
                    "-o", options.outDir + "/ITS.fasta"])
    rc.run_cmd(cmd, logger, options.verbose)
    logger.info("Concatenating ITS1 and ITS2 " + tc.OKGREEN + "(Done)" + tc.ENDC)
    '''


    # Finally move and delete tmp
    if options.remove:
        logger.info("Cleaning temporary directory")
        shutil.move(tmpDir + "/derep.summary.txt", options.outDir + "/ITSx_summary.txt")
        shutil.rmtree(tmpDir)


    logger.info(tc.OKBLUE + "PIPITS FUNITS ended successfully. \"" + "ITS.fasta" + "\" created in \"" + options.outDir + "\"" + tc.ENDC)
    logger.info(tc.OKYELLOW + "Next Step: PIPITS PROCESS [ Suggestion: pipits process -i " + options.outDir + "/" + "ITS.fasta -o out_process ]" + tc.ENDC)
    print("")
    summary_file.close()

