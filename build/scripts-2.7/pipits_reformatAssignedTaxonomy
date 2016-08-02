#!/usr/bin/python

############################################################
# Argument Options

import argparse
parser = argparse.ArgumentParser("Re-format taxonomy assignment output from RDP-CLASSIFIER.")
parser.add_argument("-i, --in",
                    action = "store",
                    dest = "input",
                    metavar = "input",
                    help = "[REQUIRED] taxonomy assignment output from RDP-CLASSIFIER",
                    required = True)
parser.add_argument("-o, --out",
                    action = "store",
                    dest = "output",
                    metavar = "output",
                    help = "[REQUIRED] reformatted taxonomy assignment file",
                    required = True)
parser.add_argument("-c",
                    action = "store",
                    dest = "confidence",
                    metavar = "confidence",
                    help = "[REQUIRED] Minimum confidence to record an assignment",
                    required = True)
options = parser.parse_args()

############################################################

import sys

THRESHOLD = float(options.confidence)

handle_input = open(options.input, "rU")
handle_output = open(options.output, "w")


for line in handle_input:

    e = line.rstrip().split("\t")
    taxonomy = []
    taxonomy_conf = []

    for i in range(5, len(e), 3):

        tax_level = e[i+1]
        tax_nomenclature = e[i].split("|")[-1].replace(" ", "_")
        tax_confidence = float(e[i+2])

        if tax_level == "domain" or tax_level == "kingdom":
            tax_nomenclature = "k__" + tax_nomenclature
        elif tax_level == "phylum":
            tax_nomenclature = "p__" + tax_nomenclature
        elif tax_level == "subphylum":
            tax_nomenclature = "subp__" + tax_nomenclature
        elif tax_level == "class":
            tax_nomenclature = "c__" + tax_nomenclature
        elif tax_level == "subclass":
            tax_nomenclature = "subc__" + tax_nomenclature
        elif tax_level == "order":
            tax_nomenclature = "o__" + tax_nomenclature
        elif tax_level == "family":
            tax_nomenclature = "f__" + tax_nomenclature
        elif tax_level == "genus":
            tax_nomenclature = "g__" + tax_nomenclature
        elif tax_level == "species":
            tax_nomenclature = "s__" + tax_nomenclature
        else:
            print(tax_level, "Error in RDP Classifier produced output: not a valid taxonomic level.")
            exit(1)

        taxonomy.append(tax_nomenclature)
        taxonomy_conf.append(tax_confidence)

    '''
    Kingdom = "k__" + e[5]
    Kingdom_conf = float(e[7])
    Phylum = "p__" + e[8].split("|")[-1]
    Phylum_conf = float(e[10])
    Class = "c__" + e[11].split("|")[-1]
    Class_conf = float(e[13])
    Order = "o__" + e[14].split("|")[-1]
    Order_conf = float(e[16])
    Family = "f__" + e[17].split("|")[-1]
    Family_conf = float(e[19])
    Genus = "g__" + e[20].split("|")[-1]
    Genus_conf = float(e[22])
    Species = "s__" + e[23].split("|")[-1]
    Species_conf = float(e[25])
    '''

    taxonomy_filtered = []
    taxonomy_conf_filtered = []

    for i in range(len(taxonomy)):
        if taxonomy_conf[i] >= THRESHOLD:
            taxonomy_filtered.append(taxonomy[i])
            taxonomy_conf_filtered.append(taxonomy_conf[i])
        else:
            break

    if len(taxonomy_conf_filtered) == 0:
        taxonomy_filtered.append("Unassignable")
        taxonomy_conf_filtered.append(1.0)

    # print(";".join(taxonomy_filtered))
    # print(";".join(map(str, taxonomy_conf_filtered)))

    handle_output.write(e[0] + "\t")
    handle_output.write("; ".join(taxonomy_filtered) + "\t")
    handle_output.write(str(taxonomy_conf_filtered[-1]) + "\n")


handle_input.close()
handle_output.close()

