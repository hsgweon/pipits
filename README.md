[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# PIPITS

### An automated pipeline for analyses of fungal internal transcribed spacer (ITS) sequences from the Illumina sequencing platform [(Gweon et al., 2015)](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.12399)

##### Shown to perform better than QIIME2 and Galaxy! - [See this paper](https://peerj.com/preprints/27019/)


## Updates

###### UPDATE (22 April 2019) - PIPITS 2.3
> - **PIPITS_PROCESS automatically downloads** UNITE database (the most recent version - UNITE version 02.02.2019), so there is no need to meddle with environment variables anymore. Just run commands and it will take care of the database issues. You can still use older database by the way using --unite option (see help by -h).
> - **PIPITS_FUNITS exploits multiple CPUs.** It's an experimental feature, so do use it with care. You can invoke to use multiple CPUs by using the usual ```-t NUMBER_OF_CPUS``` option.
> - Update PIPITS with ```conda update --channel bioconda --channel conda-forge --channel defaults pipits```
> then check you have version 2.3 installed by: ```conda list pipits```


<br>

## SYNOPSIS

### PIPITS:

- is an automated pipeline for analyses of fungal internal transcribed spacer (ITS) sequences from the Illumina sequencing platform.

- is designed to work best on POSIX systems (this essentiallly means it doesn't work in Windows).

- will need at least 4GB of RAM on your machine running 64bit Linux of mac OS, and it's been tested to be stable on Ubuntu 16.04, and macOS Mojave.

- Uses UNITE fungal db (and also comes with an option to run it against WARCUP fungal db).

- Just 4 commands, and you are good to go!

<br>

## A. Installing PIPITS

### A1. Install

It is recommended that you use a [conda](https://conda.io/) environment for running **PIPITS** to ensure that its dependencies are contained in this "sandbox". This meant that you don't mess with your existig system and you don't need to be the admin. Don't worry, it's easy - just type the following command. 

> EXPLANATION: install **PIPITS** and dependencies and create a Conda environment (here the environment is named "pipit_env" but you can choose any name you wish). PIPITS is exclusively compatible with Python3, so add "python=3.6" as below:

```shell
conda create -n pipits_env --channel bioconda --channel conda-forge --channel defaults python=3.6 pipits
```


### A2. Test PIPITS

The **PIPITS** is divided into three consequential parts:

1.  **Prepping raw sequences:** join, convert, quality filter etc.
2.  **Fungal ITS extraction:** remove conserved regions
3.  **Process** the reads to produce an OTU abundance table and the taxonomic assignment table for downstream analysis

Let's test it with a very small test dataset to ensure everything is set up correcly.
> EXPLANATION: Download & extract a test dataset

```sh
wget https://sourceforge.net/projects/pipits/files/PIPITS_TESTDATA/pipits_test.tar.gz -O pipits_test.tar.gz
tar xvfz pipits_test.tar.gz
```

> EXPLANATION: Enter into ```pipits_test``` directory

```sh
cd pipits_test
```

> EXPLANATION: Get into Conda environment you created above, and run the commands.

```sh
source activate pipits_env
pispino_createreadpairslist -i rawdata -o readpairslist.txt
pispino_seqprep -i rawdata -o out_seqprep -l readpairslist.txt
pipits_funits -i out_seqprep/prepped.fasta -o out_funits -x ITS2 -v -r
pipits_process -i out_funits/ITS.fasta -o out_process -v -r
```

<br>

## B. Running PIPITS

### B1. Sequence Preparation

Illumina reads are generally provided as demultiplexed FASTQ files where the Illumina software (BASESPACE) splits the reads into separate files, one for each barcode.

> EXPLANATION: [PISPINO](https://github.com/hsgweon/pispino) (originally part of **PIPITS**) provides a script called ```pispino_createreadpairslist``` which generates
a tab-delimited text file for all read-pairs from the directory containing your raw sequences

```sh
pispino_createreadpairslist -i rawdata -o readpairslist.txt
```

#### *Note*

1.  The command produces a tab-delimited file with three columns denoting forward and reverse read filenames and sample IDs for the pairs
2.  Prior to running the command, you need to ensure that the raw data are either uncompressed (“.fastq”), or compressed with bz2 or gz (“.fastq.bz2”, “.fastq.gz”). Sample IDs are taken from the first characters preceding an underscore (“_”) from each filename
3.  After running *pispino_createreadpairslist*, check the resulting file ("readpairslist.txt") to see correct filenames and desired sample IDs are listed in the resulting file ("readpairslist.txt"). No duplicate sample IDs are allowed.

> EXPLANATION: Once we have the list file ("readpairslist.txt"), we can then begin to "prepare" the sequences:

```sh
pispino_seqprep -i rawdata -o out_seqprep -l readpairslist.txt
```

#### *Note*

1.  Read-pairs are joined by examining the overlapping regions of sequences
2.  The resulting assembled reads are then quality filtered
3.  The header of each read is then relabelled with an index number followed by a Sample ID
4.  The resulting files are converted into FASTA and merged into a single file to produce the final output file "prepped.fasta" in the output directory


### B2. Fungal ITS extraction

The output from ***pipits_prep*** is taken as an input for this step. It is also mandatory to provide the script with which ITS subregion (i.e. ITS1 or ITS2) is to be extracted.

> EXPLANATION: the input file (indicated with "-i") is the resulting file from the previous step

```sh
pipits_funits -i out_seqprep/prepped.fasta -o out_funits -x ITS2
```

#### *Note*

1.  Selected subregion are extracted with [ITSx](http://microbiology.se/software/itsx/) and where necessary they are re-orientated to 5’ to 3’ direction. It is worth noting that ITSx uses HMMER3 (Mistry et al., 2013) to compare input sequences against a set of models built from a number of different subregions of ITS sequences found in various organisms. This makes ITSx an ideal tool for both extraction of
    desired ITS subregions as well as filtering for specific groups of
    organisms. It also means that while PIPITS has been created with the
    analysis of fungal amplicons in mind, it could be adapted for the
    analyses of other organism groups where ITS is used as a marker by
    changing the ITSx settings and reference databases
3.  Having extracted the subregion, sequences are re-inflated to reflect
    their original abundances. To date, the longest sequenceable reads
    from the Illumina technology are 300 bp x 2 which is not sufficient
    to sequence both ITS1 and ITS2 and to have an overlapping region to
    join them. For this reason the program supports only a single
    subregion extraction mode


### B3. Process sequences

> EXPLANATION: This is the final step involving clustering and assigning of taxonomy.

```sh
pipits_process -i out_funits/ITS.fasta -o out_process
```

#### *Note*

1.  Input sequences are dereplicated
2.  Short (< 100bp) and unique (singletons) are removed
3.  The sequences are clustered at 97% PID
4.  The resulting representative sequences for each cluster are
    subjected to chimera detection and removal
5.  The input sequences are mapped onto the chimera-free
    representative sequences at 97% PID
6.  The representatives are taxonomically assigned with RDP Classifier
    against the UNITE fungal ITS reference dataset
7.  The results are translated into two types of OTU abundance
    tables:
    - “**OTU abundance table**”, an OTU is defined as a cluster of reads
        with the user-defined threshold typically 97% sequence identity
        motivated by the expectation that these correspond approximately
        to species.
    - “**phylotype abundance table**”, an OTU is defined as a cluster of
        sequences binned into the same taxonomic assignments.
8.  If you have memory issues, try increasing the maximum memory with "--Xmx". For example, "--Xmx 4G".
9.	Once all finished, you can leave Conda environment by typeing
```
source deactivate
```

<br>

## C. Misc

### C1. Options

You can tweak parameters and there are several options for each of the above steps. To view them, type "-h" after each command.

```sh
pipits_prep -h
```

### C2. FUNGuild analysis

Run ***pipits_funguild.py*** on the resulting OTU table to have a reformatted version for [FUNGuild analysis](http://www.stbates.org/guilds/app.php). See their page for more detail.

```sh
pipits_funguild.py -i out_process/otu_table.txt -o out_process/otu_table_funguild.txt
```

<br>

## D. Citation
Please cite:
> Hyun S. Gweon, Anna Oliver, Joanne Taylor, Tim Booth, Melanie Gibbs, Daniel S. Read, Robert I. Griffiths and Karsten Schonrogge, PIPITS: an automated pipeline for analyses of fungal internal transcribed spacer sequences from the Illumina sequencing platform, Methods in Ecology and Evolution, DOI: 10.1111/2041-210X.12399

