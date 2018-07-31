[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# PIPITS*

### An automated pipeline for analyses of fungal internal transcribed spacer (ITS) sequences from the Illumina sequencing platform [(Gweon et al., 2015)](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.12399)

#### * Shown to perform better than QIIME2 and Galaxy! - [See this paper](https://peerj.com/preprints/27019/)


## Updates

###### UPDATE (31 July 2018)
> - PIPITS 2.2 - Changes made with the way PIPITS_FUNITS dereplicates and rereplicates sequences. Now using VSEARCH to dereplicate and PIPITS\_REREPLICATE to rereplicate. This shouldn't affect the resulting OTU table.
> - Update PIPITS with:
> ```
> conda update pipits
> ```
> then check you have version 2.2 installed by:
> ```
> conda list pipits
> ```


###### UPDATE (13 June 2018)
> - Bug fix update released: PIPITS 2.1: **There has been a major bug in PIPITS_FUNITS affecting PIPITS2** (PIPITS1 isn't affected by this bug). Due to this bug, sequences were not inflated properly after dereplication (for speed). If you have used PIPITS2, then you WILL need to update PIPITS with:
> ```
> conda update pipits
> ```
> Then check you have pipits 2.1 installed by:
> ```
> conda list pipits
> ```


###### UPDATE (29 March 2018)

> - **PIPITS 2.0** is exclusively available for Python3. This shouldn't be a problem to most of us as we will all be using Conda environment (I hope!). See the instruction for more detail.

###### UPDATE (22 March 2018)

> - **PIPITS 2.0** is now available for download from the [Bioconda](https://bioconda.github.io/index.html) channel. This means that PIPITS can be installed very easily now. See below for more information
> - It now supports Mac OS as well as Linux
> - Because of *pipits\_prep*'s generic nature, it is now called *pispino\_seqprep*, and moved out of **PIPITS**, and now is part of [**PISPINO**](https://github.com/hsgweon/pispino) package. **PISPINO** will be installed automatically by **PIPITS**
> - The recent UNITE DB (01.12.2017) has been trained and uploaded for **PIPITS**

###### UPDATE (5 July 2017)

> - PEAR is NO LONGER the default joiner method: PEAR now requires an academic licence and you cannot directly download it from its new homepage any more. I am not a big fan of tools with a restricted licence so I have decided to move on from PEAR and instead use VSEARCH's relatively new joining protocol. VSEARCH developers have openly said that their protocol is based on the algorithm used by PEAR and indeed I have tested VSEARCH on a very large ITS datasets and it was shown to perform in a very similar or almost exactly the same way as PEAR (I am considering whether it would be worth publishing a paper on these different joining methods and the implications). If you still prefer to use PEAR and have managed to get PEAR downloaded, then you can always choose to use PEAR by specifying "--joiner_method PEAR" in ***pipits\_prep***. Otherwise, VSEARCH seems to be a perfect replacement - fast, reliable and free
> - New UNITE dataset (2017-06-28) available.
> - PIPITS is purposefully designed to give limited options compared to other tools, so you have more time to analyse the OTU Table than worry about tweaking parameters and dependencies. If you have any suggestions , do let me know by commenting on Issues or email me, and I will try my best to get it implemented ASAP.

<br>

## SYNOPSIS

### PIPITS:

- is an automated pipeline for analyses of fungal internal transcribed spacer (ITS) sequences from the Illumina sequencing platform

- is designed to work best on POSIX systems (this essentiallly means it doesn't work in Windows)

- will need at least 4GB of RAM on your machine running 64bit Linux of mac OS, and it's been tested to be stable on Ubuntu 16.04, and macOS High Sierra

- is compatible on PYTHON3

<br>

## A. Setting up PIPITS


### A1. Prerequisite: set up conda channels (only for the first time)

> add the [Bioconda](https://bioconda.github.io/index.html) channel and others which [Bioconda](https://bioconda.github.io/index.html) depends on. It is important to add them in this order (this needs to be done only for the first time)

```shell
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

### A2. Install and create a Conda environment for PIPITS

It is recommended that you use a [Conda](https://conda.io/) environment for running **PIPITS** to ensure that its dependencies are contained in this "sandbox". Don't worry, it's easy - just type the following command

> install **PIPITS** and dependencies and create a Conda environment (here the environment is named "pipit_env" but you can choose any name you wish). PIPITS is exclusively compatible with Python3, so add "python=3.6" as below:

```shell
conda create -n pipits_env python=3.6 pipits
```

### A3. Reference datasets

There are three reference datasets to download:

**UNITE fungal ITS reference trained dataset. UNITE 7.2 (2017-12-01)** - [file](https://sourceforge.net/projects/pipits/files/PIPITS_DB/UNITE_retrained_01.12.2017.tar.gz)

> Create a directory for the datasets

```sh
mkdir -p $HOME/pipits/refdb
```
> We provide a trained UNITE fungal data. Download, save and extract it to an appropriate directory (e.g. $HOME/pipits/refdb). On macOS, wget is not available but you can install it - see [here](https://gist.github.com/shrayasr/8206257)

```
cd $HOME/pipits/refdb
wget https://sourceforge.net/projects/pipits/files/PIPITS_DB/UNITE_retrained_01.12.2017.tar.gz -O UNITE_retrained_01.12.2017.tar.gz
rm -rf UNITE_retrained
tar xvfz UNITE_retrained_01.12.2017.tar.gz
```

**UNITE UCHIME reference dataset** - [file](https://unite.ut.ee/sh_files/uchime_reference_dataset_28.06.2017.zip)

> Download UNITE UCHIME reference dataset for chimera removal.

```sh
cd $HOME/pipits/refdb
wget https://unite.ut.ee/sh_files/uchime_reference_dataset_28.06.2017.zip -O uchime_reference_dataset_28.06.2017.zip
rm -rf uchime_reference_dataset_28.06.2017
unzip uchime_reference_dataset_28.06.2017.zip
```

**Warcup ITS reference trained dataset** - [file](https://sourceforge.net/projects/pipits/files/warcup_retrained_V2.tar.gz)

> **PIPITS** supports Warcup ITS reference training dataset. By specifying "--warcup" when running *pipits_process*, **PIPITS** will create an additional Warcup classified OTU table.

```sh
cd $HOME/pipits/refdb
wget https://sourceforge.net/projects/pipits/files/warcup_retrained_V2.tar.gz
rm -rf warcup_retrained_V2
tar xvfz warcup_retrained_V2.tar.gz
```


### A4. Set environment variables


We need to set some environmental variables to let **PIPITS** know where these reference datasets are. Add the following lines in your system's profile file. Ubuntu's default profile file is "$HOME/.bashrc", and on mac OS, it is "$HOME/.bash_profile"


> Open the profile file with a text editor, and add:

    export PIPITS_UNITE_RETRAINED_DIR=$HOME/pipits/refdb/UNITE_retrained
    export PIPITS_UNITE_REFERENCE_DATA_CHIMERA=$HOME/pipits/refdb/uchime_reference_dataset_28.06.2017/uchime_reference_dataset_28.06.2017.fasta
    export PIPITS_WARCUP_RETRAINED_DIR=$HOME/pipits/refdb/warcup_retrained_V2

> Close and re-open the terminal, alternatively you can:

```sh
source $HOME/.bashrc
```

<br>


## B. Running PIPITS

The **PIPITS** is divided into three consequential parts:

1.  **Sequence preparation** to join, convert, quality filter etc.
2.  **Fungal ITS extraction** to remove conserved regions
3.  **Process** the reads to produce an OTU abundance table and the taxonomic assignment table for downstream analysis

Test it with a very small test dataset to ensure everything is set up correcly.
> Download & extract a test dataset

```sh
wget http://sourceforge.net/projects/pipits/files/PIPITS_TESTDATA/rawdata.tar.gz -O rawdata.tar.gz
tar xvfz rawdata.tar.gz
```
>Get into Conda environment you created above

```sh
source activate pipits_env
```


### B1. Sequence Preparation

Illumina reads are generally provided as demultiplexed FASTQ files where the Illumina software (BASESPACE) splits the reads into separate files, one for each barcode.

> [PISPINO](https://github.com/hsgweon/pispino) (originally part of **PIPITS**) provides a script called *pispino_createreadpairslist* which generates
a tab-delimited text file for all read-pairs from the directory containing your raw sequences

```sh
pispino_createreadpairslist -i rawdata -o readpairslist.txt
```

#### *Note*

1.  The command produces a tab-delimited file with three columns denoting forward and reverse read filenames and sample IDs for the pairs
2.  Prior to running the command, you need to ensure that the raw data are either uncompressed (“.fastq”), or compressed with bz2 or gz (“.fastq.bz2”, “.fastq.gz”). Sample IDs are taken from the first characters preceding an underscore (“_”) from each filename
3.  After running *pispino_createreadpairslist*, check the resulting file ("readpairslist.txt") to see correct filenames and desired sample IDs are listed in the resulting file ("readpairslist.txt"). No duplicate sample IDs are allowed.

> Once we have the list file ("readpairslist.txt"), we can then begin to "prepare" the sequences:

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

> the input file (indicated with "-i") is the resulting file from the previous step

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

> This is the final step involving clustering and assigning of taxonomy.

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

