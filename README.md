Introduction
============

- PIPITS is an automated pipeline for analyses of fungal internal transcribed spacer (ITS) sequences from the Illumina sequencing platform.

- PIPITS is designed to work best on Bio-Linux (<http://environmentalomics.org/bio-linux/>) and Ubuntu. Unfortunately, it's NOT supported on Windows or a Mac

- If you are using Bio-Linux, most of the dependencies are already on Bio-Linux. Otherwise, you will have to set up the dependencies yourself. If you are using Ubuntu, then instructions on how to set up dependencies are described below (1.8).

- If you only have a Windows machine, we recommend installing VirtualBox (free) to run Bio-Linux (<http://environmentalomics.org/bio-linux/>).


1 Setting up PIPITS
===================

1.1 Download and install
------------------------

Download the latest package from github release (<https://github.com/hsgweon/pipits/releases>) or simply copy and paste the following:

```sh
cd ~
wget https://github.com/hsgweon/pipits/releases/download/1.3.4/pipits-1.3.4.tar.gz
tar xvfz pipits-1.3.4.tar.gz
```

Then enter into the created directory and install the package with (ignore errors/warnings):

```sh
cd pipits-1.3.4
python setup.py clean --all
python setup.py install --prefix=$HOME/pipits
```

This creates a "pipits" directory in your $HOME and we will be installing pipits and some of its dependencies into this directory.
Of course if you are familiar with Linux, you are free to choose whichever method to suit your skill and taste!


1.2 Dependencies
----------------

*Download and install dependencies*

PIPITS depends on a number of external dependencies which need to be downloaded and installed. 

We advise you to use Ubuntu 16.04 (xenial) or above as all of the dependencies are now available as a Debian package! If you are, then just follow the instructions below. (If you are using Ubuntu 14.04, then we *strongly* recommend using Bio-Linux packages rather than installing dependencies yourself. See 1.8 below for the detailed instruction on how you do this.)


1. **BIOM-FORMAT version 2.x** (<http://biom-format.org/>)

   *Available as a Debian package.*

   ```
   sudo apt install biom-format-tools
   ```

2. **FAST-X tools** (<http://hannonlab.cshl.edu/fastx_toolkit>)

   *Available as a Debian package.*

   ```
   sudo apt install fastx-toolkit
   ```

3. **VSEARCH** (<https://github.com/torognes/vsearch>)

   Download a recent version of vsearch (check by visiting their GITHUB page and change the version name below accordingly. v1.10.2 should be stable for PIPITS.) Debian package is a little out of date.

   ```
   cd $HOME/pipits
   wget https://github.com/torognes/vsearch/releases/download/v1.11.2/vsearch-1.11.2-linux-x86_64.tar.gz
   tar xvfz vsearch-1.11.2-linux-x86_64.tar.gz
   ln -s $HOME/pipits/vsearch-1.11.2-linux-x86_64/bin/vsearch bin/vsearch
   ```

4. **ITSx** (<http://microbiology.se/software/itsx>) N.B. ITSx requires HMMER3

    ```sh
    cd $HOME/pipits
    wget http://microbiology.se/sw/ITSx_1.0.11.tar.gz
    tar xvfz ITSx_1.0.11.tar.gz
    ln -s $HOME/pipits/ITSx_1.0.11/ITSx bin/ITSx
    ln -s $HOME/pipits/ITSx_1.0.11/ITSx_db bin/ITSx_db
    ```
5. **PEAR** (<http://sco.h-its.org/exelixis/web/software/pear>) - N.B. PEAR prohibits commercial use of the code. See its webpage for detail.
 
    ```sh
    cd $HOME/pipits
    wget http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.6-bin-64.tar.gz
    tar xvfz pear-0.9.6-bin-64.tar.gz
    ln -s $HOME/pipits/pear-0.9.6-bin-64/pear-0.9.6-bin-64 bin/pear
    ```

6. **RDP Classifier 2.9 or above** (<http://sourceforge.net/projects/rdp-classifier>) - N.B. RDP Classifier comes with a jar file.
   
    ```sh 
    cd $HOME/pipits
    wget http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.11.zip
    unzip rdp_classifier_2.11.zip
    ln -s rdp_classifier_2.11/dist/classifier.jar ./classifier.jar
    ```

7. **HMMER3** (<http://hmmer.janelia.org/download.html>)

   *Available as a Debian package.*

   ```
   sudo apt install hmmer
   ```

8. **Java compatible Runtime**

   ```
   sudo apt install default-jre
   ```

9. **numpy**

   ```
   sudo apt install python-numpy
   ```


1.3 Reference datasets
----------------------

There are two reference datasets to download:

1. **UNITE fungal ITS reference trained dataset**

   We now provide trained UNITE fungal data (processed and trained for PIPITS).
   Please download this data (<http://sourceforge.net/projects/pipits/files/UNITE_retrained_31.01.2016.tar.gz>), save and extract it to an appropriate directory (e.g. $HOME/pipits/refdb).

   ```sh
   mkdir -p $HOME/pipits/refdb
   cd $HOME/pipits/refdb
   wget http://sourceforge.net/projects/pipits/files/UNITE_retrained_31.01.2016.tar.gz
   tar xvfz UNITE_retrained_31.01.2016.tar.gz
   ```

2. **UNITE UCHIME reference dataset**

   We also need to download UNITE UCHIME reference dataset for chimera removal. Download it from UNITE repository (<http://unite.ut.ee/repository.php>).

   ```sh
   mkdir -p $HOME/pipits/refdb
   cd $HOME/pipits/refdb
   wget https://unite.ut.ee/sh_files/uchime_reference_dataset_01.01.2016.zip
   unzip uchime_reference_dataset_01.01.2016.zip
   ```

3. **(OPTIONAL) Warcup ITS reference trained dataset**

   PIPITS supports Warcup ITS reference training dataset. By specifying "--warcup" when running pipits_process, PIPITS will also create a Warcup classified OTU table.

   ```sh
   mkdir -p $HOME/pipits/refdb
   cd $HOME/pipits/refdb
   wget http://sourceforge.net/projects/pipits/files/warcup_retrained_08.07.2014.tar.gz
   tar xvfz warcup_retrained_08.07.2014.tar.gz
   ```


1.4 Set PATH and ENVIRONMENT VARIABLE
-------------------------------------

Now we will make sure executables and modules are visible to the shell by existing in the search PATH. Also we will set some environment variables. Assuming UBUNTU is your system, this can be easily achieved by adding the following lines at the end of your profile file. The name of your profile file will depend on which shell your system is using. You can check which shell your system is using by typing *echo $SHELL* in your terminal. If it says /bin/bash, then your profile file is "~/.bashrc". N.B. UBUNTU's default shell is bash while Bio-Linux's default shell is zsh.

Open "~/.bashrc" or "~/.zshrc" (depending on which shell you are using) with a text editor such as gedit.



```sh    
gedit ~/.bashrc
(or gedit ~/.zshrc For Ubuntu 14.04 or below)
```

And then add the following lines at the end of the file:

    export PATH=$HOME/pipits/bin:$PATH
    export PYTHONPATH=$HOME/pipits/lib/python2.7/site-packages:$PYTHONPATH
    export PIPITS_UNITE_REFERENCE_DATA_CHIMERA=$HOME/pipits/refdb/uchime_reference_dataset_01.01.2016/uchime_reference_dataset_01.01.2016.fasta
    export PIPITS_UNITE_RETRAINED_DIR=$HOME/pipits/refdb/UNITE_retrained
    export PIPITS_WARCUP_RETRAINED_DIR=$HOME/pipits/refdb/warcup_retrained_08.07.2014
    export PIPITS_RDP_CLASSIFIER_JAR=$HOME/pipits/classifier.jar

Then type (or alternatively close and re-open the terminal):

```sh
source ~/.bashrc
(or source ~/.zshrc for Ubuntu 14.04 or below)
```

1.5 Re-HMMPressing
------------------

Once you downloaded and installed ITSx, we recommend re-HMMPRESSing the HMM profiles as the HMMPRESS'ed profiles may not be compatible with the version of the HMMER3 you installed:

```sh
cd $HOME/pipits/ITSx_1.0.11/ITSx_db/HMMs
rm -f *.hmm.*
echo *.hmm | xargs -n1 hmmpress
```

1.6 Test Dependencies and PIPITS
--------------------------------

When you have successfully installed these, check if they are *indeed* successfully installed by running each applications. If you get an error,
check to see if you have followed the instruction carefully.

```sh
$ biom
$ fastq_to_fasta -h
$ vsearch
$ ITSx -h
$ pear
$ ls $HOME/pipits/classifier.jar
$ hmmpress -h
```

Ok, let's test if PIPITS is all setup. Open up the very first original PIPITS which you downloaded. Please change X.X.X in the command below to the version of PIPITS you downloaded. Note that if you encounter memory issues with JAVA, try increasing the memory with "--Xmx" option.

```sh
cd $HOME/pipits-X.X.X/test_data
pipits_getreadpairslist -i rawdata
pipits_prep -i rawdata
pipits_funits -i pipits_prep/prepped.fasta -x ITS2 
pipits_process -i pipits_funits/ITS.fasta --Xmx 2G

(pipits_process -i pipits_funits/ITS.fasta --Xmx 2G -o pipits_process_with_warcup --warcup) If you want additional OTU table with Warcup classification.
```

Ensure everything works and you don't get an error message.


1.8 (Misc) For Ubuntu 14.04 users, installing dependencies using Bio-Linux packages
------------------------------------------------------------------------------------

Some of the dependencies are available as a Bio-Linux package so Bio-Linux and Ubuntu LTS users are encouraged to make use of the
Bio-Linux repository for these. To install Bio-Linux packages follow the instruction below.

You first need to add Bio-Linux repositories to your system. Add the following lines to a file "/etc/apt/sources.list" file:

```sh
deb http://nebc.nerc.ac.uk/bio-linux/ unstable bio-linux
deb http://ppa.launchpad.net/nebc/bio-linux/ubuntu trusty main
deb-src http://ppa.launchpad.net/nebc/bio-linux/ubuntu trusty main
```

Then run the command:

```sh
sudo apt-get update
sudo apt-get install bio-linux-keyring
```

Your system will then download the Bio-Linux Software Package List from the Bio-Linux server. After this you may install any of the packages available in Bio-Linux repository.

So, let's start installing the dependencies:

```sh
sudo apt-get install python-biom-format vsearch fastx-toolkit hmmer
```


1.9 (Misc) How to uninstall PIPITS
----------------------------------

You can uninstall PIPITS simply by deleting $HOME/pipits directory.


2. Getting started
==================

The PIPITS pipeline is divided into three parts:

1.  PIPITS_PREP: prepares raw reads from Illumina MiSeq sequencers for
    ITS extraction
2.  PIPITS_FUNITS: extracts fungal ITS regions from the reads
3.  PIPITS_PROCESS: analyses the reads to produce Operational Taxonomic
    Unit (OTU) abundance tables and the RDP taxonomic assignment table
    for downstream analysis

2.1 PIPITS_PREP
---------------

Illumina reads are generally provided as demultiplexed FASTQ files where
the Illumina machine software splits the reads into separate files, one
for each barcode.

PIPITS provides a script called PIPITS_GETREADPAIRSLIST which generates
a tab-delimited text file for all read-pairs from the raw sequence
directory:

```sh
pipits_getreadpairslist -i illumina_rawdata -o readpairslist.txt
```

*Note*

1.  The command produces a tab-delimited file with three columns
    denoting forward and reverse read filenames and sample IDs for the
    pairs
2.  Prior to running the command, the user needs to ensure that the raw
    sequence filenames end with one of the following extensions:
    “.fastq”, “.fastq.bz2”, “.fastq.gz”. Sample IDs are taken from the
    first characters preceding an underscore (“_”) from the filenames
3.  Before proceeding to the next step, check correct filenames and
    desired sample IDs for the pairs are listed in the resulting file
    ("readpairslist.txt")

Once we have the list file, we can begin to process the sequences:

```sh
pipits_prep -i illumina_rawdata_directory -o pipits_prep -l readpairslist.txt
```

*Note*

1.  Read-pairs are joined by examining the overlapping regions of
    sequences with PEAR
2.  The resulting assembled reads are then quality filtered with
    FASTX_FASTQ_QUALITY_FILTER
3.  The header of each read is then relabelled with an index number followed by 
    a sample ID
4.  The resulting files are converted into a FASTA format with
    FASTX_FASTQ_TO_FASTA and merged into a single file to produce the
    final output file "prepped.fasta"" in the output directory
5.  PIPITS_PREP can be run without the list file provided that the files
    in the input directory follow illumina file naming convention. It is
    generally recommended however to first make a list file of
    read-pairs prior to running the script


2.2 PIPITS FUNITS
-----------------

The output from PIPITS PREP is taken as an input for this step. It is
also mandatory to provide the script with which ITS subregion (i.e. ITS1
or ITS2) is to be extracted:

```sh
pipits_funits -i pipits_prep/prepped.fasta -o pipits_funits -x ITS2
```

*Note*

1.  Seqeunces are dereplicated (removing redundant sequences)
2.  Selected subregions of sequences of fungal origin are then extracted
    with ITSx and where necessary they are re-orientated to 5’ to 3’
    direction. It is worth noting that ITSx uses HMMER3 (Mistry et al.
    2013) to compare input sequences against a set of models built from
    a number of different subregions of ITS sequences found in various
    organisms. This makes ITSx an ideal tool for both extraction of
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


2.3 PIPITS PROCESS
------------------

This is the final process involving clustering and assigning of taxonomy
to OTUs:

```sh
pipits_process -i pipits_funits/ITS.fasta -o out_process --Xmx 2G 
```

*Note*

1.  Input sequences are dereplicated
2.  Short (< 100bp) and unique sequences are removed prior to finding
    OTUs
3.  The sequences are clustered at a user-defined threshold (97%
    sequence identity by default).
4.  The resulting representative sequences for each cluster are
    subjected to chimera detection and removal using the UNITE uchime
    reference dataset.
5.  The input sequences are then mapped onto the chimera-free
    representative sequences at the defined threshold
6.  The representatives are taxonomically assigned with RDP Classifier
    against the UNITE fungal ITS reference dataset.
7.  The results are then translated into two types of OTU abundance
    tables:
    1.  “OTU abundance table”, an OTU is defined as a cluster of reads
        with the user-defined threshold typically 97% sequence identity
        motivated by the expectation that these correspond approximately
        to species.
    2.  “phylotype abundance table”, an OTU is defined as a cluster of
        sequences binned into the same taxonomic assignments.
8.  If you have memory issues, try increasing the maximum memory with "--Xmx". For example, "--Xmx 4G".


4. Options
==========

PIPITS scripts come with a number of options for the users to alter
parameters such as distance threshold. The options can be viewed by
providing "-h" after the command, for example:

```sh
pipits_prep -h
```


5. Citation
===========

Hyun S. Gweon, Anna Oliver, Joanne Taylor, Tim Booth, Melanie Gibbs, Daniel S. Read, Robert I. Griffiths and Karsten Schonrogge, PIPITS: an automated pipeline for analyses of fungal internal transcribed spacer sequences from the Illumina sequencing platform, Methods in Ecology and Evolution, DOI: 10.1111/2041-210X.12399