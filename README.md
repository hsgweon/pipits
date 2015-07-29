Major update (29 July 2015)
---------------------------

Before start using PIPITS, it is important to note that PIPITS works on
Illumina sequences which have already been demultiplexed (i.e. each sample is already divided into different files)
Although it's most likely that your sequencing certre provided you with demultiplxed
FASTQ files, if this isn't the case, we recommend using deML (https://bioinf.eva.mpg.de/deml)
to demultiplex you files before using PIPITS.

Also, please note that PIPITS has gone through some substantial changes since its publication,
so the commands listed in the paper is going to be slightly different, so please refer to this documentation for detail.


1 PIPITS Setup
==============

1.1 Download and install
------------------------

Download the latest package:

    $ wget https://github.com/hsgweon/pipits/archive/master.zip

Then enter into the created directory and install the package with:

    $ cd pipits-master
    $ python setup.py install --prefix=$HOME/pipits

This creates a "pipits" directory in your $HOME and we will be installing pipits and all of its dependencies into this directory.
Of course if you are familiar with Linux, you are free to choose other ways to install.


1.2 Dependencies
----------------

*Download and install dependencies*

PIPITS depends on a number of external dependencies which need to be
downloaded and installed. (If you are a Ubuntu user, it is much much easier to make use 
of Bio-Linux packages rather than installing dependencies yourself. See 1.8 below)

Unless you are much familiar with Linux, I recommend installing dependencies this way. 
All we need to do is download, install and make them "visible" to PIPITS.
We will install every dependencies in the following directory.


-   BIOM-FORMAT v.1.3.x (<https://pypi.python.org/pypi/biom-format/1.3.1>)

    $ cd $HOME/pipits
    $ wget https://pypi.python.org/packages/source/b/biom-format/biom-format-1.3.1.tar.gz
    $ tar xfz biom-format-1.3.1.tar.gz
    $ cd biom-format-1.3.1
    $ python setup.py install --prefix=$HOME/pipits

-   FAST-X tools (<http://hannonlab.cshl.edu/fastx_toolkit>)

    $ cd $HOME/pipits
    $ wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
    $ tar xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2

-   VSEARCH (<https://github.com/torognes/vsearch>)

    $ cd $HOME/pipits
    $ wget https://github.com/torognes/vsearch/releases/download/v1.1.3/vsearch-1.1.3-linux-x86_64
    $ chmod +x vsearch-1.1.3-linux-x86_64
    $ ln -s $HOME/pipits/vsearch-1.1.3-linux-x86_64 bin/vsearch

-   ITSx (<http://microbiology.se/software/itsx>) N.B. ITSx requires
    HMMER3

    $ cd $HOME/pipits
    $ wget http://microbiology.se/sw/ITSx_1.0.11.tar.gz
    $ tar xvfz ITSx_1.0.11.tar.gz
    $ ln -s $HOME/pipits/ITSx_1.0.11/ITSx bin/ITSx
    $ ln -s $HOME/pipits/ITSx_1.0.11/ITSx_db bin/ITSx_db

-   PEAR (<http://sco.h-its.org/exelixis/web/software/pear>) - N.B. PEAR
    prohibits commercial use of the code. See its page for detail.

    $ cd $HOME/pipits
    $ wget http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.6-bin-64.tar.gz
    $ tar xvfz pear-0.9.6-bin-64.tar.gz
    $ ln -s $HOME/pipits/pear-0.9.6-bin-64/pear-0.9.6-bin-64 bin/pear

-   RDP Classifier 2.9 or above
    (<http://sourceforge.net/projects/rdp-classifier>) - N.B. RDP
    Classifier comes with a jar file.
    
    $ cd $HOME/pipits
    $ wget http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.10.2.zip
    $ unzip rdp_classifier_2.10.2.zip
    $ ln -s rdp_classifier_2.10.2/dist/classifier.jar ./classifier.jar

-   HMMER3 (<http://hmmer.janelia.org/download.html>) - This is needed for ITSx. 
    Choose "with Linux/Intel x86_64 binaries" unless you are using an "old" 32-bit PC (unlikely for most people I presume).

    $ cd $HOME/pipits
    $ wget http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz
    $ tar xfz hmmer-3.1b2-linux-intel-x86_64.tar.gz
    $ cd hmmer-3.1b2-linux-intel-x86_64
    $ ./configure --prefix $HOME/pipits
    $ make
    $ cd ..
    $ ln -s $HOME/pipits/hmmer-3.1b2-linux-intel-x86_64/binaries/* bin/


1.3 Reference datasets
----------------------

There are two reference datasets which need to be downloaded:

1. UNITE fungal ITS reference training dataset

Download the most recent version of UNITE training data from their
sourceforge webpage
(<http://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData>
), save it to an appropriate directory (e.g. $HOME/pipits/refdb) and
extract the file.

For example:

    $ mkdir -p $HOME/pipits/refdb
    $ cd $HOME/pipits/refdb
    $ wget http://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/fungalits_UNITE_trainingdata_07042014.zip
    $ unzip fungalits_UNITE_trainingdata_07042014.zip

The extracted directory contains UNITE training files namely (a) FASTA
file and (b) taxonomy file with lineage. We will use these files for
retraining RDP Classifier a moment later.


2. UNITE UCHIME reference dataset

We also need to download UNITE UCHIME reference dataset for chimera
removal. Download it from UNITE repository
(<http://unite.ut.ee/repository.php>).

For example:

    $ mkdir -p $HOME/pipits/refdb
    $ cd $HOME/pipits/refdb
    $ wget https://unite.ut.ee/sh_files/uchime_reference_dataset_26.07.2014.zip
    $ unzip uchime_reference_dataset_26.07.2014.zip


1.4 Set PATH and ENVIRONMENT VARIABLE
-------------------------------------

Make sure executables and modules are visible to the shell by existing in the search
PATH. Also set some environment variables as shown below. 
Assuming UBUNTU is your system, this can be achieved by adding the following line in "~/.zshrc" file:

    export PATH=$HOME/pipits/bin:$PATH
    export PYTHONPATH=$HOME/pipits/lib/python2.7/site-packages:$PYTHONPATH
    export PIPITS_UNITE_REFERENCE_DATA_CHIMERA=$HOME/pipits/refdb/final_release_version/uchime_sh_refs_dynamic_original_985_03.07.2014.fasta
    export PIPITS_UNITE_RETRAINED_DIR=$HOME/pipits/refdb/unite_retrained
    export PIPITS_RDP_CLASSIFIER_JAR=$HOME/pipits/classifier.jar

Then type (or alternatively close and re-open the terminal):

    $ source ~/.zshrc


1.5 Re-HMMPressing
------------------

Also once you downloaded and installed ITSx, we recommend re-HMMPRESSing the HMM profiles as the HMMPRESS'ed profiles may not be compatible with the version of the HMMER3 you installed:

    $ cd $HOME/pipits/ITSx_1.0.11/ITSx_db/HMMs
    $ rm *.hmm.*
    $ echo *.hmm | xargs -n1 hmmpress


1.6 Retrain RDP Classifier
--------------------------

Lastly we need to re-train RDP Classifier with the downloaded "UNITE fungal ITS
reference training dataset". PIPITS provides a script called
"retrain_rdp" for this task. To run the command, you need to give (i,
ii) the files from "UNITE fungal ITS reference training dataset"; (iii)
output directory name; and (iv) the location of the RDP Classifier .jar
file. Note that this step does not need to be repeated until a new set
of training data is available to retrain the classifier. For example:

    $ cd $HOME/pipits/refdb
    $ pipits_retrain_rdp -f fungalits_UNITE_trainingdata_07042014/UNITE.RDP_04.07.14.rmdup.fasta -t fungalits_UNITE_trainingdata_07042014/UNITE.RDP_04.07.14.tax -j $HOME/pipits/classifier.jar -o unite_retrained


1.7 Test Dependencies and PIPITS
--------------------------------

When you have successfully installed these, check if they are *indeed* successfully installed by running each applications. If you get an error,
check to see if you have followed the instruction carefully.

    $ biom
    $ fastq_to_fasta -h
    $ vsearch
    $ ITSx -h
    $ pear
    $ ls $HOME/pipits/classifier.jar
    $ hmmpress -h


Ok, let's test if PIPITS is all setup. Open up the very first original PIPITS which you downloaded. 

    $ cd pipits-master
    $ pipits_getreadpairslist -i test_data
    $ pipits_prep -i test_data
    $ pipits_funits -i pipits_prep/prepped.fasta -x ITS2 
    $ pipits_process -i pipits_funits/ITS.fasta -l readpairslist.txt

Ensure everything works and you don't get an error message.


1.8 (Misc) For Ubuntu users, installing dependencies using Bio-Linux packages
------------------------------------------------------------------------------------

Some of the dependencies are available as a Bio-Linux package so
Bio-Linux and Ubuntu LTS users are encouraged to make use of the
Bio-Linux repository for these. To install Bio-Linux packages follow the
instruction below.

You first need to add Bio-Linux repositories to your system. Add the
following lines to a file "/etc/apt/sources.list file":

    $ deb http://nebc.nerc.ac.uk/bio-linux/ unstable bio-linux
    $ deb http://ppa.launchpad.net/nebc/bio-linux/ubuntu precise main
    $ deb-src http://ppa.launchpad.net/nebc/bio-linux/ubuntu precise main

Then run the command:

    $ sudo apt-get update
    $ sudo apt-get install bio-linux-keyring

Your system will then download the Bio-Linux Software Package List from
our server. After the download you may install any of the packages from
our repository.

Then install the packages, for example:

    $ sudo apt-get install python-biom-format vsearch fastx-toolkit


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

    pipits_getreadpairslist -i illumina_rawdata/ -o readpairslist.txt

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

    pipits_prep -i illumina_rawdata_directory -o out_prep -l readpairslist.txt

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

    pipits_funits -i pipits_prep/prepped.fasta -o out_funits -x ITS2

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

    pipits_process -i pipits_funits/ITS.fasta -o out_process

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


4. Options
==========

PIPITS scripts come with a number of options for the users to alter
parameters such as distance threshold. The options can be viewed by
providing "-h" after the command, for example:

    $ pipits_prep -h


5. Citation
===========

Hyun S. Gweon, Anna Oliver, Joanne Taylor, Tim Booth, Melanie Gibbs, Daniel S. Read, Robert I. Griffiths and Karsten Schonrogge, PIPITS: an automated pipeline for analyses of fungal internal transcribed spacer sequences from the Illumina sequencing platform, Methods in Ecology and Evolution, DOI: 10.1111/2041-210X.12399