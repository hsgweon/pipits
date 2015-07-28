1 PIPITS Setup
==============

1.1 Download and install
------------------------

Download the latest package:

    $ git clone https://github.com/hsgweon/pipits

Then enter into the created directory and install the package with:

    $ cd pipits
    $ sudo python setup.py install

Alternatively if you want to install to a location other than the
standard location, use "prefix".

    $ python setup.py install --prefix=$HOME/.local

This is recommended if you don't have the root access. Make sure
executables are visible to the shell by existing in the search
path, by adding "$HOME/.local/bin" to your PATH variable.
Assuming you are using UBUNTU, this can be achieved by adding the
following line in "~/.zshrc" file:

    export PATH=$HOME/.local/bin:$PATH

Then type (or alternatively close and re-open the terminal. Basic Linux!):

    $ source ~/.zshrc



1.2 Dependencies
----------------

*Download and install dependencies*

PIPITS depends on a number of external dependencies which need to be
downloaded and installed.

-   BIOM-FORMAT v.1.3.x
    (<https://pypi.python.org/pypi/biom-format/1.3.1>)
-   FAST-X tools (<http://hannonlab.cshl.edu/fastx_toolkit>).
-   VSEARCH (<https://github.com/torognes/vsearch>)
-   ITSx (<http://microbiology.se/software/itsx>) N.B. ITSx requires
    HMMER3
-   PEAR (<http://sco.h-its.org/exelixis/web/software/pear>) - N.B. PEAR
    prohibits commercial use of the code. See its page for detail.
-   RDP Classifier 2.9 or above
    (<http://sourceforge.net/projects/rdp-classifier>) - N.B. RDP
    Classifier comes with a jar file.
-   HMMER3 (<http://hmmer.janelia.org/download.html>) - Choose "with Linux/Intel x86_64 binaries" unless you are using an "old" 32-bit PC (unlikely for most people I presume).

Once you downloaded and installed ITSx, we recommend you to re-HMMPRESS the HMM profiles as the HMMPRESS'ed profiles may not be compatible with the version of the HMMER you installed.
So assumming you downloaded ITSx in "$HOME/Software" directory and installed it there:

    $ cd $HOME/Software/ITSx_1.0.11/ITSx_db/HMMs
    $ rm *.hmm.*
    $ echo *.hmm | xargs -n1 hmmpress


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

    $ mkdir -p $HOME/pipits_refdb
    $ cd $HOME/pipits_refdb
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

    $ cd $HOME/pipits_refdb
    $ wget https://unite.ut.ee/sh_files/uchime_reference_dataset_26.07.2014.zip
    $ unzip uchime_reference_dataset_26.07.2014.zip


1.4 Retrain RDP Classifier
--------------------------

We need to re-train RDP Classifier with the downloaded "UNITE fungal ITS
reference training dataset". PIPITS provides a script called
"retrain_rdp" for this task. To run the command, you need to give (i,
ii) the files from "UNITE fungal ITS reference training dataset"; (iii)
output directory name; and (iv) the location of the RDP Classifier .jar
file. Note that this step does not need to be repeated until a new set
of training data is available to retrain the classifier. For example:

    $ cd $HOME/pipits_refdb
    $ retrain_rdp -f fungalits_UNITE_trainingdata_07042014/UNITE.RDP_04.07.14.rmdup.fasta -t fungalits_UNITE_trainingdata_07042014/UNITE.RDP_04.07.14.tax -j $HOME/Software/rdp_classifier_2.9/dist/classifier.jar -o unite_retrained


1.5 Configuration file
----------------------

Good! Now everything is set up, but before we get started, we need to
let PIPITS know where all these dependencies, datasets, trained datasets
are by editing a PIPITS configuration file. You can find a template
"pipits_config" file in the downloaded PIPITS directory.

It looks like:

    ########################
    # PIPITS configuration #
    ########################

    # Lines beginning with "#" is ignored.

    [DEPENDENCIES]
    # If the dependencies are not visible to the shell, then their locations need to be defined here, for example:

    FASTX_FASTQ_QUALITY_FILTER = /usr/bin/fastq_quality_filter
    FASTX_FASTQ_TO_FASTA =       /usr/bin/fastq_to_fasta
    BIOM =                       /usr/bin/biom
    PEAR =                       $HOME/Software/pear-0.9.5-bin-64/pear-0.9.5-64
    VSEARCH =                    $HOME/Software/vsearch
    ITSx =                       $HOME/Software/ITSx_1.0.10/ITSx


    # Location of RDP Classifier jar file must be specified, for example:

    RDP_CLASSIFIER_JAR =            $HOME/Software/rdp_classifier_2.9/dist/classifier.jar

    [DB]
    UNITE_REFERENCE_DATA_CHIMERA =  $HOME/pipits/refdb/final_release_version/uchime_sh_refs_dynamic_original_985_03.07.2014.fasta
    UNITE_RETRAINED_DIR =           $HOME/pipits/refdb/unite_retrained

In the configuration file:

First, the dependency executables (except RDP Classifier as this comes
with a .jar file) need to be specified in the "pipits_config" file
unless the executables are visible to the shell by existing in
executable search path (i.e. listed in the $PATH environment variable).

Second, you need to let PIPITS know where RDP Classifier is by adding
the location of RDP Classifier ".jar"" file.

Third,the location of the UNITE UCHIME reference data, and the
re-trained UNITE directory need to be specified.

Last but not least, once this file has been edited copy the file to your
home directory and rename it as ".pipits_config" by typing the
following:

    $ cd $HOME/pipits-1.0.0
    $ cp pipits_config  $HOME/.pipits_config

Ok, now you are all good to go! Let's go back to the downloaded PIPITS
directory and test to see if everything is set up correctly by running
PIPITS on a test dataset provided:

    $ cd $HOME/pipits-1.0.0
    $ pipits all -i test_data -x ITS2 --prefix pipits_test -v

Ensure you see "PIPITS_PROCESS ended successfully."

1.3 (Misc) For Ubuntu LTS users, installing dependencies using Bio-Linux packages
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

2. Getting started
==================

The PIPITS pipeline is divided into three parts:

1.  PIPITS PREP: prepares raw reads from Illumina MiSeq sequencers for
    ITS extraction
2.  PIPITS FUNITS: extracts fungal ITS regions from the reads
3.  PIPITS PROCESS: analyses the reads to produce Operational Taxonomic
    Unit (OTU) abundance tables and the RDP taxonomic assignment table
    for downstream analysis

2.1 PIPITS PREP
---------------

Illumina reads are generally provided as demultiplexed FASTQ files where
the Illumina machine software splits the reads into separate files, one
for each barcode.

PIPITS provides a script called PIPITS_GETREADPAIRSLIST which generates
a tab-delimited text file for all read-pairs from the raw sequence
directory:

    pipits getreadpairslist -i illumina_rawdata/ -o readpairslist.txt

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

    pipits prep -i illumina_rawdata -o out_prep -l readpairslist.txt

*Note*

1.  Read-pairs are joined by examining the overlapping regions of
    sequences with PEAR
2.  The resulting assembled reads are then quality filtered with
    FASTX_FASTQ_QUALITY_FILTER
3.  The header of each read is then relabelled with an index number and
    a sample ID
4.  The resulting files are converted into a FASTA format with
    FASTX_FASTQ_TO_FASTA and merged into a single file to produce the
    final output file "prepped.fasta"" in the output directory
5.  PIPITS PREP can be run without the list file provided that the files
    in the input directory follow illumina file naming convention. It is
    generally recommended however to first make a list file of
    read-pairs prior to running the script

2.2 PIPITS FUNITS
-----------------

The output from PIPITS PREP is taken as an input for this step. It is
also mandatory to provide the script with which ITS subregion (i.e. ITS1
or ITS2) is to be extracted:

    pipits funits -i out_prep/prepped.fasta -o out_funits/ -x ITS2

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

    pipits process -i out_funits/ITS.fasta -o out_process

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

3. All in one!
==============

Provided the user is aware of the three steps above, it is possible to
run the entire PIPITS pipeline with a single command (which was used for
testing the whole pipeline above):

    pipits all -i illumina_rawdata -x ITS2 --prefix mypipits

4. Options
==========

PIPITS scripts come with a number of options for the users to alter
parameters such as distance threshold. The options can be viewed by
providing "-h" after the command, for example:

    pipits prep -h
