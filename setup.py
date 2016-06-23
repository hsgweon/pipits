from distutils.core import setup
import os

__version__ = os.environ.get("VERSION", "1.3.4")

setup(name = "pipits",
      version = __version__,
      description = "PIPITS: An automated pipeline for analyses of fungal ITS sequences from the Illumina sequencing platform",
      author = "Hyun Soon Gweon",
      author_email = "hyugwe@ceh.ac.uk",
      url = "http://sourceforge.net/projects/pipits",
      packages = ["pipits"],
      scripts = ["bin/pipits_getreadpairslist",
                 "bin/pipits_getsamplelistfromfasta",
                 "bin/pipits_prep", 
                 "bin/pipits_funits", 
                 "bin/pipits_dereplicate_fasta", 
                 "bin/pipits_fasta_filter_by_length", 
                 "bin/pipits_inflate_fasta", 
                 "bin/pipits_process", 
                 "bin/pipits_uc2otutable",
                 "bin/pipits_reformatAssignedTaxonomy", 
                 "bin/pipits_phylotype_biom"],
      long_description = "An open source stand-alone suite of software for automated processing of Illumina MiSeq sequences for fungal community analysis. PIPITS exploits a number of state of the art applications including manipulating paired end reads; automated ITS subregion filtering; OTU picking; and notably is the first pipeline to employ the sensitive RDP Classifier to taxonomically assign sequences against the taxonomically curated fungal ITS UNITE database,. We provide detailed descriptions of the pipeline and show its utility in the analysis of fungal ITS sequences generated on the MiSeq platform."
      )
