from distutils.core import setup
import os

__version__ = os.environ.get("VERSION", "1.0.0")

setup(name = "pipits",
      version = __version__,
      description = "PIPITS: An automated pipeline for analyses of fungal ITS sequences from the Illumina sequencing platform",
      author = "Hyun Soon Gweon",
      author_email = "hyugwe@ceh.ac.uk",
      url = "http://sourceforge.net/projects/pipits",
      packages = ["pipits", "pipits.pipits_uc"],
      scripts = ["bin/pipits", "bin/retrain_rdp"],
      long_description = "An open source stand-alone suite of software for automated processing of Illumina MiSeq sequences for fungal community analysis. PIPITS exploits a number of state of the art applications including manipulating paired end reads; automated ITS subregion filtering; OTU picking; and notably is the first pipeline to employ the sensitive RDP Classifier to taxonomically assign sequences against the taxonomically curated fungal ITS UNITE database,. We provide detailed descriptions of the pipeline and show its utility in the analysis of fungal ITS sequences generated on the MiSeq platform.  from Illumina MiSeq. "
      )
