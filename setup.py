from setuptools import setup
import os

__version__ = os.environ.get("VERSION", "2.0")

setup(
	name = "pipits",
	version = __version__,
	packages = ["pipits"],
	scripts = [	"bin/pipits_getsamplelistfromfasta",
				"bin/pipits_funits", 
				"bin/pipits_process", 
				"bin/pipits_uc2otutable",
				"bin/pipits_reformatAssignedTaxonomy", 
				"bin/pipits_phylotype_biom",
				"bin/pipits_funguild.py",
				"bin/pipits_retrain_rdp"],
	description = "PIPITS: An automated pipeline for analyses of fungal ITS sequences from the Illumina sequencing platform",
	long_description = "An open source stand-alone suite of software for automated processing of Illumina MiSeq sequences for fungal community analysis. PIPITS exploits a number of state of the art applications including manipulating paired end reads; automated ITS subregion filtering; OTU picking; and notably is the first pipeline to employ the sensitive RDP Classifier to taxonomically assign sequences against the taxonomically curated fungal ITS UNITE database,. We provide detailed descriptions of the pipeline and show its utility in the analysis of fungal ITS sequences generated on the MiSeq platform.",
	author = "Hyun Soon Gweon",
	author_email = "h.s.gweon@reading.ac.uk",
	url = "https://github.com/hsgweon/pipits",
	license = "GPA"
)
