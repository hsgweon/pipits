from setuptools import setup, find_packages
import os

__version__ = os.environ.get("VERSION", "4.0")

setup(
    name="pipits",
    version=__version__,
    packages=['pipits', 'pipits.bin'],
    # scripts=[], # Keep scripts empty if using entry_points
    description="PIPITS: An automated pipeline for fungal ITS analyses",
    long_description="An automated pipeline for analyses of fungal internal transcribed spacer (ITS) sequences from the Illumina sequencing platform.", # Added a placeholder description
    author="Hyun Soon Gweon",
    author_email="h.s.gweon@reading.ac.uk",
    url="https://github.com/hsgweon/pipits",
    license="GPL-3.0-or-later",
    entry_points={
        'console_scripts': [
            # Corrected entry points to point to the main function
            'pipits_createreadpairslist = pipits.bin.pipits_createreadpairslist:main',
            'pipits_prep = pipits.bin.pipits_prep:main',
            'pipits_funits = pipits.bin.pipits_funits:main',
            'pipits_process = pipits.bin.pipits_process:main',
        ],
    },
    package_data={'pipits.bin': []}, # Usually empty if only using scripts/entry points
    classifiers=[ # Optional: Add classifiers for better package metadata
        'Development Status :: 4 - Beta', # Or 5 - Production/Stable if appropriate
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10', # Specify supported version
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.10', # Specify minimum Python version
    # Add dependencies here if needed, e.g.:
    # install_requires=[
    #     'biom-format',
    #     # Add other direct dependencies if PIPITS code imports them directly
    # ],
)
