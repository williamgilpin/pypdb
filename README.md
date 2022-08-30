# PyPDB

A Python 3 toolkit for performing searches with the RCSB Protein Data Bank (PDB). This can be used to perform advanced searches for PDB IDs matching various criteria, as well as to look up information associated with specific PDB IDs. This tool allows standard operations that can be perfomed from within the PDB website (BLAST, PFAM lookup, etc.) to be performed from within Python scripts.

Examples of each function and its associated output can be found in [`demos/demos.ipynb`](demos/demos.ipynb).

If you use this module for any published work, please consider citing the accompanying paper

      Gilpin, W. "PyPDB: A Python API for the Protein Data Bank." 
      Bioinformatics, Oxford Journals, 2016.

**Help wanted: Currently, pypdb is in maintainence-only mode. We would very much appreciate help with adding features and developing a long-term roadmap.**

## Installation

Install using pip:

	$ pip install pypdb

To install the development version, which contains the latest features and fixes, install directly from GitHub using

   	$ pip install git+git://github.com/williamgilpin/pypdb

If you need to  install directly from setup.py,

    $ python setup.py install

Test the installation, and check that the code successfully connects to the PDB, navigate to the root directory and run

	$ pytest 

This code has been designed and tested for Python 3.

## Usage

This package can be used to get lists of PDB IDs associated with specific search terms, experiment types, structures, and other common criteria.

Given a list of PDBs, this package can be used to fetch any data associated with those PDBs, including their dates of deposition, lists of authors and associated publications, their sequences or structures, their top BLAST matches, and other query-specific attributes like lists of a ligands or chemical structure.

A set of demos is included in the iPython notebook **demos.ipynb**. A static version of this notebook (for viewing) is available as **demos.html**

## Issues and Feature Requests

If you run into an issue, or if you find a workaround for an existing issue, we would very much appreciate it if you could post your question or code as a GitHub issue.

If posting a feature request, please check that your request is possible using [the current GUI on current RCSB website](https://www.rcsb.org/search/advanced). If so, please perform your search, and then click the link that says `JSON` in the upper right hand corner of the Advanced Search box. Please post that JSON code with your feature request.




