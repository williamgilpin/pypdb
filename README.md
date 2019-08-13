# PyPDB

A Python 3 toolkit for performing searches with the RCSB Protein Data Bank (PDB) using its XML-based API. This can be used to perform advanced searches for PDB IDs matching various criteria, as well as to look up information associated with specific PDB IDs. This tool allows standard operations that can be perfomed from within the PDB website (BLAST, PFAM lookup, etc.) to be performed within Python scripts, allowing it to supplement existing tools (i.e. Biopython) that are designed for manipulating .pdb files.

Examples of each function and its associated output can be found in `demos/demos.nb`. Full documentation for each function can also be accessed via `help(function_name)`
Slightly out-of-date documentation may also be viewed in the separate repository [pypdb_docs](http://www.wgilpin.com/pypdb_docs/html/)

If you use this module for any published work, please consider citing the accompanying paper

      Gilpin, W. "PyPDB: A Python API for the Protein Data Bank." 
      Bioinformatics, Oxford Journals, 2015.


## Installation

Install using pip:

	$ pip install pypdb

or, using conda:

	$ conda install -c williamgilpin pypdb

To install the development version, which contains the latest bug fixes, install directly from GitHub using

   	$ pip install git+git://github.com/williamgilpin/pypdb

If you need to  install directly from setup.py,

    $ python setup.py install

This code is designed for Python 3. The installation adds the libraries [**xmltodict**](https://github.com/martinblech/xmltodict) and [**BeautifulSoup**](http://www.crummy.com/software/BeautifulSoup/). Please follow the licenses associated with these libraries.
BeautifulSoup is only necessary to use one function, get_blast2.

## Usage

This package can be used to get lists of PDB IDs associated with specific search terms, experiment types, structures, and other common criteria.

Given a list of PDBs, this package can be used to fetch any data associated with those PDBs, including their dates of deposition, lists of authors and associated publications, their sequences or structures, their top BLAST matches, and other query-specific attributes like lists of a ligands or chemical structure.

A set of demos is included in the iPython notebook **demos.ipynb**. A static version of this notebook (for viewing) is available as **demos.html**

## Development

This package has not been refined for every use case (there are many elaborate request types allowed through the PDB's xml-based API), but I believe that enough of the basic structure is present that any additional features may be easily implemented. I encourage forks and pull requests.

More information about the PDB's API and request format can be found at this page:

[The RCSB PDB RESTful Web Service interface](http://www.rcsb.org/pdb/software/rest.do)


