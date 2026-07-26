# PyPDB

A Python 3 toolkit for performing searches with the RCSB Protein Data Bank (PDB). This can be used to perform advanced searches for PDB IDs matching various criteria, as well as to look up information associated with specific PDB IDs. This tool allows standard operations that can be perfomed from within the PDB website (BLAST, PFAM lookup, etc.) to be performed from within Python scripts.

If you use this module for any published work, please consider citing the accompanying paper

      Gilpin, W. "PyPDB: A Python API for the Protein Data Bank."
      Bioinformatics, Oxford Journals, 2016.

## Installation

Install using pip:

    $ pip install pypdb

Or using conda, from [conda-forge](https://anaconda.org/conda-forge/pypdb):

    $ conda install -c conda-forge pypdb

To install the development version, which contains the latest features and fixes, install directly from GitHub using

    $ pip install git+https://github.com/williamgilpin/pypdb

To install from a local checkout,

    $ pip install .

For development and test dependencies,

    $ pip install -e ".[dev]"

Test the installation, and check that the code successfully connects to the PDB, navigate to the root directory and run

    $ python -m pytest

This code has been designed and tested for Python 3.

## Usage

### PDB Text Search

Search the PDB for entries matching a term, and get back a list of PDB IDs:

```python
from pypdb import Query

found_pdbs = Query("ribosome").search()
print(found_pdbs[:5])
```

The `query_type` argument searches a specific field instead of the full text. A few of the most common:

```python
# By source organism
Query("Dictyostelium", query_type="OrganismQuery").search()

# By Enzyme Classification number (partial numbers match the whole subtree)
Query("2.3.2.5", query_type="ec_number").search()

# By UniProt accession
Query("P68871", query_type="uniprot").search()

# By experimental method
Query("SOLID-STATE NMR", query_type="ExpTypeQuery").search()

# By ligand, either as a chemical component ID or a SMILES string
Query("NAG", query_type="chemical").search()
Query("Clc1nc(Br)nc2nc[nH]c12", query_type="chemical").search()
```

Other query types include `PubmedIdQuery`, `TreeEntityQuery` (NCBI TaxID), `AdvancedAuthorQuery`, `pfam`, `sequence`, `seqmotif`, and `structure`. See `help(Query)` for the full list.

For queries combining several conditions with AND/OR logic, negation, ranges, or comparisons, see [`search/EXAMPLES.md`](pypdb/clients/search/EXAMPLES.md).

### PDB Data Fetch

Look up the information associated with a single PDB ID:

```python
from pypdb import get_info

info = get_info("4HHB")
print(info["struct"]["title"])     # THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN...
print(info["exptl"][0]["method"])  # X-RAY DIFFRACTION
```

Fetch the ligands bound to an entry, or the description of a chemical component:

```python
from pypdb import get_ligands, describe_chemical

ligands = get_ligands("4HHB")["ligandInfo"]["ligand"]
print([ligand["@chemicalID"] for ligand in ligands])  # ['HEM', 'PO4']

print(describe_chemical("NAG")["chem_comp"]["name"])
```

Download a structure file:

```python
from pypdb import get_pdb_file

cif = get_pdb_file("4lza", filetype="cif", compression=True)
```

### Ligand coordinates

Ligand coordinates can be downloaded in either of the two forms the RCSB website offers, both in Kekule form. *Ideal* coordinates come from the Chemical Component Dictionary and are independent of any structure; *instance* coordinates are the ligand as modelled inside one particular entry:

```python
from pypdb import get_ideal_ligand_file, get_ligand_instances, get_ligand_instance_file

# Idealized coordinates for a chemical component (SDF only)
sdf = get_ideal_ligand_file("ATP")

# Coordinates of each copy as modelled in an entry (SDF or MOL2)
for instance in get_ligand_instances("4HHB"):
    coords = get_ligand_instance_file("4HHB",
                                      instance["auth_asym_id"],
                                      instance["auth_seq_id"],
                                      filetype="mol2")
```

`get_ligand_instances` reports the chain and residue number of every ligand copy in an entry, which is what identifies a specific instance to download.

Integrated PubMed, UniProt, and structural interface data can be fetched as well. For fetching many properties across many entries at once, see [`data/EXAMPLES.md`](pypdb/clients/data/EXAMPLES.md).

### Counting and paginating results

Searches return every matching entry by default, which can be slow for broad queries. To retrieve just a tally or the highest-scoring hits:

```python
from pypdb import count_results, get_top_results

print(count_results("ribosome"))            # 9560
print(get_top_results("crispr", max_results=5))
```

Searches can also be restricted to a page of results, scored with a particular strategy, or extended to include computed structure models (such as AlphaFold predictions) via `RequestOptions`. See [`search/EXAMPLES.md`](pypdb/clients/search/EXAMPLES.md).

More worked examples of every function live in [`demos/demos.ipynb`](demos/demos.ipynb).

## Releasing a new version

The git tag is the only place the version is written. To publish a release to PyPI, tag and push:

```bash
git tag v2.9 && git push origin v2.9
```

That's the whole process — there is no version file to edit. [setuptools-scm](https://setuptools-scm.readthedocs.io/) derives the package version from the tag at build time and writes it into `pypdb/_version.py`, which is generated rather than checked in (it is listed in `.gitignore`).

Pushing a `v*` tag runs the tests, builds the distributions, and uploads them to PyPI. This uses [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/), so no API token is stored in the repository — it requires a one-time publisher registration on PyPI pointing at the `python-publish.yml` workflow and the `pypi` environment.

Between releases, `pypdb.__version__` reports a development version derived from the most recent tag (e.g. `2.9.dev1+g1a2b3c4`).

## Issues and Feature Requests

If you run into an issue, or if you find a workaround for an existing issue, please post your question or code as a GitHub issue.

If posting a feature request, please check that your request is possible using [the current GUI on current RCSB website](https://www.rcsb.org/search/advanced). If so, please perform your search, and then click the link that says `JSON` in the upper right hand corner of the Advanced Search box. Please post that JSON code with your feature request.



