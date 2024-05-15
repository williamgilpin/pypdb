from setuptools import setup

# read the contents of the README file so that PyPI can use it as the long description                               
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

modules_list = [
    "pypdb",
    "pypdb.util",
    "pypdb.clients",
    "pypdb.clients.search",
    "pypdb.clients.search.operators",
    "pypdb.clients.data",
    "pypdb.clients.data.graphql",
    "pypdb.clients.fasta",
    "pypdb.clients.pdb",
]

setup(
    name='pypdb',
    packages=modules_list,  # same as 'name'
    py_modules=modules_list,
    version='2.04',
    install_requires=['requests'],
    description='A Python wrapper for the RCSB Protein Data Bank (PDB) API',
    author='William Gilpin',
    author_email='firstnamelastname@gmail.com',
    url='https://github.com/williamgilpin/pypdb',
    download_url='https://github.com/williamgilpin/pypdb/tarball/0.6',
    keywords=['protein', 'data', 'RESTful', 'api'],
    classifiers=[],
    long_description=long_description,
    long_description_content_type='text/markdown'
)
