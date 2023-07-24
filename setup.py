from setuptools import setup

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
    version='2.03',
    install_requires=['requests', 'pandas'],
    description='A Python wrapper for the RCSB Protein Data Bank (PDB) API',
    author='William Gilpin',
    author_email='firstnamelastname@gmail.com',
    url='https://github.com/williamgilpin/pypdb',
    download_url='https://github.com/williamgilpin/pypdb/tarball/0.6',
    keywords=['protein', 'data', 'RESTful', 'api'],
    classifiers=[],
)
