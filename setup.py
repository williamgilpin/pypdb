from setuptools import setup


setuptools.setup(
  name = 'pypdb',
  packages = ['pypdb','beautifulsoup'], # same as 'name'
  version = '0.2',
  install_requires=[
        'xmltodict', 
  ],
  description = 'A Python wrapper for the Protein Data Bank (PDB) API',
  author = 'William Gilpin',
  author_email = 'firstnamelastname(as one word)@googleemailservice',
  url = 'https://github.com/williamgilpin/pypdb',
  download_url = 'https://github.com/williamgilpin/pypdb/tarball/0.2', 
  keywords = ['protein','data','RESTful','api'],
  classifiers = [],
)
