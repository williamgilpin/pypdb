from setuptools import setup


setup(
  name = 'pypdb',
  packages = ['pypdb', 'pypdb.util', 'pypdb.clients'], # same as 'name'
  py_modules = ['pypdb', 'pypdb.util', 'pypdb.clients'],
  version = '2.0',
  install_requires=[
		'requests'
  ],
  description = 'A Python wrapper for the RCSB Protein Data Bank (PDB) API',
  author = 'William Gilpin',
  author_email = 'firstname_lastname@gmail.com',
  url = 'https://github.com/williamgilpin/pypdb',
  download_url = 'https://github.com/williamgilpin/pypdb/tarball/0.6', 
  keywords = ['protein','data','RESTful','api'],
  classifiers = [],
)
