try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
  name = 'pypdb',
  packages = ['pypdb'], # same as 'name'
  version = '0.1',
  install_requires=[
        'xmltodict', 
        'dicttoxml',
  ],
  description = 'A Python wrapper for the Protein Data Bank API',
  author = 'William Gilpin',
  author_email = 'firstnamelastname(as one word)@googleemailservice',
  url = 'https://github.com/williamgilpin/pypdb',
  download_url = 'https://github.com/williamgilpin/pypdb/tarball/0.1', 
  keywords = ['protein','data','RESTful','api'],
  classifiers = [],
)
