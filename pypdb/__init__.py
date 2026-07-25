from ._version_helper import get_version

# Derived from the git tag at build time; see pyproject.toml
__version__ = get_version()

from .pypdb import *
# from .pypdb.util import *
