"""Resolves the installed package version.

The version is derived from the git tag by setuptools-scm, which writes it to
`pypdb/_version.py` at build time. That file is generated rather than checked
in, so this module falls back to the installed package metadata for source
checkouts that have not been built.
"""

__all__ = ["get_version"]


def get_version() -> str:
    """Returns the version of the installed `pypdb` package."""
    try:
        from pypdb._version import __version__

        return __version__
    except ImportError:
        pass

    try:
        from importlib.metadata import PackageNotFoundError, version

        return version("pypdb")
    except PackageNotFoundError:
        # Running from an unbuilt, uninstalled source checkout
        return "0.0.0.dev0"
