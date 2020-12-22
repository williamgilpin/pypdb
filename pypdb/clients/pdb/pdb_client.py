"""File containing logic to download PDB file entries from the RCSB Database."""

from enum import Enum
import gzip
from typing import Optional
import warnings

from pypdb.util import http_requests

PDB_DOWNLOAD_BASE_URL = "https://files.rcsb.org/download/"


class PDBFileType(Enum):
    PDB = "pdb"  # Older file format.
    CIF = "cif"  # Newer file format (replacing PDB file type)
    XML = "xml"  # Another alternative representation.
    STRUCTFACT = "structfact"  # For structural factors (only populated for some entries)


def get_pdb_file(pdb_id: str,
                 filetype=PDBFileType.PDB,
                 compression=False) -> Optional[str]:
    '''Get the full PDB file associated with a PDB_ID

    Parameters
    ----------

    pdb_id : A 4 character string giving a pdb entry of interest

    filetype: The file type.
        PDB is the older file format,
        CIF is the newer replacement.
        XML an also be obtained and parsed using the various xml tools included in PyPDB
        STRUCTFACT retrieves structure factors (only available for certain PDB entries)

    compression : Whether or not to request the data as a compressed (gz) version of the file
        (note that the compression is undone by this function)

    Returns
    -------

    result : string
        The string representing the full PDB file as an uncompressed string.
        (returns None if the request to RCSB failed)

    Examples
    --------
    >>> pdb_file = get_pdb_file('4lza', filetype='cif', compression=True)
    >>> print(pdb_file[:200])
    data_4LZA
    #
    _entry.id   4LZA
    #
    _audit_conform.dict_name       mmcif_pdbx.dic
    _audit_conform.dict_version    4.032
    _audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx

    '''

    if filetype is PDBFileType.CIF and not compression:
        warnings.warn("Consider using `get_pdb_file` with compression=True "
                      "for CIF files (it makes the file download faster!)")

    pdb_url_builder = [PDB_DOWNLOAD_BASE_URL, pdb_id]

    if filetype is PDBFileType.STRUCTFACT:
        pdb_url_builder.append("-sf.cif")
    else:
        pdb_url_builder += [".", filetype.value]

    if compression:
        pdb_url_builder += ".gz"

    pdb_url = "".join(pdb_url_builder)

    print(
        "Sending GET request to {} to fetch {}'s {} file as a string.".format(
            pdb_url, pdb_id, filetype.value))

    response = http_requests.request_limited(pdb_url)

    if response is None or not response.ok:
        warnings.warn("Retrieval failed, returning None")
        return None

    if compression:
        return gzip.decompress(response.content)
    return response.text
