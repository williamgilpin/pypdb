"""Logic to download ligand coordinate files from the RCSB Database.

RCSB serves ligand coordinates in two distinct flavours, mirroring the
"Download Ideal CCD" / "Download Instance Coordinates" options on the website:

* *Ideal* coordinates come from the Chemical Component Dictionary (CCD), and
  describe an idealized, structure-independent conformer of the component.
* *Instance* coordinates describe one particular copy of the ligand as it was
  modelled inside a given PDB entry, in that entry's coordinate frame.

Both are written in Kekule form, with explicit bond orders.
"""

from enum import Enum
from typing import Optional
import warnings

from pypdb.util import http_requests

# Serves idealized coordinates from the Chemical Component Dictionary
LIGAND_DOWNLOAD_BASE_URL = "https://files.rcsb.org/ligands/download/"
# ModelServer serves coordinates of a ligand as modelled within an entry
MODEL_SERVER_BASE_URL = "https://models.rcsb.org/v1/"


class LigandFileType(Enum):
    """Format to download ligand coordinates in.

    Note that RCSB only publishes idealized CCD coordinates as SDF; MOL2 is
    available for instance coordinates.
    """
    SDF = "sdf"
    MOL2 = "mol2"


def get_ideal_ligand_file(chem_id: str,
                          filetype: LigandFileType = LigandFileType.SDF,
                          verbosity: bool = True) -> Optional[str]:
    """Fetches idealized coordinates for a chemical component from the CCD.

    These coordinates are independent of any particular structure, which makes
    them well suited to building reference libraries of small molecules.

    Args:
        chem_id: Chemical component ID of the ligand (e.g. `"ATP"`).
        filetype: Format to download. RCSB only publishes idealized
            coordinates as SDF.
        verbosity: Print out the download URL to the console (default: True)

    Returns:
        The coordinate file as a string, or None if the request failed.

    Example:
        >>> sdf = get_ideal_ligand_file("ATP")
        >>> print(sdf.splitlines()[0])
        ATP
    """
    if filetype is not LigandFileType.SDF:
        warnings.warn(
            "RCSB only publishes idealized CCD coordinates in SDF format. "
            "Use `get_ligand_instance_file` for MOL2 coordinates.")
        return None

    ligand_url = "{}{}_ideal.{}".format(LIGAND_DOWNLOAD_BASE_URL,
                                        chem_id.upper(), filetype.value)

    if verbosity:
        print("Sending GET request to {} to fetch {}'s idealized {} file.".
              format(ligand_url, chem_id, filetype.value))

    response = http_requests.request_limited(ligand_url)

    if response is None or not response.ok:
        warnings.warn("Retrieval failed, returning None")
        return None

    return response.text


def get_ligand_instance_file(pdb_id: str,
                             auth_asym_id: str,
                             auth_seq_id: int,
                             filetype: LigandFileType = LigandFileType.SDF,
                             verbosity: bool = True) -> Optional[str]:
    """Fetches the coordinates of one ligand as modelled within a PDB entry.

    Unlike the idealized CCD coordinates, these are the experimentally
    determined coordinates, in the entry's own frame of reference.

    The chain and residue number identify which copy of the ligand to fetch,
    since the same component often appears several times in one structure.
    `get_ligand_instances` lists the available copies.

    Args:
        pdb_id: A 4 character string giving a pdb entry of interest.
        auth_asym_id: Author-assigned chain ID of the ligand (e.g. `"A"`).
        auth_seq_id: Author-assigned residue number of the ligand.
        filetype: Format to download (SDF or MOL2).
        verbosity: Print out the download URL to the console (default: True)

    Returns:
        The coordinate file as a string, or None if the request failed.

    Example:
        >>> sdf = get_ligand_instance_file("4HHB", "A", 142)
        >>> print(sdf.splitlines()[0])
        HEM
    """
    ligand_url = ("{}{}/ligand?auth_asym_id={}&auth_seq_id={}&encoding={}"
                  .format(MODEL_SERVER_BASE_URL, pdb_id.lower(), auth_asym_id,
                          auth_seq_id, filetype.value))

    if verbosity:
        print("Sending GET request to {} to fetch {}'s ligand {} file.".format(
            ligand_url, pdb_id, filetype.value))

    response = http_requests.request_limited(ligand_url)

    if response is None or not response.ok:
        warnings.warn("Retrieval failed, returning None")
        return None

    # ModelServer answers a query that matches no atoms with HTTP 200 and a
    # body containing only diagnostics, so the atom count is what actually
    # distinguishes a hit from a miss.
    if "<model_server_stats.element_count>\n0" in response.text:
        warnings.warn(
            "No ligand found in {} at chain {} residue {}, returning None. "
            "Use `get_ligand_instances` to list the available copies.".format(
                pdb_id, auth_asym_id, auth_seq_id))
        return None

    return response.text
