"""Operators associated with RCSB structure motif search.

A structure motif search finds spatial arrangements of residues, such as a
catalytic triad, regardless of the sequence or fold they sit in. The motif is
defined by pointing at a handful of residues within an existing PDB entry.

For details, see: https://search.rcsb.org/index.html#search-services
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

# RCSB requires a motif to be defined by between 2 and 10 residues
MIN_MOTIF_RESIDUES = 2
MAX_MOTIF_RESIDUES = 10


class AtomPairingScheme(Enum):
    """Which atoms are used when superimposing a hit onto the motif."""
    ALL = "ALL"
    BACKBONE = "BACKBONE"
    SIDE_CHAIN = "SIDE_CHAIN"
    PSEUDO_ATOMS = "PSEUDO_ATOMS"


@dataclass
class StructureMotifResidue:
    """One residue of a structure motif.

    Residues are addressed by their mmCIF `label_asym_id` and `label_seq_id`,
    which are not always the same as the author-assigned chain and residue
    numbers shown on the RCSB website.
    """
    label_asym_id: str
    label_seq_id: int
    # Residue types permitted at this position in addition to the original one
    # (e.g. `["LYS", "HIS"]` to also allow those residues to match here)
    exchanges: Optional[List[str]] = None

    def _to_dict(self) -> Dict[str, Any]:
        return {
            "label_asym_id": self.label_asym_id,
            "label_seq_id": self.label_seq_id
        }


class InvalidMotifResidueCountError(Exception):
    """Raised when a motif has too few or too many residues."""


@dataclass
class StrucMotifOperator:
    """Operator to perform a 3D structure motif search.

    The motif is defined by residues within `pdb_entry_id`. For example, to
    search for the enolase superfamily's catalytic residues:

    ```
    StrucMotifOperator(
        pdb_entry_id="2MNR",
        residues=[
            StructureMotifResidue(label_asym_id="A", label_seq_id=162),
            StructureMotifResidue(label_asym_id="A", label_seq_id=193),
            StructureMotifResidue(label_asym_id="A", label_seq_id=219),
        ])
    ```
    """
    # Entry the motif is defined within
    pdb_entry_id: str
    # Between 2 and 10 residues making up the motif
    residues: List[StructureMotifResidue] = field(default_factory=list)
    # Maximum allowed RMSD (in Angstroms) between a hit and the motif
    rmsd_cutoff: Optional[float] = None
    # Which atoms to use when superimposing a hit onto the motif
    atom_pairing_scheme: Optional[AtomPairingScheme] = None

    def __post_init__(self):
        if not (MIN_MOTIF_RESIDUES <= len(self.residues) <=
                MAX_MOTIF_RESIDUES):
            raise InvalidMotifResidueCountError(
                "A structure motif must be defined by between {} and {} "
                "residues, but {} were given.".format(MIN_MOTIF_RESIDUES,
                                                      MAX_MOTIF_RESIDUES,
                                                      len(self.residues)))

    def _to_dict(self) -> Dict[str, Any]:
        value: Dict[str, Any] = {
            "entry_id": self.pdb_entry_id,
            "residue_ids": [residue._to_dict() for residue in self.residues]
        }

        params: Dict[str, Any] = {"value": value}

        if self.rmsd_cutoff is not None:
            params["rmsd_cutoff"] = self.rmsd_cutoff

        if self.atom_pairing_scheme is not None:
            params["atom_pairing_scheme"] = self.atom_pairing_scheme.value

        exchanges = [{
            "residue_id": residue._to_dict(),
            "allowed": residue.exchanges
        } for residue in self.residues if residue.exchanges]
        if exchanges:
            params["exchanges"] = exchanges

        return params
