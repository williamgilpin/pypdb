"""Operators associated with RCSB structural search."""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict


class StructureSearchMode(Enum):
    """Mode to search structures with. See:
    https://github.com/biocryst/biozernike/
    """
    STRICT_SHAPE_MATCH = "strict_shape_match"
    RELAXED_SHAPE_MATCH = "relaxed_shape_match"


@dataclass
class StructureOperator:
    """Operator to perform 3D Structural search using:
    https://github.com/biocryst/biozernike/

    Will return similar 3D structures using default search options.
    """
    # Entry and Assembly # for the chainstructure you want to use for search.
    # (results will show other PDB entities with similiar 3D Structures)
    pdb_entry_id: str
    assembly_id: int = 1
    # Structure search mode
    search_mode: StructureSearchMode = StructureSearchMode.STRICT_SHAPE_MATCH

    def _to_dict(self) -> Dict[str, Any]:
        return {
            "value": {
                "entry_id": self.pdb_entry_id,
                "assembly_id": str(self.assembly_id)
            },
            "operator": self.search_mode.value
        }
