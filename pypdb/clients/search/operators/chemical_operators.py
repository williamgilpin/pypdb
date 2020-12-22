"""Search operators corresponding to Chemical search using SMILES or InChI."""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict


class DescriptorMatchingCriterion(Enum):
    """Criterion describing what constitutes a chemical 'match' in RCSB search.

    For definitions of these criteria, see:
    https://search.rcsb.org/#search-services
    """
    GRAPH_STRICT = "graph-strict"
    GRAPH_RELAXED = "graph-relaxed"
    GRAPH_RELAXED_STEREO = "graph-relaxed-stereo"
    FINGERPRINT_SIMILARITY = "fingerprint-similarity"


@dataclass
class ChemicalOperator:
    """Search operator for Chemical searches using SMILES / InChI."""
    # Descriptor for matching (i.e. a valid SMILES or InChI string)
    descriptor: str
    # Criterion for what constitutes a match ("graph-strict" by default)
    matching_criterion: DescriptorMatchingCriterion = DescriptorMatchingCriterion.GRAPH_STRICT

    def __post_init__(self):
        """Derives whether the chemical descriptor string is SMILES or InChI."""
        # All InChI strings definitionally start with "InChI="
        if self.descriptor.startswith("InChI="):
            self.descriptor_type = "InChI"
        else:
            # Otherwise, assume SMILES string by default
            self.descriptor_type = "SMILES"

    def _to_dict(self) -> Dict[str, Any]:
        return {
            "value": self.descriptor,
            "type": "descriptor",
            "descriptor_type": self.descriptor_type,
            "match_type": self.matching_criterion.value
        }
