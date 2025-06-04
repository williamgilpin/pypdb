"""Operators associated with SeqMotif searching using RCSB Search API."""

from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict


class SequenceType(Enum):
    """Type of sequence being searched for motifs."""
    DNA = "pdb_dna_sequence"
    RNA = "pdb_rna_sequence"
    PROTEIN = "pdb_protein_sequence"


class PatternType(Enum):
    """Type of pattern being used for SeqMotif search."""
    SIMPLE = "simple"
    PROSITE = "prosite"
    REGEX = "regex"


@dataclass
class SeqMotifOperator:
    # Pattern to search with
    pattern: str
    sequence_type: SequenceType
    pattern_type: PatternType

    def _to_dict(self) -> Dict[str, Any]:
        return {
            "value": self.pattern,
            "pattern_type": self.pattern_type.value,
            "target": self.sequence_type.value
        }


# DO NOT APPROVE: DO NOT APPROVE THIS CL UNTIL ADDED TO VALIDATION
