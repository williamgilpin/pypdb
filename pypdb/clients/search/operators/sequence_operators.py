"""Search operator for searching sequences using MMseqs2 (BLAST-like)."""
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Optional, Union


class SequenceType(Enum):
    """Type of sequence being searched."""
    DNA = "pdb_dna_sequence"
    RNA = "pdb_rna_sequence"
    PROTEIN = "pdb_protein_sequence"


class CannotAutoresolveSequenceTypeError(Exception):
    """Raised when a sequence is ambiguous as to its `SequenceType`."""


@dataclass
class SequenceOperator:
    """Default search operator; searches across available fields search,
    and returns a hit if a match happens in any field."""
    sequence: str
    # If the sequence type is not specified, tries to autoresolve the type from
    # the sequence itself
    sequence_type: Optional[SequenceType] = None
    # Maximum E Value allowed for results
    # (see: https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html)
    evalue_cutoff: float = 100
    # Minimum identity cutoff allowed for results
    # (see: https://www.ncbi.nlm.nih.gov/books/NBK62051/def-item/identity/)
    identity_cutoff: float = 0.95

    def __post_init__(self):
        if self.sequence_type is None:
            self._autoresolve_sequence_type()

    def _autoresolve_sequence_type(self):
        unique_letters = set(list(self.sequence))

        dna_letter_set = set(["A", "T", "C", "G"])
        rna_letter_set = set(["A", "U", "C", "G"])
        protein_letter_set = set(list("ABCDEFGHIKLMNPQRSTVWXYZ"))
        protein_fingerprint_set = set(list("BDEFHIKLMNPQRSVWXYZ"))
        if unique_letters.issubset(dna_letter_set) and "T" in unique_letters:
            self.sequence_type = SequenceType.DNA
        elif unique_letters.issubset(rna_letter_set) and "U" in unique_letters:
            self.sequence_type = SequenceType.RNA
        elif (unique_letters.issubset(protein_letter_set)
              and protein_fingerprint_set & unique_letters):
            self.sequence_type = SequenceType.PROTEIN
        else:
            raise CannotAutoresolveSequenceTypeError(
                "Sequence is ambiguous as to its SequenceType: `{}`".format(
                    self.sequence))

    def _to_dict(self) -> Dict[str, Any]:
        return {
            "evalue_cutoff": self.evalue_cutoff,
            "identity_cutoff": self.identity_cutoff,
            "target": self.sequence_type.value,  # type: ignore
            "value": self.sequence
        }
