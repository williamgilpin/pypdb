"""Operators for searching molecular definitions (chemical components / BIRD).

These operate against the `text_chem` search service, which searches the
textual annotations of chemical definitions rather than of PDB structures.

For details, see: https://search.rcsb.org/index.html#search-services
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Union

from pypdb.clients.search.operators.text_operators import _maybe_add_negation


@dataclass
class ChemicalExactMatchOperator:
    """Matches a molecular-definition attribute exactly.

    For example, to find the chemical component with ID "NAG":

    ```
    ChemicalExactMatchOperator(
        attribute="rcsb_chem_comp_container_identifiers.comp_id",
        value="NAG")
    ```
    """
    attribute: str
    value: Any
    negation: bool = False

    def _to_dict(self) -> Dict[str, Any]:
        return _maybe_add_negation(
            {
                "attribute": self.attribute,
                "operator": "exact_match",
                "value": self.value
            }, self.negation)


@dataclass
class ChemicalInOperator:
    """Matches a molecular-definition attribute against a list of values."""
    attribute: str
    values: List[Any]
    negation: bool = False

    def _to_dict(self) -> Dict[str, Any]:
        return _maybe_add_negation(
            {
                "attribute": self.attribute,
                "operator": "in",
                "value": self.values
            }, self.negation)


@dataclass
class ChemicalContainsWordsOperator:
    """Searches a molecular-definition attribute for any of `value`'s words."""
    attribute: str
    value: str
    negation: bool = False

    def _to_dict(self) -> Dict[str, str]:
        return _maybe_add_negation(
            {
                "attribute": self.attribute,
                "operator": "contains_words",
                "value": self.value
            }, self.negation)


# An object of type `ChemicalTextSearchOperator` can be any of the following:
ChemicalTextSearchOperator = Union[ChemicalExactMatchOperator,
                                   ChemicalInOperator,
                                   ChemicalContainsWordsOperator]

# List of all molecular-definition text operator classes, used to infer the
# `text_chem` search service (please update alongside the tuple above).
CHEMICAL_TEXT_SEARCH_OPERATORS = [
    ChemicalExactMatchOperator, ChemicalInOperator,
    ChemicalContainsWordsOperator
]
