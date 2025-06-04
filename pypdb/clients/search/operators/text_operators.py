"""Implementation of SearchOperators for text queries against RCSB API."""
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, Union, List

# --- Implementations of RCSB Queries for each SearchOperators ---
# See: https://search.rcsb.org/index.html#search-operators for details

# For information on available RCSB search attributes, see:
# https://search.rcsb.org/search-attributes.html


@dataclass
class DefaultOperator:
    """Default search operator; searches across available fields search,
    and returns a hit if a match happens in any field."""
    value: str

    def _to_dict(self) -> Dict[str, str]:
        return {"value": self.value}


@dataclass
class ExactMatchOperator:
    """Exact match operator indicates that the input value should match a field
    value exactly (including whitespaces, special characters and case)."""
    attribute: str
    value: Any

    def _to_dict(self) -> Dict[str, Any]:
        return {
            "attribute": self.attribute,
            "operator": "exact_match",
            "value": self.value
        }


@dataclass
class InOperator:
    """The in operator allows you to specify multiple values in a single search
    expression. It returns results if any value in a list of input values
    matches. It can be used instead of multiple OR conditions."""
    attribute: str
    values: List[Any]  # List of strings, numbers or date strings

    def _to_dict(self) -> Dict[str, Any]:
        return {
            "attribute": self.attribute,
            "operator": "in",
            "value": self.values
        }


@dataclass
class ContainsWordsOperator:
    """Searches attribute field to check if any words within `value` are found.

    For example, "actin-binding protein" will return results containing
    "actin" OR "binding" OR "protein" within the attribute.
    """
    attribute: str
    value: str

    def _to_dict(self) -> Dict[str, str]:
        return {
            "attribute": self.attribute,
            "operator": "contains_words",
            "value": self.value
        }


@dataclass
class ContainsPhraseOperator:
    """Searches attribute, and returns hits if-and-only-if all words in the
    value are in the attribute field, in that order.

    For example, "actin-binding protein" will be interpreted as
    "actin" AND "binding" AND "protein" occurring in a given order."""
    attribute: str
    value: str

    def _to_dict(self) -> Dict[str, str]:
        return {
            "attribute": self.attribute,
            "operator": "contains_phrase",
            "value": self.value
        }


class ComparisonType(Enum):
    GREATER = "greater"
    GREATER_OR_EQUAL = "greater_or_equal"
    EQUAL = "equals"
    NOT_EQUAL = "not_equal"
    LESS_OR_EQUAL = "less_or_equal"
    LESS = "less"


# TODO(lacoperon): Add support for initializing this, and RangeOperator, from
#                  datetime.datetime objects for ease of use.


@dataclass
class ComparisonOperator:
    """Searches attribute, returns hits if the attribute field comparison to the
    value is True.

    For example, to get structures after a certain date, you could use the
    following:

    ```
    date_filter_operator = ComparisonOperator(
          value="2019-01-01T00:00:00Z",
          attribute="rcsb_accession_info.initial_release_date",
          comparison_type=ComparisonType.GREATER)
    ```
    """

    attribute: str
    value: Any
    comparison_type: ComparisonType

    def _to_dict(self) -> Dict[str, Any]:
        if self.comparison_type is ComparisonType.NOT_EQUAL:
            param_dict = {"operator": "equals", "negation": True}
        else:
            param_dict = {"operator": self.comparison_type.value}

        param_dict["attribute"] = self.attribute
        param_dict["value"] = self.value

        return param_dict


# @dataclass
# class RangeOperator:
#     """Returns results with attributes within range."""
#     attribute: str
#     from_value: Any
#     to_value: Any
#     include_lower: bool = True  # Default inclusive
#     include_upper: bool = True  # Default inclusive

#     def _to_dict(self) -> Dict[str, Any]:
#         return {
#             "operator": "range",
#             "attribute": self.attribute,
#             "value": {
#                 "from": self.from_value,
#                 "to": self.to_value,
#                 "include_lower": self.include_lower,
#                 "include_upper": self.include_upper
#             }
#         }


@dataclass
class RangeOperator:
    """Returns results with attributes within range.."""
    attribute: str
    from_value: Any
    to_value: Any
    include_lower: bool = True  # Default inclusive
    include_upper: bool = True  # Default inclusive
    negation: bool = False

    def _to_dict(self) -> Dict[str, Any]:
        return {
            "operator": "range",
            "attribute": self.attribute,
            "negation": self.negation,
            "value": {"from": self.from_value,
                      "to": self.to_value},
        }


@dataclass
class ExistsOperator:
    attribute: str

    def _to_dict(self) -> Dict[str, str]:
        return {"operator": "exists", "attribute": self.attribute}


# An object of type `TextSearchOperator` can be any of the following classes:
TextSearchOperator = Union[DefaultOperator, ExactMatchOperator, InOperator,
                           ContainsWordsOperator, ContainsPhraseOperator,
                           ComparisonOperator, RangeOperator, ExistsOperator]

# List of all TextSearchOperator-associated classes, for backwards compatability
# in terms of checking SearchOperator validity
# (please change this when you change the `Union` definition)
TEXT_SEARCH_OPERATORS = [
    DefaultOperator, ExactMatchOperator, InOperator, ContainsWordsOperator,
    ContainsPhraseOperator, ComparisonOperator, RangeOperator, ExistsOperator
]
