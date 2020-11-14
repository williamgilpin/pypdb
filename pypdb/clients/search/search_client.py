"""Barebones Python API Wrapper implementation around RCSB Search API.

This file contains Python dataclasses that formalize the API within native
Python objects. This aims to completely implement the Search API in Python.

For RCSB API docs, see: https://search.rcsb.org/index.html
"""

# TODO(lacoperon): Implement request options

from dataclasses import dataclass
from enum import Enum
import json
import requests
from typing import Any, Dict, List, Optional, Union
import warnings

from pypdb.clients.search.operators import text_operators
from pypdb.clients.search.operators.text_operators import TextSearchOperator

SEARCH_URL_ENDPOINT : str = "https://search.rcsb.org/rcsbsearch/v1/query"

class SearchService(Enum):
    """Which type of field is being searched."""
    TEXT = "text"
    SEQUENCE = "sequence"
    SEQMOTIF = "seqmotif"
    STRUCTURE = "structure"
    CHEMICAL = "chemical"


class LogicalOperator(Enum):
    """Operation used to combine `QueryGroup` results."""
    AND = "and"
    OR = "or"

# SearchOperators correspond to individual search operations, that can be
# aggregated using a `QueryGroup`.
#
# Currently, the only available search operators are associated with the
# 'text' service. For the list of available RCSB services,
# see: https://search.rcsb.org/index.html#search-services
SearchOperator = TextSearchOperator


@dataclass
class QueryNode:
    """Individual query node, associated with a search using `search_service`
    using logic defined in the `search_operator`.
    """
    search_service: SearchService
    search_operator: SearchOperator

    def to_dict(self):
        return {
            "type": "terminal",
            "service": self.search_service.value,
            "parameters": self.search_operator.to_dict()
        }


@dataclass
class QueryGroup:
    """Group of search operators against RCSB Search API,
    whose independent results are aggregated with `logical_operator`.

    For example, for searches with `query_nodes=[n1,n2,n3]`,
    and `logical_operator=LogicalOperator.AND`, results will only be
    returned for hits that match all of n1, n2 and n3's queries.

    `logical_operator=LogicalOperator.OR` would return results that match any
    of n1, n2 or n3's queries.
    """
    # Elements within the list of `queries` can either be `QueryNode` instances
    # (corresponding to individual queries)
    # or `QueryGroup` instances (corresponding to groups of queries).
    #
    # This allows building arbitrarily complex query logic in the search tree.
    queries: List[Union[QueryNode, "QueryGroup"]]

    # Boolean to aggregate the results of `queries`.
    logical_operator: LogicalOperator

    def to_dict(self):
        return {
            "type": "group",
            "logical_operator": self.logical_operator.value,
            "nodes": [query.to_dict() for query in self.queries]
        }


class ReturnType(Enum):
    """For details, see: https://search.rcsb.org/index.html#return-type"""
    ENTRY = "entry"
    ASSEMBLY = "assembly"
    POLYMER_ENTITY = "polymer_entity"
    NON_POLYMER_ENTITY = "non_polymer_entity"
    POLYMER_INSTANCE = "polymer_instance"


class InappropriateSearchOperatorException(Exception):
    """Raised when the provided SearchService and SearchOperator are
    mutually incompatible.

    For example, you can't search against the
    SEQMOTIF service using the RangeOperator, as that's not supported by the
    RCSB Search API."""

RawJSONDictResponse = Dict[str, Any]

def perform_search(search_service: SearchService,
                   search_operator: SearchOperator,
                   return_type: ReturnType,
                   return_raw_json_dict: bool=False
                   ) -> Union[List[str],
                              RawJSONDictResponse]:
    """Performs search specified by `search_operator`, against `search_service`.
    Returns entity strings of type `return_type` that match the resulting hits.

    Strictly a subset of the functionality exposed in
    `perform_search_with_graph`, this function does not support searching on
    multiple conditions at once.

    If you require this functionality, please use `perform_search_with_graph`
    instead.

    Args:
        search_service: What type of RCSB Search Service to query.
        search_operator: Parameters defining the search condition.
        return_type: What type of RCSB entity to return.
        return_raw_json_dict: If True, this function returns the raw JSON
            response from RCSB, instead of a

    Returns:
        List of entity ids, corresponding to entities that match the given
        query.

    Example usage to search for PDB entries that are from 'Mus musculus':
    ```
    from pypdb.clients.search. import perform_search
    from pypdb.clients.search. import SearchService, ReturnType
    from pypdb.clients.search.operators.text_operators import ExactMatchOperator
    pdb_ids = perform_search(
               search_service=SearchService.TEXT,
               search_operator=text_operators.ExactMatchOperator(
                 attribute="rcsb_entity_source_organism.taxonomy_lineage.name",
                 value="Mus musculus"
               ),
               return_type=ReturnType.ENTRY)
    print(pdb_ids)
    )
    ```
    """

    if search_service != SearchService.TEXT:
        raise NotImplementedError(
            "Not currently implemented (but watch this space)")

    if search_service is SearchService.TEXT:
        # Uses undocumented `__args__` to check if operator is TEXT
        # (see: https://stackoverflow.com/questions/45957615/ )
        if not type(search_operator) in text_operators.TEXT_SEARCH_OPERATORS:
            raise InappropriateSearchOperatorException(
                "Searches against the 'text' service should have SearchOperators "
                "that are TextSearchOperators, as in `text_search_operators.py`.")

    rcsb_query_dict = {
        "query": QueryNode(search_service=search_service,
                           search_operator=search_operator).to_dict(),
        "request_options": {"return_all_hits": True},
        "return_type": return_type.value
    }

    print("Querying RCSB Search using the following parameters:\n %s" %
    json.dumps(rcsb_query_dict))

    response = requests.post(url=SEARCH_URL_ENDPOINT,
                             data=json.dumps(rcsb_query_dict))

    if not response.ok:
        warnings.warn("It appears request failed with:" + response.text)
        response.raise_for_status()

    # If specified, returns raw JSON response from RCSB as Dict
    # (rather than entity IDs as a string list)
    if return_raw_json_dict:
        return response.json()

    # Converts RCSB result to list of identifiers corresponding to
    # the `return_type`.
    identifiers = []
    for query_hit in response.json()["result_set"]:
        identifiers.append(query_hit["identifier"])

    return identifiers

def perform_search_with_graph(query_object) -> List[str]:
    """Performs specified search using RCSB's search node logic.

    See https://search.rcsb.org/index.html#building-search-request under
    "Terminal node" and "Group node" for details.

    Args:
        query_object: Fully-specified QueryNode or QueryGroup
            object corresponding to the desired search.

    Returns:
        List of strings, corresponding to hits in the database. Will be of the
        format specified by the `return_type`.
    """
    raise NotImplementedError("Not currently implemented.")

    # TODO(lacoperon): Implement this.
