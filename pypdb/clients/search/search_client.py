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
from pypdb.clients.search.operators import sequence_operators
from pypdb.clients.search.operators.text_operators import TextSearchOperator
from pypdb.clients.search.operators.sequence_operators import SequenceSearchOperator

SEARCH_URL_ENDPOINT: str = "https://search.rcsb.org/rcsbsearch/v1/query"


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
SearchOperator = Union[
    TextSearchOperator,
    SequenceSearchOperator]


@dataclass
class QueryNode:
    """Individual query node, associated with a search using `search_service`
    using logic defined in the `search_operator`.
    """
    search_service: SearchService
    search_operator: SearchOperator

    def _to_dict(self):
        return {
            "type": "terminal",
            "service": self.search_service.value,
            "parameters": self.search_operator._to_dict()
        }

    def _validate(self) -> None:
        """Validates queries to SearchService use a supporting SearchOperator.

        Used to raise Errors notifying users of invalid RCSB queries before
        those queries hit RCSB's Search servers."""

        if self.search_service not in [SearchService.TEXT,
                                       SearchService.SEQUENCE]:
            raise NotImplementedError(
                "This service isn't yet implemented in the RCSB 2.0 API "
                "(but watch this space)")

        # Each SearchService is assocaited with a list of valid search operators
        if self.search_service is SearchService.TEXT:
            appropriate_operator_list = text_operators.TEXT_SEARCH_OPERATORS # type: ignore
            operator_file="pypdb/clients/search/operators/text_operators.py"
        elif self.search_service is SearchService.SEQUENCE:
            appropriate_operator_list = sequence_operators.SEQUENCE_SEARCH_OPERATORS # type: ignore
            operator_file="pypdb/clients/search/operators/sequence_operators.py"
        else:
            # Default to search being OK if there's no validation for
            # this operator defined yet
            return

        if not type(self.search_operator) in appropriate_operator_list:
            raise InappropriateSearchOperatorException(
                "Searches against the '{}' service should only use ".format(self.search_service),
                " SearchOperators defined in in `{}`.".format(operator_file))

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

    def _to_dict(self):
        return {
            "type": "group",
            "logical_operator": self.logical_operator.value,
            "nodes": [query._to_dict() for query in self.queries]
        }

    def _validate(self):
        """Validates nodes within the QueryGroup contain valid queries

        Used to raise Errors notifying users of invalid RCSB queries before
        those queries hit RCSB's Search servers."""
        for query in self.queries:
            query._validate()


class ReturnType(Enum):
    """For details, see: https://search.rcsb.org/index.html#return-type"""
    ENTRY = "entry"
    ASSEMBLY = "assembly"
    POLYMER_ENTITY = "polymer_entity"
    NON_POLYMER_ENTITY = "non_polymer_entity"
    POLYMER_INSTANCE = "polymer_instance"

@dataclass
class RequestOptions:
    """Options to configure which results are returned, and in what order."""
    # Returns `num_results` results starting at`result_start_index` (pagination)
    # If these indices are not defined, defaults to return all results.
    # (returning all results can be slow for compute-intensive searches)
    result_start_index: Optional[int] = None
    num_results: Optional[int] = None
    # What attribute to sort by.
    # This should either be "score"  (to sort by score),
    # or a valid RCSB attribute value
    # (e.g. "rcsb_accession_info.initial_release_date")
    sort_by: Optional[str] = "score"
    # Whether to sort by score ascending, or score descending
    desc: Optional[bool] = True

    def _to_dict(self):
        result_dict = {}
        if self.result_start_index != None and self.num_results != None:
            result_dict["pager"] = {
                "start": self.result_start_index,
                "rows": self.num_results
            }

        if self.sort_by != None and self.desc != None:
            result_dict["sort"] = [
                {
                "sort_by": self.sort_by,
                "direction": "desc" if self.desc else "asc"
                }
            ]

        return result_dict

@dataclass
class ScoredResult:
    entity_id: str  # PDB Entity ID (e.g. 5JUP for the entry return type)
    score: float

class InappropriateSearchOperatorException(Exception):
    """Raised when the provided SearchService and SearchOperator are
    mutually incompatible.

    For example, you can't search against the
    SEQMOTIF service using the RangeOperator, as that's not supported by the
    RCSB Search API."""


RawJSONDictResponse = Dict[str, Any]


def perform_search(search_service: SearchService,
                   search_operator: SearchOperator,
                   return_type: ReturnType = ReturnType.ENTRY,
                   request_options: Optional[RequestOptions] = None,
                   return_with_scores: bool = False,
                   return_raw_json_dict: bool = False
                   ) -> Union[List[str],
                              List[ScoredResult],
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
        request_options: Object containing information for result pagination
          and sorting functionality.
        return_with_scores: Whether or not to return the entity results with
            their associated scores. For example, you might want to do this to
            get
            the top X hits that are similar to a certain protein sequence.
            (if this is true, returns List[ScoredResult] instead of List[str])
        return_raw_json_dict: If True, this function returns the raw JSON
            response from RCSB, instead of a

    Returns:
        List of entity ids, corresponding to entities that match the given
        query.

        If `return_with_scores=True`, returns a list of ScoredResult instead.
        If `return_raw_json_dict=True`, returns the raw JSON response from RCSB.

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

    query_node = QueryNode(search_service=search_service,
                           search_operator=search_operator)

    return perform_search_with_graph(query_object=query_node,
                                     return_type=return_type,
                                     request_options=request_options,
                                     return_with_scores=return_with_scores,
                                     return_raw_json_dict=return_raw_json_dict)


def perform_search_with_graph(query_object: Union[QueryNode, QueryGroup],
                              return_type: ReturnType = ReturnType.ENTRY,
                              request_options: Optional[RequestOptions] = None,
                              return_with_scores: bool = False,
                              return_raw_json_dict: bool = False
                              ) -> Union[List[str],
                                         RawJSONDictResponse,
                                         List[ScoredResult]]:
    """Performs specified search using RCSB's search node logic.

    Essentially, this allows you to ask multiple questions in one RCSB query.

    For example, you can ask for structures that satisfy all of the following
    conditions at once:
        * Are either from Mus musculus or from Homo sapiens lineage
        * Are both under 4 angstroms of resolution, and published after 2019
        * Are labelled as "actin-binding protein" OR
            contain "actin" AND "calmodulin" in their titles.

    See https://search.rcsb.org/index.html#building-search-request under
    "Terminal node" and "Group node" for more details.

    Args:
        query_object: Fully-specified QueryNode or QueryGroup
            object corresponding to the desired search.
        return_type: Type of entities to return.
        return_with_scores: Whether or not to return the entity results with
            their associated scores. For example, you might want to do this to
            get the top X hits that are similar to a certain protein sequence.
        return_raw_json_dict: Whether to return raw JSON response.
            (for example, to analyze the scores of various matches)

    Returns:
        List of strings, corresponding to hits in the database. Will be of the
        format specified by the `return_type`.

        If `return_with_scores=True`, returns a list of ScoredResult instead.
        If `return_raw_json_dict=True`, returns the raw JSON response from RCSB.
    """

    # Validates that, to the best of our knowledge, the `query_object`
    # is a valid query against the RCSB Search API.
    query_object._validate()

    if request_options is not None:
        request_options_dict = request_options._to_dict()
    else:
        request_options_dict = {'return_all_hits': True}

    rcsb_query_dict = {
        "query": query_object._to_dict(),
        "request_options": request_options_dict,
        "return_type": return_type.value
    }

    print("Querying RCSB Search using the following parameters:\n %s \n" %
          json.dumps(rcsb_query_dict))

    response = requests.post(url=SEARCH_URL_ENDPOINT,
                             data=json.dumps(rcsb_query_dict))


    # If your search queries are failing here, it could be that your attribute
    # doesn't support the SearchOperator you're using.
    # See: https://search.rcsb.org/search-attributes.html
    if not response.ok:
        warnings.warn("It appears request failed with:" + response.text)
        response.raise_for_status()

    # If specified, returns raw JSON response from RCSB as Dict
    # (rather than entity IDs as a string list)
    if return_raw_json_dict:
        return response.json()

    # Converts RCSB result to list of identifiers corresponding to
    # the `return_type`. Annotated with score if `return_with_scores`.
    results = []
    for query_hit in response.json()["result_set"]:
        if return_with_scores:
            results.append(ScoredResult(
                entity_id=query_hit["identifier"],
                score=query_hit["score"]
            ))
        else:
            results.append(query_hit["identifier"])

    return results
