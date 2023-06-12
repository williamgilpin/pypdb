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

from pypdb.clients.search.operators import sequence_operators
from pypdb.clients.search.operators import text_operators
from pypdb.clients.search.operators.chemical_operators import ChemicalOperator
from pypdb.clients.search.operators.seqmotif_operators import SeqMotifOperator
from pypdb.clients.search.operators.sequence_operators import SequenceOperator
from pypdb.clients.search.operators.structure_operators import StructureOperator
from pypdb.clients.search.operators.text_operators import TextSearchOperator

SEARCH_URL_ENDPOINT: str = "https://search.rcsb.org/rcsbsearch/v2/query"
"""SearchOperators correspond to individual search operations.

These can be used to search on their own using `perform_search`, or they can be
aggregated together into a `QueryGroup` to search using multiple operators at
once using `perform_search_with_graph`.
"""
SearchOperator = Union[TextSearchOperator, SequenceOperator, StructureOperator,
                       SeqMotifOperator]


class LogicalOperator(Enum):
    """Operation used to combine `QueryGroup` results."""
    AND = "and"
    OR = "or"


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
    # Elements within the list of `queries` can either be `SearchOperator`
    # instances (corresponding to individual queries)
    # or `QueryGroup` instances (corresponding to groups of queries).
    #
    # This allows building arbitrarily complex query logic in the search tree.
    queries: List[Union[SearchOperator, "QueryGroup"]]

    # Boolean to aggregate the results of `queries`.
    logical_operator: LogicalOperator

    def _to_dict(self):
        return {
            "type":
            "group",
            "logical_operator":
            self.logical_operator.value,
            "nodes": [
                _QueryNode(query)._to_dict()
                if type(query) is not QueryGroup else query._to_dict()
                for query in self.queries
            ]
        }


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
            result_dict["paginate"] = {
                "start": self.result_start_index,
                "rows": self.num_results
            }

        if self.sort_by != None and self.desc != None:
            result_dict["sort"] = [{
                "sort_by": self.sort_by,
                "direction": "desc" if self.desc else "asc"
            }]

        return result_dict


@dataclass
class ScoredResult:
    entity_id: str  # PDB Entity ID (e.g. 5JUP for the entry return type)
    score: float


RawJSONDictResponse = Dict[str, Any]


def perform_search(
    search_operator: SearchOperator,
    return_type: ReturnType = ReturnType.ENTRY,
    request_options: Optional[RequestOptions] = None,
    return_with_scores: bool = False,
    return_raw_json_dict: bool = False,
    verbosity: bool = True,
) -> Union[List[str], List[ScoredResult], RawJSONDictResponse]:
    """Performs search specified by `search_operator`.
    Returns entity strings of type `return_type` that match the resulting hits.

    Strictly a subset of the functionality exposed in
    `perform_search_with_graph`, this function does not support searching on
    multiple conditions at once.

    If you require this functionality, please use `perform_search_with_graph`
    instead.

    Args:
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
        verbosity: Print out the search query to the console (default: True)

    Returns:
        List of entity ids, corresponding to entities that match the given
        query.

        If `return_with_scores=True`, returns a list of ScoredResult instead.
        If `return_raw_json_dict=True`, returns the raw JSON response from RCSB.

    Example usage to search for PDB entries that are from 'Mus musculus':
    ```
    from pypdb.clients.search. import perform_search
    from pypdb.clients.search. import ReturnType
    from pypdb.clients.search.operators.text_operators import ExactMatchOperator
    pdb_ids = perform_search(
               search_operator=text_operators.ExactMatchOperator(
                 attribute="rcsb_entity_source_organism.taxonomy_lineage.name",
                 value="Mus musculus"
               ),
               return_type=ReturnType.ENTRY)
    print(pdb_ids)
    )
    ```
    """

    return perform_search_with_graph(query_object=search_operator,
                                     return_type=return_type,
                                     request_options=request_options,
                                     return_with_scores=return_with_scores,
                                     return_raw_json_dict=return_raw_json_dict,
                                     verbosity=verbosity)


_SEARCH_OPERATORS = text_operators.TEXT_SEARCH_OPERATORS + [
    SequenceOperator, StructureOperator, SeqMotifOperator
]


def perform_search_with_graph(
    query_object: Union[SearchOperator, QueryGroup],
    return_type: ReturnType = ReturnType.ENTRY,
    request_options: Optional[RequestOptions] = None,
    return_with_scores: bool = False,
    return_raw_json_dict: bool = False,
    verbosity: bool = True,
) -> Union[List[str], RawJSONDictResponse, List[ScoredResult]]:
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
        query_object: Fully-specified SearchOperator or QueryGroup
            object corresponding to the desired search.
        return_type: Type of entities to return.
        return_with_scores: Whether or not to return the entity results with
            their associated scores. For example, you might want to do this to
            get the top X hits that are similar to a certain protein sequence.
        return_raw_json_dict: Whether to return raw JSON response.
            (for example, to analyze the scores of various matches)
        verbosity: Print out the search query to the console (default: True)

    Returns:
        List of strings, corresponding to hits in the database. Will be of the
        format specified by the `return_type`.

        If `return_with_scores=True`, returns a list of ScoredResult instead.
        If `return_raw_json_dict=True`, returns the raw JSON response from RCSB.
    """

    if type(query_object) in _SEARCH_OPERATORS:
        cast_query_object = _QueryNode(query_object)  # type: ignore
    else:
        cast_query_object = query_object  # type: ignore

    if request_options is not None:
        request_options_dict = request_options._to_dict()
    else:
        request_options_dict = {'return_all_hits': True}

    rcsb_query_dict = {
        "query": cast_query_object._to_dict(),
        "request_options": request_options_dict,
        "return_type": return_type.value
    }

    if verbosity:
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
            results.append(
                ScoredResult(entity_id=query_hit["identifier"],
                             score=query_hit["score"]))
        else:
            results.append(query_hit["identifier"])

    return results


class SearchService(Enum):
    """Which type of field is being searched.

    Auto-inferred from search operator."""
    BASIC_SEARCH = "full_text"
    TEXT = "text"
    SEQUENCE = "sequence"
    SEQMOTIF = "seqmotif"
    STRUCTURE = "structure"
    CHEMICAL = "chemical"


class CannotInferSearchServiceException(Exception):
    """Raised when the RCSB Search API Service cannot be inferred."""


def _infer_search_service(search_operator: SearchOperator) -> SearchService:

    if isinstance(search_operator, text_operators.DefaultOperator):
        return SearchService.BASIC_SEARCH
    elif type(search_operator) in text_operators.TEXT_SEARCH_OPERATORS:
        return SearchService.TEXT
    elif type(search_operator) is SequenceOperator:
        return SearchService.SEQUENCE
    elif type(search_operator) is StructureOperator:
        return SearchService.STRUCTURE
    elif type(search_operator) is SeqMotifOperator:
        return SearchService.SEQMOTIF
    elif type(search_operator) is ChemicalOperator:
        return SearchService.CHEMICAL
    else:
        raise CannotInferSearchServiceException(
            "Cannot infer Search Service for {}".format(type(search_operator)))


@dataclass
class _QueryNode:
    """Individual query node, performing a query defined by the provided
    `search_operator`
    """
    search_operator: SearchOperator

    def _to_dict(self):
        return {
            "type": "terminal",
            "service": _infer_search_service(self.search_operator).value,
            "parameters": self.search_operator._to_dict()
        }
