"""Contains logic to perform arbitrary GraphQL searches against RCSB.

For the differences between the GraphQL and RESTful searches, see:
https://data.rcsb.org/index.html#gql-vs-rest
"""
import requests
import  warnings
from typing import Any  # DO NOT APPROVE: fix this to actual type

RSCB_GRAPHQL_URL = "https://data.rcsb.org/graphql?query="


def search_graphql(graphql_json_query: str) -> Any:
    """Performs RCSB search with JSON query using GraphQL.

    For details on what the RCSB GraphQL interface is, see:
        https://data.rcsb.org/index.html#gql-api

    This function should return the equivalent information as this site:
        https://data.rcsb.org/graphql/index.html

    Args:
        graphql_json_query: GraphQL JSON query, as a string. Whitespace doesn't
            matter. e.g. "{entry(entry_id:"4HHB"){exptl{method}}}"
    """

    response = requests.post(url=RSCB_GRAPHQL_URL,
                             json=graphql_json_query)

    if not response.ok:
        warnings.warn(f"It appears request failed with: {response.text}")
        response.raise_for_status()


    return response.json()