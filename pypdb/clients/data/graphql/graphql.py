"""Contains logic to perform arbitrary GraphQL searches against RCSB.

For the differences between the GraphQL and RESTful searches, see:
https://data.rcsb.org/index.html#gql-vs-rest
"""

from typing import Any  # DO NOT APPROVE: fix this to actual type

RSCB_GRAPHQL_URL = "https://data.rcsb.org/graphql?query="

def search_graphql(graphql_json_query : str) -> Any:
    """Performs RCSB search with JSON query using GraphQL.

    For details on what the RCSB GraphQL interface is, see:
        https://data.rcsb.org/index.html#gql-api

    This function should return the equivalent information as this site:
        https://data.rcsb.org/graphql/index.html

    Args:
        graphql_json_query: GraphQL JSON query, as a string. Whitespace doesn't
            matter. e.g. "{entry(entry_id:"4HHB"){exptl{method}}}"
    """

    # Strips all whitespace from JSON string, for URL encode
    stripped_json_query = "".join(graphql_json_query.split())
    URL_REQUEST = f"{RSCB_GRAPHQL_URL}{stripped_json_query}"
    print(URL_REQUEST)
