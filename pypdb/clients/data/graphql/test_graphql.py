"""Unit tests for RCSB DATA API Python wrapper."""
import unittest
from unittest import mock
import requests

from pypdb.clients.data.graphql import graphql

class TestGraphQL(unittest.TestCase):
    @mock.patch.object(requests, "post")
    def test_simple_search(self, mock_post):
        json_query = {'query': '{ entry(entry_id: "4HHB"){struct {title}} }'}
        expected_return_json_as_dict = {'data': {'entry': {'struct': {'title': 'THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 ANGSTROMS RESOLUTION'}}}}

        mock_response = mock.create_autospec(requests.Response, instance=True)
        mock_response.json.return_value = expected_return_json_as_dict
        mock_post.return_value = mock_response

        results = graphql.search_graphql(json_query)

        mock_post.assert_called_once_with(url=graphql.RSCB_GRAPHQL_URL, json=json_query)
        self.assertEqual(results, expected_return_json_as_dict)


if __name__ == '__main__':
    unittest.main()