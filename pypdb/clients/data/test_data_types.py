"""
Unit tests for the data_types classes.
"""
import unittest
from unittest import mock
import requests
from pypdb.clients.data.graphql.graphql import RSCB_GRAPHQL_URL

from pypdb.clients.data.data_types import Entry

class TestEntry(unittest.TestCase):
    def test_create(self):
        entry = Entry("4HHB")

        self.assertIsNone(entry.properties)
        self.assertIsNone(entry.json_query)
        self.assertIsNone(entry.response)

        self.assertEqual(entry.id, "4HHB")

        entry.check_pdb_id()

        self.assertEqual(entry.id, ["4HHB"])

    def test_generate_json_query(self):
        entry = Entry("4HHB")

        property = {"exptl":["method", "details"]}

        entry.add_property(property)

        self.assertIsNotNone(entry.properties)

        entry.generate_json_query()

        self.assertTrue(isinstance(entry.json_query, dict))
        self.assertTrue("query" in entry.json_query)

    @mock.patch.object(requests, "post")
    def test_fetch(self, mock_post):
        entry = Entry("4HHB")
        property = {"exptl": ["method"]}
        entry.add_property(property)

        expected_return = {'data': {'entry': {'exptl': {"method": "X-RAY DIFFRACTION"}}}}
        mock_response = mock.create_autospec(requests.Response, instance=True)
        mock_response.json.return_value = expected_return
        mock_post.return_value = mock_response

        entry.fetch_data()

        mock_post.assert_called_once_with(url=RSCB_GRAPHQL_URL, json=entry.json_query)
        self.assertEqual(entry.response, expected_return)

if __name__ == '__main__':
    unittest.main()