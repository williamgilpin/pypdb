import time
import requests
import unittest
from unittest import mock
import warnings

from pypdb.clients.data.graphql import graphql


class TestGraphQLRequests(unittest.TestCase):
    def test_fails_with_invalid_request(self):
        self.assertTrue(False)

if __name__ == '__main__':
    unittest.main()
