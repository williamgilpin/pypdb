"""Tests for RCSB Text Search Service Operators
(admittedly, a lot is tested in `search_client_test.py` too)
"""

import unittest

from pypdb.clients.search.operators import text_operators


class TestTextOperators(unittest.TestCase):
    def test_not_equals_operator(self):
        not_equals_operator = text_operators.ComparisonOperator(
            attribute="struct.favourite_marvel_movie",
            value="Thor: Ragnarok",
            comparison_type=text_operators.ComparisonType.NOT_EQUAL)

        self.assertEqual(
            not_equals_operator._to_dict(), {
                "attribute": "struct.favourite_marvel_movie",
                "value": "Thor: Ragnarok",
                "operator": "equals",
                "negation": True
            })
