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

    def test_negated_contains_words_operator(self):
        negated_operator = text_operators.ContainsWordsOperator(
            attribute="struct.title",
            value="kinase inhibitor",
            negation=True)

        self.assertEqual(
            negated_operator._to_dict(), {
                "attribute": "struct.title",
                "operator": "contains_words",
                "value": "kinase inhibitor",
                "negation": True
            })

    def test_negated_default_operator(self):
        negated_operator = text_operators.DefaultOperator(
            value="ribosome",
            negation=True)

        self.assertEqual(
            negated_operator._to_dict(), {
                "value": "ribosome",
                "negation": True
            })

    def test_negated_exists_operator(self):
        negated_operator = text_operators.ExistsOperator(
            attribute="rcsb_primary_citation.pdbx_database_id_PubMed",
            negation=True)

        self.assertEqual(
            negated_operator._to_dict(), {
                "operator": "exists",
                "attribute": "rcsb_primary_citation.pdbx_database_id_PubMed",
                "negation": True
            })
