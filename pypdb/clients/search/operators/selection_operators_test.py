"""Tests for RCSB molecular-definition (text_chem) Search Service Operators."""

import unittest

from pypdb.clients.search.operators import selection_operators


class TestSelectionOperators(unittest.TestCase):
    def test_chemical_exact_match_operator_to_dict(self):
        operator = selection_operators.ChemicalExactMatchOperator(
            attribute="rcsb_chem_comp_container_identifiers.comp_id",
            value="NAG")

        self.assertEqual(
            operator._to_dict(), {
                "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                "operator": "exact_match",
                "value": "NAG"
            })

    def test_chemical_exact_match_operator_with_negation(self):
        operator = selection_operators.ChemicalExactMatchOperator(
            attribute="rcsb_chem_comp_container_identifiers.comp_id",
            value="NAG",
            negation=True)

        self.assertEqual(
            operator._to_dict(), {
                "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                "operator": "exact_match",
                "value": "NAG",
                "negation": True
            })

    def test_chemical_in_operator_to_dict(self):
        operator = selection_operators.ChemicalInOperator(
            attribute="rcsb_chem_comp_container_identifiers.comp_id",
            values=["NAG", "GAL"])

        self.assertEqual(
            operator._to_dict(), {
                "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                "operator": "in",
                "value": ["NAG", "GAL"]
            })

    def test_chemical_contains_words_operator_to_dict(self):
        operator = selection_operators.ChemicalContainsWordsOperator(
            attribute="chem_comp.name", value="glucose")

        self.assertEqual(
            operator._to_dict(), {
                "attribute": "chem_comp.name",
                "operator": "contains_words",
                "value": "glucose"
            })


if __name__ == '__main__':
    unittest.main()
