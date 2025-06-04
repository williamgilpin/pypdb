"""Tests for Structural searches against RCSB Search API."""

import unittest

from pypdb.clients.search.operators import structure_operators


class TestStructureOperators(unittest.TestCase):
    def test_not_equals_operator(self):
        structure_operator = structure_operators.StructureOperator(
            pdb_entry_id="HK97",
            assembly_id=4,
            search_mode=structure_operators.StructureSearchMode.
            STRICT_SHAPE_MATCH)

        self.assertEqual(
            structure_operator._to_dict(), {
                "value": {
                    "entry_id": "HK97",
                    "assembly_id": "4"
                },
                "operator": "strict_shape_match"
            })

        structure_operator_two = structure_operators.StructureOperator(
            pdb_entry_id="CP77",
            assembly_id=7,
            search_mode=structure_operators.StructureSearchMode.
            RELAXED_SHAPE_MATCH)

        self.assertEqual(
            structure_operator_two._to_dict(), {
                "value": {
                    "entry_id": "CP77",
                    "assembly_id": "7"
                },
                "operator": "relaxed_shape_match"
            })
