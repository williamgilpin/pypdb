"""Tests for RCSB Text Search Service Operators
(admittedly, a lot is tested in `search_client_test.py` too)
"""

import pytest
import unittest

from pypdb.clients.search.operators import sequence_operators


class TestSequenceOperators(unittest.TestCase):
    def test_sequence_operator(self):
        search_operator = sequence_operators.SequenceOperator(
            sequence="AUGAUUCGGCGCUAAAAAAAA",
            sequence_type=sequence_operators.SequenceType.RNA,
            evalue_cutoff=100,
            identity_cutoff=0.95)

        self.assertEqual(
            search_operator._to_dict(), {
                "evalue_cutoff": 100,
                "identity_cutoff": 0.95,
                "target": "pdb_rna_sequence",
                "value": "AUGAUUCGGCGCUAAAAAAAA",
            })

    def test_autoresolve_sequence_type(self):
        self.assertEqual(
            sequence_operators.SequenceOperator("ATGGGGTAA").sequence_type,
            sequence_operators.SequenceType.DNA)
        self.assertEqual(
            sequence_operators.SequenceOperator("AUGGGGCCCUAA").sequence_type,
            sequence_operators.SequenceType.RNA)
        self.assertEqual(
            sequence_operators.SequenceOperator(
                "MAETREGGQSGAAS").sequence_type,
            sequence_operators.SequenceType.PROTEIN)
        with pytest.raises(
                sequence_operators.CannotAutoresolveSequenceTypeError):
            sequence_operators.SequenceOperator("AAAAAAAA")
