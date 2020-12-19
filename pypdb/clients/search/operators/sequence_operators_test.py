"""Tests for RCSB Text Search Service Operators
(admittedly, a lot is tested in `search_client_test.py` too)
"""

import unittest

from pypdb.clients.search.operators import sequence_operators

class TestSequenceOperators(unittest.TestCase):

    def test_sequence_operator(self):
        search_operator = sequence_operators.SequenceOperator(
            sequence="AUGAUUCGGCGCUAAAAAAAA",
            sequence_type = sequence_operators.SequenceType.RNA,
            evalue_cutoff=100,
            identity_cutoff=0.95
        )

        self.assertEqual(
            search_operator._to_dict(),
            {
                "evalue_cutoff": 100,
                "identity_cutoff": 0.95,
                "target": "pdb_rna_sequence",
                "value": "AUGAUUCGGCGCUAAAAAAAA",
            }
        )
