"""Tests for RCSB SeqMotif Search Service Operators."""

import unittest

from pypdb.clients.search.operators import seqmotif_operators


class TestSeqMotifOperators(unittest.TestCase):
    def test_seqmotif_operator_to_dict(self):
        seqmotif_operator = seqmotif_operators.SeqMotifOperator(
            pattern_type=seqmotif_operators.PatternType.PROSITE,
            sequence_type=seqmotif_operators.SequenceType.PROTEIN,
            pattern="C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H.")

        self.assertEqual(
            seqmotif_operator._to_dict(), {
                "value": "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H.",
                "pattern_type": "prosite",
                "target": "pdb_protein_sequence"
            })
