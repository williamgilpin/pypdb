"""Tests for RCSB SeqMotif Search Service Operators."""

import unittest

from pypdb.clients.search.operators import chemical_operators


class TestChemicalOperators(unittest.TestCase):
    def test_chemical_operator_to_dict(self):
        # InChI
        inchi_operator = chemical_operators.ChemicalOperator(
            # Panadol
            descriptor=
            "InChI=1S/C8H9NO2/c1-6(10)9-7-2-4-8(11)5-3-7/h2-5,11H,1H3,(H,9,10)",
            matching_criterion=chemical_operators.DescriptorMatchingCriterion.
            GRAPH_RELAXED_STEREO)
        self.assertEqual(inchi_operator.descriptor_type, "InChI")
        self.assertEqual(
            inchi_operator._to_dict(), {
                "value":
                "InChI=1S/C8H9NO2/c1-6(10)9-7-2-4-8(11)5-3-7/h2-5,11H,1H3,(H,9,10)",
                "type": "descriptor",
                "descriptor_type": "InChI",
                "match_type": "graph-relaxed-stereo"
            })

        # SMILES
        smiles_operator = chemical_operators.ChemicalOperator(
            descriptor=
            "CC(C)C[C@H](NC(=O)OCC1CCC(F)(F)CC1)C(=O)N[C@@H](C[C@@H]2CCNC2=O)[C@@H](O)[S](O)(=O)=O"
        )
        self.assertEqual(
            smiles_operator._to_dict(), {
                "value":
                "CC(C)C[C@H](NC(=O)OCC1CCC(F)(F)CC1)C(=O)N[C@@H](C[C@@H]2CCNC2=O)[C@@H](O)[S](O)(=O)=O",
                "type": "descriptor",
                "descriptor_type": "SMILES",
                "match_type": "graph-strict"
            })
