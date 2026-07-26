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

    def test_chemical_formula_operator_to_dict(self):
        # Thiamine
        formula_operator = chemical_operators.ChemicalFormulaOperator(
            formula="C12H17N4OS")

        self.assertEqual(formula_operator._to_dict(), {
            "value": "C12H17N4OS",
            "type": "formula",
            "match_subset": False
        })

    def test_chemical_formula_operator_with_match_subset(self):
        formula_operator = chemical_operators.ChemicalFormulaOperator(
            formula="C12H17N4OS", match_subset=True)

        self.assertEqual(formula_operator._to_dict(), {
            "value": "C12H17N4OS",
            "type": "formula",
            "match_subset": True
        })

    def test_chemical_similarity_operator_with_smiles(self):
        similarity_operator = chemical_operators.ChemicalSimilarityOperator(
            descriptor="CC(=O)NC1C(O)OC(CO)C(O)C1O")

        self.assertEqual(
            similarity_operator._to_dict(), {
                "value": "CC(=O)NC1C(O)OC(CO)C(O)C1O",
                "type": "descriptor",
                "descriptor_type": "SMILES",
                "match_type": "fingerprint-similarity"
            })

    def test_chemical_similarity_operator_with_inchi(self):
        similarity_operator = chemical_operators.ChemicalSimilarityOperator(
            descriptor="InChI=1S/C10H26N4/c11-5-3-9-13-7-1-2-8-14-10-4-6-12")

        self.assertEqual(similarity_operator.descriptor_type, "InChI")
