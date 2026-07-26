"""Tests for RCSB structure motif search operators."""

import unittest

from pypdb.clients.search.operators import strucmotif_operators


class TestStrucMotifOperator(unittest.TestCase):

    def _residues(self, count=3):
        return [
            strucmotif_operators.StructureMotifResidue(label_asym_id="A",
                                                       label_seq_id=160 + i)
            for i in range(count)
        ]

    def test_minimal_motif(self):
        operator = strucmotif_operators.StrucMotifOperator(
            pdb_entry_id="2MNR", residues=self._residues())

        self.assertEqual(
            operator._to_dict(), {
                "value": {
                    "entry_id":
                    "2MNR",
                    "residue_ids": [{
                        "label_asym_id": "A",
                        "label_seq_id": 160
                    }, {
                        "label_asym_id": "A",
                        "label_seq_id": 161
                    }, {
                        "label_asym_id": "A",
                        "label_seq_id": 162
                    }]
                }
            })

    def test_motif_with_rmsd_and_pairing_scheme(self):
        operator = strucmotif_operators.StrucMotifOperator(
            pdb_entry_id="2MNR",
            residues=self._residues(),
            rmsd_cutoff=2,
            atom_pairing_scheme=strucmotif_operators.AtomPairingScheme.BACKBONE)

        result = operator._to_dict()

        self.assertEqual(result["rmsd_cutoff"], 2)
        self.assertEqual(result["atom_pairing_scheme"], "BACKBONE")

    def test_motif_with_exchanges(self):
        residues = [
            strucmotif_operators.StructureMotifResidue(
                label_asym_id="A", label_seq_id=162, exchanges=["LYS", "HIS"]),
            strucmotif_operators.StructureMotifResidue(label_asym_id="A",
                                                       label_seq_id=193),
        ]
        operator = strucmotif_operators.StrucMotifOperator(
            pdb_entry_id="2MNR", residues=residues)

        result = operator._to_dict()

        self.assertEqual(result["exchanges"], [{
            "residue_id": {
                "label_asym_id": "A",
                "label_seq_id": 162
            },
            "allowed": ["LYS", "HIS"]
        }])

    def test_too_few_residues_is_rejected(self):
        with self.assertRaises(
                strucmotif_operators.InvalidMotifResidueCountError):
            strucmotif_operators.StrucMotifOperator(pdb_entry_id="2MNR",
                                                    residues=self._residues(1))

    def test_too_many_residues_is_rejected(self):
        with self.assertRaises(
                strucmotif_operators.InvalidMotifResidueCountError):
            strucmotif_operators.StrucMotifOperator(
                pdb_entry_id="2MNR", residues=self._residues(11))


if __name__ == '__main__':
    unittest.main()
