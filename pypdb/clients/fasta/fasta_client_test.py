"""Tests for RCSB FASTA fetching logic."""
import pytest
import requests
import unittest
from unittest import mock

from pypdb.clients.fasta import fasta_client


class TestFastaLogic(unittest.TestCase):
    @mock.patch.object(requests, "get")
    @mock.patch.object(fasta_client, "_parse_fasta_text_to_list")
    def test_get_fasta_file(self, mock_parse_fasta, mock_get):
        mock_response = mock.Mock()
        mock_response.ok = True
        mock_response.text = "fake_fasta_response"
        mock_get.return_value = mock_response

        fasta_client.get_fasta_from_rcsb_entry("6TML", verbosity=True)
        mock_get.assert_called_once_with(
            "https://www.rcsb.org/fasta/entry/6TML")
        mock_parse_fasta.assert_called_once_with("fake_fasta_response")

    def test_parse_fasta_file(self):

        test_fasta_raw_text = """
>6TML_1|Chains Q7,Q8,Q9,q7,q8,q9|ATPTG11|Toxoplasma gondii (strain ATCC 50853 / GT1) (507601)
MVRNQRYPASPVQEIFLPEPVPFVQFDQTAPSPNSPPAPLPSPSLSQCEEQKDRYR
>6TML_2|Chain i9|ATPTG7|Toxoplasma gondii (strain ATCC 50853 / GT1) (507601)
MPSSSSEDAQGGNRFECVSNSTSPRRKNATKDEAACLQPRRSAVSGPREDVLCIR
>6TML_32|Chains H1,H2,H3,H4|subunit c|Toxoplasma gondii (strain ATCC 50853 / GT1) (507601)
MFFSRLSLSALKAAPAREAL"""

        self.assertEqual(
            fasta_client._parse_fasta_text_to_list(test_fasta_raw_text), [
                fasta_client.FastaSequence(
                    entity_id="6TML_1",
                    chains=["Q7", "Q8", "Q9", "q7", "q8", "q9"],
                    sequence=
                    "MVRNQRYPASPVQEIFLPEPVPFVQFDQTAPSPNSPPAPLPSPSLSQCEEQKDRYR",
                    fasta_header=
                    "6TML_1|Chains Q7,Q8,Q9,q7,q8,q9|ATPTG11|Toxoplasma gondii (strain ATCC 50853 / GT1) (507601)"
                ),
                fasta_client.FastaSequence(
                    entity_id="6TML_2",
                    chains=["i9"],
                    sequence=
                    "MPSSSSEDAQGGNRFECVSNSTSPRRKNATKDEAACLQPRRSAVSGPREDVLCIR",
                    fasta_header=
                    "6TML_2|Chain i9|ATPTG7|Toxoplasma gondii (strain ATCC 50853 / GT1) (507601)"
                ),
                fasta_client.FastaSequence(
                    entity_id="6TML_32",
                    chains=["H1", "H2", "H3", "H4"],
                    sequence="MFFSRLSLSALKAAPAREAL",
                    fasta_header=
                    "6TML_32|Chains H1,H2,H3,H4|subunit c|Toxoplasma gondii (strain ATCC 50853 / GT1) (507601)"
                )
            ])

    def test_get_fasta_from_rcsb_entry_integration(self):
        # This test makes a network call to the RCSB PDB.
        pdb_id = '1EHZ'
        results = fasta_client.get_fasta_from_rcsb_entry(pdb_id)

        self.assertTrue(results) # Check that the list is not empty
        self.assertIsInstance(results, list)

        # Assuming 1EHZ has at least one sequence, typically protein.
        # If 1EHZ changes structure significantly, this test might need updating.
        result = results[0]

        self.assertIsInstance(result, fasta_client.FastaSequence)
        self.assertTrue(result.entity_id.startswith(f'{pdb_id}_'))

        self.assertIsInstance(result.chains, list)
        self.assertTrue(len(result.chains) > 0)
        self.assertIsInstance(result.chains[0], str)

        self.assertIsInstance(result.sequence, str)
        self.assertTrue(len(result.sequence) > 10)

        # Check if all characters in the sequence are valid protein characters
        # (standard amino acids; some PDBs might have non-standard ones, e.g. U for Selenocysteine)
        valid_protein_chars = set('ACDEFGHIKLMNPQRSTVWYU') # Added 'U'
        sequence_chars = set(result.sequence.upper())
        self.assertTrue(sequence_chars.issubset(valid_protein_chars),
                        f"Sequence contains non-standard protein characters: {sequence_chars - valid_protein_chars}")

        self.assertIsInstance(result.fasta_header, str)
        self.assertTrue(result.fasta_header.startswith(f'{pdb_id}_'))
