"""Tests for RCSB ligand coordinate download logic."""

import unittest
from unittest import mock
import warnings

from pypdb.clients.ligand import ligand_client
from pypdb.util import http_requests


class TestIdealLigandDownloading(unittest.TestCase):

    @mock.patch.object(http_requests, "request_limited")
    def test_ideal_sdf_request(self, mock_request):
        mock_response = mock.MagicMock()
        mock_response.ok = True
        mock_response.text = "ATP\n  CCTOOLS\n"
        mock_request.return_value = mock_response

        result = ligand_client.get_ideal_ligand_file("atp", verbosity=False)

        mock_request.assert_called_once_with(
            "https://files.rcsb.org/ligands/download/ATP_ideal.sdf")
        self.assertEqual(result, "ATP\n  CCTOOLS\n")

    def test_ideal_mol2_is_unsupported(self):
        with self.assertWarns(UserWarning):
            result = ligand_client.get_ideal_ligand_file(
                "ATP",
                filetype=ligand_client.LigandFileType.MOL2,
                verbosity=False)

        self.assertIsNone(result)

    @mock.patch.object(http_requests, "request_limited")
    def test_failed_request_returns_none(self, mock_request):
        mock_request.return_value = None

        with self.assertWarns(UserWarning):
            result = ligand_client.get_ideal_ligand_file("ATP",
                                                         verbosity=False)

        self.assertIsNone(result)


class TestLigandInstanceDownloading(unittest.TestCase):

    @mock.patch.object(http_requests, "request_limited")
    def test_instance_sdf_request(self, mock_request):
        mock_response = mock.MagicMock()
        mock_response.ok = True
        mock_response.text = "HEM\n  ModelServer\n"
        mock_request.return_value = mock_response

        result = ligand_client.get_ligand_instance_file("4HHB",
                                                        "A",
                                                        142,
                                                        verbosity=False)

        mock_request.assert_called_once_with(
            "https://models.rcsb.org/v1/4hhb/ligand"
            "?auth_asym_id=A&auth_seq_id=142&encoding=sdf")
        self.assertEqual(result, "HEM\n  ModelServer\n")

    @mock.patch.object(http_requests, "request_limited")
    def test_instance_mol2_request(self, mock_request):
        mock_response = mock.MagicMock()
        mock_response.ok = True
        mock_response.text = "# Name: HEM\n"
        mock_request.return_value = mock_response

        ligand_client.get_ligand_instance_file(
            "4HHB",
            "A",
            142,
            filetype=ligand_client.LigandFileType.MOL2,
            verbosity=False)

        self.assertIn("encoding=mol2", mock_request.call_args[0][0])

    @mock.patch.object(http_requests, "request_limited")
    def test_empty_result_returns_none(self, mock_request):
        # ModelServer returns HTTP 200 with no atoms when nothing matches
        mock_response = mock.MagicMock()
        mock_response.ok = True
        mock_response.text = (
            "> <model_server_stats.element_count>\n0\n\n$$$$\n")
        mock_request.return_value = mock_response

        with self.assertWarns(UserWarning):
            result = ligand_client.get_ligand_instance_file("4HHB",
                                                            "Z",
                                                            9999,
                                                            verbosity=False)

        self.assertIsNone(result)


class TestLigandDownloadIntegration(unittest.TestCase):
    """Tests that hit the live RCSB service."""

    def test_fetches_ideal_ccd_coordinates(self):
        result = ligand_client.get_ideal_ligand_file("ATP", verbosity=False)

        self.assertIsNotNone(result)
        self.assertEqual(result.splitlines()[0].strip(), "ATP")
        # The counts line of an SDF records the number of atoms and bonds
        self.assertIn("47 49", result.splitlines()[3])

    def test_fetches_instance_coordinates(self):
        result = ligand_client.get_ligand_instance_file("4HHB",
                                                        "A",
                                                        142,
                                                        verbosity=False)

        self.assertIsNotNone(result)
        self.assertEqual(result.splitlines()[0].strip(), "HEM")
        # Heme has 43 atoms and 50 bonds
        self.assertIn("43 50", result.splitlines()[3])

    def test_fetches_instance_coordinates_as_mol2(self):
        result = ligand_client.get_ligand_instance_file(
            "4HHB",
            "A",
            142,
            filetype=ligand_client.LigandFileType.MOL2,
            verbosity=False)

        self.assertIsNotNone(result)
        self.assertIn("HEM", result.splitlines()[0])


if __name__ == '__main__':
    unittest.main()
