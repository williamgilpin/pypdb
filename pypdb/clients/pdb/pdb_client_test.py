import gzip
import unittest
from unittest import mock

from pypdb.clients.pdb import pdb_client
from pypdb.util import http_requests


class TestPDBFileDownloading(unittest.TestCase):
    @mock.patch.object(http_requests, "request_limited", autospec=True)
    def test_unsuccessful_test_returns_none(self, mock_http_requests):

        mock_return_value = mock.Mock()
        mock_return_value.ok = False
        mock_http_requests.return_value = mock_return_value

        self.assertIsNone(pdb_client.get_pdb_file("5TML"))
        mock_http_requests.assert_called_once_with(
            "https://files.rcsb.org/download/5TML.pdb")

    @mock.patch.object(http_requests, "request_limited", autospec=True)
    @mock.patch.object(gzip, "decompress")
    def test_compressed_cif_file(self, mock_decompress, mock_http_requests):
        mock_return_value_cif = mock.Mock()
        mock_return_value_cif.ok = True
        mock_return_value_cif.content = "fake_compressed_cif"
        mock_http_requests.return_value = mock_return_value_cif
        mock_decompress.return_value = "fake_decompressed_cif"

        self.assertEqual(
            "fake_decompressed_cif",
            pdb_client.get_pdb_file("1A2B",
                                    pdb_client.PDBFileType.CIF,
                                    compression=True))
        mock_http_requests.assert_called_once_with(
            "https://files.rcsb.org/download/1A2B.cif.gz")
        mock_decompress.assert_called_once_with("fake_compressed_cif")

    @mock.patch.object(http_requests, "request_limited", autospec=True)
    def test_umcompressed_pdb(self, mock_http_requests):
        mock_return_value_pdb = mock.Mock()
        mock_return_value_pdb.text = "fake_uncompressed_pdb"
        mock_return_value_pdb.ok = True
        mock_http_requests.return_value = mock_return_value_pdb

        self.assertEqual("fake_uncompressed_pdb",
                         pdb_client.get_pdb_file("1234"))
        mock_http_requests.assert_called_once_with(
            "https://files.rcsb.org/download/1234.pdb")

    @mock.patch.object(http_requests, "request_limited", autospec=True)
    @mock.patch.object(gzip, "decompress")
    def test_compressed_structfact(self, mock_decompress, mock_http_requests):
        mock_return_value_pdb = mock.Mock()
        mock_return_value_pdb.content = "fake_compressed_structfact"
        mock_return_value_pdb.ok = True
        mock_http_requests.return_value = mock_return_value_pdb
        mock_decompress.return_value = "fake_decompressed_structfact"

        self.assertEqual(
            "fake_decompressed_structfact",
            pdb_client.get_pdb_file("HK97",
                                    pdb_client.PDBFileType.STRUCTFACT,
                                    compression=True))
        mock_http_requests.assert_called_once_with(
            "https://files.rcsb.org/download/HK97-sf.cif.gz")
        mock_decompress.assert_called_once_with("fake_compressed_structfact")

    @mock.patch.object(http_requests, "request_limited", autospec=True)
    def test_uncompressed_xml(self, mock_http_requests):
        mock_return_value_pdb = mock.Mock()
        mock_return_value_pdb.text = "fake_uncompressed_xml"
        mock_return_value_pdb.ok = True
        mock_http_requests.return_value = mock_return_value_pdb

        self.assertEqual(
            "fake_uncompressed_xml",
            pdb_client.get_pdb_file("MI17",
                                    pdb_client.PDBFileType.XML,
                                    compression=False))
        mock_http_requests.assert_called_once_with(
            "https://files.rcsb.org/download/MI17.xml")


if __name__ == '__main__':
    unittest.main()
