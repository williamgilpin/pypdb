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

    def test_get_pdb_file_pdb_format_integration(self):
        # Test PDB format, uncompressed
        pdb_content_uncompressed = pdb_client.get_pdb_file('1EHZ', filetype=pdb_client.PDBFileType.PDB, compression=False)
        self.assertTrue(pdb_content_uncompressed)
        self.assertIsInstance(pdb_content_uncompressed, str)
        self.assertTrue(pdb_content_uncompressed.startswith("HEADER"))

        # Test PDB format, compressed (library handles decompression)
        pdb_content_compressed = pdb_client.get_pdb_file('1EHZ', filetype=pdb_client.PDBFileType.PDB, compression=True)
        self.assertTrue(pdb_content_compressed)
        self.assertIsInstance(pdb_content_compressed, str)
        self.assertTrue(pdb_content_compressed.startswith("HEADER"))

    def test_get_pdb_file_cif_format_integration(self):
        # Test CIF format, uncompressed
        cif_content_uncompressed = pdb_client.get_pdb_file('1EHZ', filetype=pdb_client.PDBFileType.CIF, compression=False)
        self.assertTrue(cif_content_uncompressed)
        self.assertIsInstance(cif_content_uncompressed, str)
        self.assertIn("_entry.id   1EHZ", cif_content_uncompressed)

        # Test CIF format, compressed (library handles decompression)
        cif_content_compressed = pdb_client.get_pdb_file('1EHZ', filetype=pdb_client.PDBFileType.CIF, compression=True)
        self.assertTrue(cif_content_compressed)
        self.assertIsInstance(cif_content_compressed, str)
        self.assertIn("_entry.id   1EHZ", cif_content_compressed)

    def test_get_pdb_file_xml_format_integration(self):
        # Test XML format, uncompressed
        xml_content_uncompressed = pdb_client.get_pdb_file('1EHZ', filetype=pdb_client.PDBFileType.XML, compression=False)
        self.assertTrue(xml_content_uncompressed)
        self.assertIsInstance(xml_content_uncompressed, str)
        self.assertIn("<?xml version=", xml_content_uncompressed)
        self.assertIn("<pdbx:datablock ", xml_content_uncompressed) # Check for common PDBx/XML namespace

        # Test XML format, compressed (library handles decompression)
        xml_content_compressed = pdb_client.get_pdb_file('1EHZ', filetype=pdb_client.PDBFileType.XML, compression=True)
        self.assertTrue(xml_content_compressed)
        self.assertIsInstance(xml_content_compressed, str)
        self.assertIn("<?xml version=", xml_content_compressed)
        self.assertIn("<pdbx:datablock ", xml_content_compressed)

    def test_get_pdb_file_structfact_integration(self):
        pdb_id = '1AKI' # This PDB ID is known to have structure factors
        # Test STRUCTFACT format, uncompressed
        sf_content_uncompressed = pdb_client.get_pdb_file(pdb_id, filetype=pdb_client.PDBFileType.STRUCTFACT, compression=False)
        self.assertTrue(sf_content_uncompressed)
        self.assertIsInstance(sf_content_uncompressed, str)
        self.assertIn(f"_entry.id   {pdb_id}", sf_content_uncompressed)
        self.assertIn("loop_", sf_content_uncompressed) # Structure factor files are in mmCIF format
        self.assertIn("_refln.", sf_content_uncompressed) # Common reflection data items

        # Test STRUCTFACT format, compressed (library handles decompression)
        sf_content_compressed = pdb_client.get_pdb_file(pdb_id, filetype=pdb_client.PDBFileType.STRUCTFACT, compression=True)
        self.assertTrue(sf_content_compressed)
        self.assertIsInstance(sf_content_compressed, str)
        self.assertIn(f"_entry.id   {pdb_id}", sf_content_compressed)
        self.assertIn("loop_", sf_content_compressed)
        self.assertIn("_refln.", sf_content_compressed)


if __name__ == '__main__':
    unittest.main()
