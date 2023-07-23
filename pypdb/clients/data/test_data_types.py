"""
Unit tests for the data_types classes.
"""
import unittest
from unittest import mock
import requests
from pypdb.clients.data.graphql.graphql import RSCB_GRAPHQL_URL

from pypdb.clients.data.data_types import DataFetcher, DataType

class TestEntry(unittest.TestCase):
    def test_create(self):
        entry = DataFetcher("4HHB", DataType.ENTRY)

        self.assertTrue(isinstance(entry.properties, dict))
        self.assertTrue(not entry.properties)

        self.assertTrue(isinstance(entry.json_query, dict))
        self.assertTrue(not entry.json_query)

        self.assertTrue(isinstance(entry.response, dict))
        self.assertTrue(not entry.response)

        self.assertEqual(entry.id, ["4HHB"])

    def test_generate_json_query(self):
        entry = DataFetcher("4HHB", DataType.ENTRY)

        property = {"exptl":["method", "details"]}

        entry.add_property(property)

        self.assertIsNotNone(entry.properties)

        entry.generate_json_query()

        self.assertTrue(isinstance(entry.json_query, dict))
        self.assertTrue("query" in entry.json_query)

    def test_fetch_entry(self):
        entry = DataFetcher("4HHB", DataType.ENTRY)
        property = {"exptl":["method", "details"]}

        entry.add_property(property)

        entry.fetch_data()

        self.assertTrue(entry.response)

    def test_return_data_as_pandas_df(self):
        entry = DataFetcher(["4HHB", "12CA", "3PQR"], DataType.ENTRY)
        property = {"exptl":["method", "details"]}

        entry.add_property(property)

        entry.fetch_data()
        df = entry.return_data_as_pandas_df()

        self.assertTrue(df is not None)

    def test_assembly_fetch(self):
        assembly = DataFetcher(["4HHB-1", "12CA-1", "3PQR-1"], DataType.ASSEMBLY)
        property = {"rcsb_assembly_info": ["entry_id", "assembly_id", "polymer_entity_instance_count"]}

        assembly.add_property(property)
        assembly.fetch_data()

        self.assertFalse(not assembly.response)

    def test_polymer_entity_fetch(self):
        fetcher = DataFetcher(["2CPK_1","3WHM_1","2D5Z_1"], DataType.POLYMER_ENTITY)

        property = {"rcsb_id": [], 
                    "rcsb_entity_source_organism": ["ncbi_taxonomy_id", "ncbi_scientific_name"],
                    "rcsb_cluster_membership": ["cluster_id", "identity"]}

        fetcher.add_property(property)
        fetcher.fetch_data()

        self.assertFalse(not fetcher.response)

        df = fetcher.return_data_as_pandas_df()
        self.assertFalse(df is None)

    def test_polymer_instance_fetch(self):
        fetcher = DataFetcher(["4HHB.A", "12CA.A", "3PQR.A"], DataType.POLYMER_ENTITY_INSTANCE)
        property = {"rcsb_id": [],
                    "rcsb_polymer_instance_annotation": ["annotation_id", "name", "type"]}
        fetcher.add_property(property)
        fetcher.fetch_data()

        self.assertFalse(not fetcher.response)

        df = fetcher.return_data_as_pandas_df()
        self.assertFalse(df is None)

    def test_branched_entity_fetch(self):
        fetcher = DataFetcher(["5FMB_2", "6L63_3"], DataType.BRANCHED_ENTITY)
        property = {"pdbx_entity_branch": ["type"],
                    "pdbx_entity_branch_descriptor": ["type", "descriptor"]}

        fetcher.add_property(property)
        fetcher.fetch_data()

        self.assertFalse(not fetcher.response)

        df = fetcher.return_data_as_pandas_df()
        self.assertFalse(df is None)

    def test_chem_comps_fetch(self):
        fetcher = DataFetcher(["NAG","EBW"], DataType.CHEMICAL_COMPONENT)
        property = {"rcsb_id":[], "chem_comp": ["type", "formula_weight","name","formula"],
                    "rcsb_chem_comp_info":["initial_release_date"]}
        fetcher.add_property(property)
        fetcher.fetch_data()
        self.assertFalse(not fetcher.response)

        df = fetcher.return_data_as_pandas_df()
        self.assertFalse(df is None)

if __name__ == '__main__':
    unittest.main()