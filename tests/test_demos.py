"""Tests covering every case demonstrated in `demos/demos.ipynb`.

Each test here corresponds to a cell of the demos notebook, so that the
documented usage of the package stays verified. These tests hit the live RCSB
API; result counts change as the PDB grows, so assertions check for
well-known entries and result shapes rather than exact result sets.
"""

import unittest

## Import from local directory
import sys
sys.path.insert(0, '../pypdb')
from pypdb import *

from pypdb.clients.search.search_client import LogicalOperator
from pypdb.clients.search.search_client import QueryGroup
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import perform_search_with_graph
from pypdb.clients.search.operators import text_operators


class TestDemoSearchQueries(unittest.TestCase):
    """Covers the "Search functions that return lists of PDB IDs" section."""

    def test_search_by_keyword(self):
        found_pdbs = Query("ribosome").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("4V5A", found_pdbs)

    def test_search_by_pubmed_id(self):
        found_pdbs = Query(27499440, "PubmedIdQuery").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("5IMT", found_pdbs)

    def test_search_by_taxonomy_id(self):
        found_pdbs = Query('6239', 'TreeEntityQuery').search()  # C. elegans
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("1D4X", found_pdbs)

    def test_search_by_experimental_method(self):
        found_pdbs = Query('SOLID-STATE NMR',
                           query_type='ExpTypeQuery').search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("1CEK", found_pdbs)

    def test_search_by_structure_similarity(self):
        found_pdbs = Query('2E8D', query_type="structure").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("2E8D", found_pdbs)

    def test_search_by_author(self):
        found_pdbs = Query('Perutz, M.F.',
                           query_type='AdvancedAuthorQuery').search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("2HHB", found_pdbs)

    def test_search_by_organism(self):
        found_pdbs = Query("Dictyostelium", query_type="OrganismQuery").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("2H84", found_pdbs)

    def test_search_by_uniprot_id(self):
        found_pdbs = Query("A0A023GPI8", query_type="uniprot").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("4K1Y", found_pdbs)

    def test_search_by_pfam_number(self):
        found_pdbs = Query("PF00008", query_type="pfam").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("1A3P", found_pdbs)

    def test_search_by_ligand_smiles(self):
        smiles = "Clc1nc(Br)nc2nc[nH]c12"
        found_pdbs = Query(smiles, query_type="chemical").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertTrue(all(isinstance(pdb_id, str) for pdb_id in found_pdbs))

    def test_search_by_ligand_component_id(self):
        # Short chemical terms are chemical component (ligand) IDs
        found_pdbs = Query("NAG", query_type="chemical").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("5FYJ", found_pdbs)

    def test_negated_text_operator(self):
        negated_title_search = text_operators.ContainsWordsOperator(
            attribute="struct.title",
            value="kinase inhibitor",
            negation=True,
        )

        found_pdbs = perform_search(negated_title_search,
                                    ReturnType.ENTRY,
                                    verbosity=False)

        self.assertTrue(len(found_pdbs) > 0)
        # A structure whose title is about kinase inhibition should be excluded
        kinase_inhibitor_pdbs = perform_search(
            text_operators.ContainsPhraseOperator(
                attribute="struct.title", value="kinase inhibitor"),
            ReturnType.ENTRY,
            verbosity=False)
        self.assertTrue(len(kinase_inhibitor_pdbs) > 0)
        self.assertTrue(
            set(found_pdbs).isdisjoint(set(kinase_inhibitor_pdbs)),
            "Negated search should exclude all non-negated search hits")


class TestDemoInformationSearch(unittest.TestCase):
    """Covers the "Information Search functions" section."""

    def test_find_papers(self):
        matching_papers = find_papers('crispr', max_results=10)
        self.assertTrue(len(matching_papers) > 0)
        self.assertTrue(
            all(isinstance(paper, str) for paper in matching_papers))


class TestDemoSinglePDBFunctions(unittest.TestCase):
    """Covers "Functions that return information about single PDB IDs"."""

    def test_get_pdb_file_as_cif(self):
        pdb_file = get_pdb_file('4lza', filetype='cif', compression=False)
        self.assertIsNotNone(pdb_file)
        self.assertTrue(pdb_file.startswith("data_4LZA"))

    def test_get_info(self):
        all_info = get_info('4LZA')
        self.assertIsNotNone(all_info)
        self.assertEqual(all_info['rcsb_id'], '4LZA')
        self.assertIn('struct', all_info)

    def test_sequence_search(self):
        q = Query(
            "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK"
            "VADALTAVAHVDDMPNAL",
            query_type="sequence",
            return_type="polymer_entity")

        results = q.search()

        self.assertIsNotNone(results)
        # `polymer_entity` searches return the raw JSON response
        self.assertIn("result_set", results)
        self.assertTrue(len(results["result_set"]) > 0)


class TestDemoAdvancedSearchAPI(unittest.TestCase):
    """Covers the "New API for advanced search" section."""

    def test_default_operator_entry_search(self):
        results = perform_search(text_operators.DefaultOperator(
            value="ribosome"),
                                 ReturnType.ENTRY,
                                 verbosity=False)

        self.assertTrue(len(results) > 0)
        self.assertIn("4V5A", results)

    def test_exact_match_polymer_entity_search(self):
        results = perform_search(text_operators.ExactMatchOperator(
            value="Mus musculus",
            attribute="rcsb_entity_source_organism.taxonomy_lineage.name"),
                                 ReturnType.POLYMER_ENTITY,
                                 verbosity=False)

        self.assertTrue(len(results) > 0)
        # Polymer entity IDs look like `10EN_4`
        self.assertTrue(all("_" in result for result in results[:10]))

    def test_in_operator_non_polymer_entity_search(self):
        results = perform_search(text_operators.InOperator(
            values=["Mus musculus", "Homo sapiens"],
            attribute="rcsb_entity_source_organism.taxonomy_lineage.name"),
                                 ReturnType.NON_POLYMER_ENTITY,
                                 verbosity=False)

        self.assertTrue(len(results) > 0)
        self.assertTrue(all("_" in result for result in results[:10]))

    def test_contains_words_polymer_instance_search(self):
        results = perform_search(text_operators.ContainsWordsOperator(
            value="actin-binding protein", attribute="struct.title"),
                                 ReturnType.POLYMER_INSTANCE,
                                 verbosity=False)

        self.assertTrue(len(results) > 0)
        # Polymer instance IDs look like `1HQZ.B`
        self.assertTrue(all("." in result for result in results[:10]))

    def test_contains_phrase_assembly_search(self):
        results = perform_search(text_operators.ContainsPhraseOperator(
            value="actin-binding protein", attribute="struct.title"),
                                 ReturnType.ASSEMBLY,
                                 verbosity=False)

        self.assertTrue(len(results) > 0)
        # Assembly IDs look like `1HQZ-1`
        self.assertTrue(all("-" in result for result in results[:10]))

    def test_comparison_operator_release_date(self):
        results = perform_search(text_operators.ComparisonOperator(
            value="2019-01-01T00:00:00Z",
            attribute="rcsb_accession_info.initial_release_date",
            comparison_type=text_operators.ComparisonType.GREATER),
                                 ReturnType.ENTRY,
                                 verbosity=False)

        self.assertTrue(len(results) > 0)

    def test_range_operator_release_date(self):
        results = perform_search(text_operators.RangeOperator(
            from_value="2019-01-01T00:00:00Z",
            to_value="2020-01-01T00:00:00Z",
            include_lower=True,
            include_upper=False,
            attribute="rcsb_accession_info.initial_release_date"),
                                 ReturnType.ENTRY,
                                 verbosity=False)

        self.assertTrue(len(results) > 0)

    def test_range_operator_cell_length_suppressed_output(self):
        cell_a_operator = text_operators.RangeOperator(
            attribute='cell.length_a',
            from_value=80,
            to_value=84,
            include_upper=True)

        results = perform_search_with_graph(
            query_object=cell_a_operator,
            return_type=ReturnType.ENTRY,
            verbosity=False,
        )

        self.assertTrue(len(results) > 0)

    def test_comparison_operator_resolution(self):
        results = perform_search(text_operators.ComparisonOperator(
            value=4,
            attribute="rcsb_entry_info.resolution_combined",
            comparison_type=text_operators.ComparisonType.LESS),
                                 ReturnType.ENTRY,
                                 verbosity=False)

        self.assertTrue(len(results) > 0)

    def test_exists_operator(self):
        results = perform_search(text_operators.ExistsOperator(
            attribute="rcsb_accession_info.initial_release_date"),
                                 ReturnType.ENTRY,
                                 verbosity=False)

        self.assertTrue(len(results) > 0)

    def test_graph_search_with_nested_query_groups(self):
        under_4A_resolution_operator = text_operators.ComparisonOperator(
            value=4,
            attribute="rcsb_entry_info.resolution_combined",
            comparison_type=text_operators.ComparisonType.GREATER)

        is_mus_operator = text_operators.ExactMatchOperator(
            value="Mus musculus",
            attribute="rcsb_entity_source_organism.taxonomy_lineage.name")

        is_human_operator = text_operators.ExactMatchOperator(
            value="Homo sapiens",
            attribute="rcsb_entity_source_organism.taxonomy_lineage.name")

        is_human_or_mus_group = QueryGroup(
            queries=[is_mus_operator, is_human_operator],
            logical_operator=LogicalOperator.OR)

        is_under_4A_and_human_or_mus_group = QueryGroup(
            queries=[is_human_or_mus_group, under_4A_resolution_operator],
            logical_operator=LogicalOperator.AND)

        results = perform_search_with_graph(
            query_object=is_under_4A_and_human_or_mus_group,
            return_type=ReturnType.ENTRY,
            verbosity=False)

        self.assertTrue(len(results) > 0)


if __name__ == '__main__':
    unittest.main()
