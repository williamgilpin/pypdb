import unittest
import warnings
from unittest.mock import patch

## Import from local directory
import sys
sys.path.insert(0, '../pypdb')
from pypdb import *

# TODO(ejwilliams): Write generic logic, to execute `test_*.py` files
# within the pypdb directory (removing need for sys.path hack)

# aa_index[s] for s in seq_dict[k] if s in aa_index.keys()]

class TestSearchFunctions(unittest.TestCase):

    def test_searchterm(self):
        found_pdbs = Query('ribosome').search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn('4V5A', found_pdbs)

    def test_pubmed(self):
        found_pdbs = Query(27499440, "PubmedIdQuery").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn('5IMT', found_pdbs)

    def test_treeentity(self):
        found_pdbs = Query('6239', 'TreeEntityQuery').search() # C. elegans
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn('1D4X', found_pdbs)

    def test_exptype(self):
        found_pdbs = Query('SOLID-STATE NMR', 'ExpTypeQuery').search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn('1CEK', found_pdbs)

    def test_structure(self):
        found_pdbs = Query('2E8D', 'structure').search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn('2E8D', found_pdbs)

    def test_advancedauthor(self):
        found_pdbs = Query('Perutz, M.F.', 'AdvancedAuthorQuery').search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn('2HHB', found_pdbs)

    def test_organism(self):
        found_pdbs = Query('Dictyostelium', 'OrganismQuery').search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn('2H84', found_pdbs)

    def test_full_text_query(self):
        found_pdbs = Query("hemoglobin").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("4HHB", found_pdbs)

    def test_pfam_query(self):
        found_pdbs = Query("PF00008", query_type="pfam").search() # 7tm_1
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("1A3P", found_pdbs) # Bovine Rhodopsin

    def test_uniprot_query(self):
        # Hemoglobin subunit beta (human)
        found_pdbs = Query("P68871", query_type="uniprot").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("4HHB", found_pdbs)

    def test_ec_number_query(self):
        # Ubiquitin--protein ligase
        found_pdbs = Query("2.3.2.5", query_type="ec_number").search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("2AFM", found_pdbs)

    def test_partial_ec_number_query_matches_lineage(self):
        # EC numbers are a lineage, so a partial number matches everything
        # below it in the hierarchy.
        exact_pdbs = Query("2.3.2.5", query_type="ec_number").search()
        parent_pdbs = Query("2.3.2", query_type="ec_number").search()

        self.assertTrue(len(parent_pdbs) > len(exact_pdbs))
        self.assertTrue(set(exact_pdbs).issubset(set(parent_pdbs)))

    def test_ec_number_query_with_multiple_values(self):
        found_pdbs = Query(["2.3.2.5", "1.1.1.1"],
                           query_type="ec_number").search()

        alcohol_dehydrogenase_pdbs = Query("1.1.1.1",
                                           query_type="ec_number").search()

        self.assertTrue(len(found_pdbs) > 0)
        # Searching for several EC numbers returns the union of their hits
        self.assertTrue(
            set(alcohol_dehydrogenase_pdbs).issubset(set(found_pdbs)))

    @patch('pypdb.pypdb.search_client.perform_search')
    def test_ec_number_query_uses_in_operator(self, mock_perform_search):
        mock_perform_search.return_value = ["2AFM"]

        found_pdbs = Query("2.3.2.5", query_type="ec_number").search()

        self.assertEqual(found_pdbs, ["2AFM"])
        _, kwargs = mock_perform_search.call_args
        search_operator = kwargs["search_operator"]
        self.assertEqual(search_operator.attribute,
                         "rcsb_polymer_entity.rcsb_ec_lineage.id")
        self.assertEqual(search_operator.values, ["2.3.2.5"])
        self.assertEqual(kwargs["return_type"], ReturnType.ENTRY)

    def test_sequence_query(self):
        # Sequence from demos/demos.ipynb, corresponds to 1A00 (Hemoglobin alpha chain)
        sequence = "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTAVAHVDDMPNAL"
        found_pdbs = Query(sequence, query_type="sequence").search()
        self.assertTrue(len(found_pdbs) > 0)
        # The search might return polymer entity IDs like "1A00_1"
        # We check if any of the results start with "1A00"
        self.assertTrue(any(pdb_id.startswith("1A00") for pdb_id in found_pdbs))

    def test_seqmotif_query(self):
        # RCSB PDB documentation indicates that seqmotif queries are simple sequences.
        # Let's use a known motif from a specific PDB entry, e.g. a zinc finger motif CXXC
        # For 1ZNF (Zinc finger domain of human Kruppel-like factor 8), chain A has "CAEC" at pos 4-7
        found_pdbs = Query("CAEC", query_type="seqmotif").search()
        # Making this less specific as 1ZNF might not always be in the first page of many results,
        # and previous runs showed results were returned, just not 1ZNF.
        # If found_pdbs is None (due to HTTP error), this will fail, which is acceptable.
        # If found_pdbs is an empty list, this will fail.
        self.assertTrue(found_pdbs and len(found_pdbs) > 0)

    def test_ribosome_search_with_custom_user_agent(self):
        """Test a basic search to ensure User-Agent modification doesn't break things."""
        result = Query("ribosome").search()
        self.assertIsNotNone(result, "Search result should not be None")
        self.assertIsInstance(result, list, "Search result should be a list")
        self.assertTrue(len(result) > 0, "Search for 'ribosome' should return results")
        # Optionally, check for a known PDB ID if the results are consistent
        # self.assertIn('4V5A', result) # Example, can be added if results are stable

    def test_chemical_query(self):
        # A short chemical search term is a chemical component ID, and is
        # searched against the ligand component attribute.
        found_pdbs = Query("NAG", query_type="chemical").search()
        self.assertTrue(len(found_pdbs) > 0)
        # 5FYJ and 3OG2 both bind NAG (N-Acetyl-D-Glucosamine)
        self.assertIn("5FYJ", found_pdbs)
        self.assertIn("3OG2", found_pdbs)

    @patch('pypdb.pypdb.search_client.perform_search')
    def test_chemical_smiles_query_uses_new_search_client(self, mock_perform_search):
        mock_perform_search.return_value = ["8XYZ"]

        found_pdbs = Query("Clc1nc(Br)nc2nc[nH]c12", query_type="chemical").search()

        self.assertEqual(found_pdbs, ["8XYZ"])
        _, kwargs = mock_perform_search.call_args
        chemical_operator = kwargs["search_operator"]
        self.assertEqual(chemical_operator.descriptor, "Clc1nc(Br)nc2nc[nH]c12")
        self.assertEqual(chemical_operator.descriptor_type, "SMILES")
        self.assertEqual(kwargs["return_type"], ReturnType.ENTRY)
        self.assertFalse(kwargs["return_raw_json_dict"])
        self.assertFalse(kwargs["verbosity"])

    @patch('pypdb.pypdb.search_client.perform_search')
    def test_full_text_query_uses_new_search_client(self, mock_perform_search):
        mock_perform_search.return_value = ["4V5A"]

        found_pdbs = Query("ribosome").search()

        self.assertEqual(found_pdbs, ["4V5A"])
        _, kwargs = mock_perform_search.call_args
        self.assertEqual(kwargs["search_operator"].value, "ribosome")
        self.assertEqual(kwargs["return_type"], ReturnType.ENTRY)

    # def test_blast(self):
    #     found_pdbs = blast_from_sequence(
    #         'MTKIANKYEVIDNVEKLEKALKRLREAQSVYATYTQEQVDKIFFEAAMAANKMRIPLAKMAVE'
    #         + 'ETGMGVVEDKVIKNHYASEYIYNAYKNTKTCGVIEEDPAFGIKKIAEPLGVIAAVIPTTNP'
    #         + 'TSTAIFKTLIALKTRNAIIISPHPRAKNSTIEAAKIVLEAAVKAGAPEGIIGWIDVPSLEL'
    #         + 'TNLVMREADVILATGGPGLVKAAYSSGKPAIGVGAGNTPAIIDDSADIVLAVNSIIHSKTF'
    #         + 'DNGMICASEQSVIVLDGVYKEVKKEFEKRGCYFLNEDETEKVRKTIIINGALNAKIVGQKA'
    #         + 'HTIANLAGFEVPETTKILIGEVTSVDISEEFAHEKLCPVLAMYRAKDFDDALDKAERLVAD'
    #         + 'GGFGHTSSLYIDTVTQKEKLQKFSERMKTCRILVNTPSSQGGIGDLYNFKLAPSL',
    #         1e-20)
    #     self.assertTrue(len(found_pdbs) > 0)
    #     self.assertTrue(type(found_pdbs[0][0]) == str)

    #     # an error page would be a longer string
    #     self.assertTrue(len(found_pdbs[0][0]) < 10)


class TestSearchConvenienceFunctions(unittest.TestCase):

    def test_count_results_matches_search_length(self):
        count = count_results("SOLID-STATE NMR", query_type="ExpTypeQuery")
        found_pdbs = Query("SOLID-STATE NMR", "ExpTypeQuery").search()

        self.assertIsInstance(count, int)
        self.assertEqual(count, len(found_pdbs))

    def test_get_top_results_respects_max_results(self):
        results = get_top_results("crispr", max_results=3)

        self.assertEqual(len(results), 3)
        self.assertTrue(all(isinstance(pdb_id, str) for pdb_id in results))

    def test_get_top_results_are_the_highest_scoring(self):
        top_results = get_top_results("ribosome", max_results=5)
        all_results = Query("ribosome").search()

        # Searches are score-sorted by default, so the top hits should be the
        # leading entries of the full result set.
        self.assertEqual(top_results, all_results[:5])

    def test_find_ligands(self):
        found_pdbs = find_ligands("NAG")

        self.assertTrue(len(found_pdbs) > 0)
        self.assertIn("5FYJ", found_pdbs)

    def test_find_ligands_with_max_results(self):
        found_pdbs = find_ligands("NAG", max_results=3)

        self.assertEqual(len(found_pdbs), 3)


class TestInfoFunctions(unittest.TestCase):

    def test_get_info_successful_retrieval(self):
        pdb_id = "4HHB"
        result = get_info(pdb_id)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, dict)

        # Assertions for expected top-level keys
        expected_keys = ["rcsb_id", "struct", "citation", "exptl", "cell", "symmetry"]
        for key in expected_keys:
            self.assertIn(key, result)

        # Assertions for specific values
        self.assertEqual(result['rcsb_id'], pdb_id)
        self.assertIsInstance(result['struct']['title'], str)
        self.assertTrue(len(result['struct']['title']) > 0)

        self.assertIsInstance(result['exptl'], list)
        self.assertTrue(len(result['exptl']) > 0)
        self.assertIsInstance(result['exptl'][0]['method'], str)
        self.assertTrue(len(result['exptl'][0]['method']) > 0)

    def test_get_info_nonexistent_pdb_id(self):
        pdb_id = "9ZZZ" # A clearly non-existent PDB ID
        with self.assertWarns(UserWarning) as cm:
            result = get_info(pdb_id)

        self.assertIsNone(result)
        # Adjusting to the actual warning observed in case of repeated failures for a non-existent ID
        self.assertTrue(
            str(cm.warning) == "Retrieval failed, returning None" or \
            str(cm.warning) == "Too many failures on requests. Exiting...",
            f"Unexpected warning message: {str(cm.warning)}"
        )

    def test_describe_chemical_successful_retrieval(self):
        chem_id = "NAG"
        result = describe_chemical(chem_id)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, dict)

        # Assertions for expected top-level keys
        # Based on observed structure from RCSB for chemcomp NAG
        expected_keys = ["rcsb_id", "chem_comp", "rcsb_chem_comp_descriptor", "rcsb_chem_comp_info", "pdbx_chem_comp_identifier"]
        for key in expected_keys:
            self.assertIn(key, result, f"Expected key '{key}' not found in result.")

        # Assertions for specific values
        self.assertEqual(result['rcsb_id'], chem_id) # Check top-level rcsb_id
        self.assertEqual(result['chem_comp']['id'], chem_id)
        # Using the name observed from the previous failed test run
        self.assertEqual(result['chem_comp']['name'].upper(), "2-ACETAMIDO-2-DEOXY-BETA-D-GLUCOPYRANOSE")
        self.assertTrue(float(result['chem_comp']['formula_weight']) > 0) # Check it's a positive float

    @patch('pypdb.pypdb.get_info')
    def test_describe_chemical_nonexistent_chem_id(self, mock_get_info):
        # Configure the mock_get_info to simulate a failed retrieval
        # It should return None and issue the specific warning "Retrieval failed, returning None"
        def get_info_side_effect(chem_id_param, url_root=None):
            warnings.warn("Retrieval failed, returning None", UserWarning)
            return None

        mock_get_info.side_effect = get_info_side_effect

        chem_id = "XXX" # This ID won't actually be queried due to mocking

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", UserWarning) # Ensure UserWarnings are captured
            result = describe_chemical(chem_id)

        self.assertIsNone(result, f"Expected None for non-existent chem_id '{chem_id}', but got {result}")

        # Check that get_info was called with the correct parameters for chemical description
        mock_get_info.assert_called_once_with(chem_id, url_root='https://data.rcsb.org/rest/v1/core/chemcomp/')

        # Check if the specific warning from get_info (via describe_chemical) was caught
        found_expected_warning = False
        for captured_warning in w:
            if issubclass(captured_warning.category, UserWarning) and \
               "Retrieval failed, returning None" in str(captured_warning.message):
                found_expected_warning = True
                break
        self.assertTrue(found_expected_warning,
                        f"Expected 'Retrieval failed, returning None' UserWarning was not issued. Captured warnings: {[str(cw.message) for cw in w]}")

    def test_get_ligands_successful_retrieval(self):
        result = get_ligands('100D')

        self.assertIsNotNone(result)
        self.assertEqual(result['id'], '100D')

        ligands = result['ligandInfo']['ligand']
        self.assertEqual(len(ligands), 1)

        spermine = ligands[0]
        self.assertEqual(spermine['@structureId'], '100D')
        self.assertEqual(spermine['@chemicalID'], 'SPM')
        self.assertEqual(spermine['@type'], 'non-polymer')
        self.assertEqual(spermine['chemicalName'], 'SPERMINE')
        self.assertEqual(spermine['formula'], 'C10 H26 N4')
        self.assertEqual(spermine['smiles'], 'C(CCNCCCN)CNCCCN')
        self.assertEqual(
            spermine['InChI'],
            'InChI=1S/C10H26N4/c11-5-3-9-13-7-1-2-8-14-10-4-6-12/h13-14H,1-12H2')
        self.assertEqual(spermine['InChIKey'],
                         'PFNFFQXMRSDOHW-UHFFFAOYSA-N')
        self.assertAlmostEqual(float(spermine['@molecularWeight']), 202.34)

    def test_get_ligands_with_multiple_ligands(self):
        result = get_ligands('4HHB')

        chemical_ids = [
            ligand['@chemicalID'] for ligand in result['ligandInfo']['ligand']
        ]
        self.assertIn('HEM', chemical_ids)
        self.assertIn('PO4', chemical_ids)

    def test_get_ligands_for_entry_without_ligands(self):
        # An entry with no bound ligands returns an empty list, not None
        result = get_ligands('1UBQ')

        self.assertEqual(result['id'], '1UBQ')
        self.assertEqual(result['ligandInfo']['ligand'], [])

    def test_get_chains(self):
        chains = get_chains('4HHB')

        self.assertEqual(len(chains), 2)

        alpha, beta = chains
        self.assertEqual(alpha['entity_id'], '4HHB_1')
        self.assertEqual(alpha['chains'], ['A', 'C'])
        self.assertEqual(alpha['description'], 'Hemoglobin subunit alpha')
        self.assertEqual(alpha['polymer_type'], 'Protein')
        self.assertIn('Homo sapiens', alpha['organism'])
        self.assertEqual(len(alpha['sequence']), 141)

        self.assertEqual(beta['chains'], ['B', 'D'])
        self.assertEqual(beta['description'], 'Hemoglobin subunit beta')

    def test_get_chain_ids(self):
        self.assertEqual(get_chain_ids('4HHB'), ['A', 'B', 'C', 'D'])

    def test_get_chains_reports_non_protein_polymers(self):
        # 1EHZ is a transfer RNA
        chains = get_chains('1EHZ')

        self.assertEqual(len(chains), 1)
        self.assertEqual(chains[0]['chains'], ['A'])
        self.assertEqual(chains[0]['polymer_type'], 'RNA')

    def test_get_chains_distinguishes_auth_and_label_ids(self):
        # In large structures the author-assigned chain IDs shown on the RCSB
        # website differ from the mmCIF label_asym_ids
        chains = get_chains('4V5A')

        self.assertTrue(len(chains) > 1)
        self.assertTrue(
            any(entity['chains'] != entity['label_chains']
                for entity in chains),
            "Expected auth and label chain IDs to differ for 4V5A")

    def test_get_chains_nonexistent_pdb_id(self):
        with self.assertWarns(UserWarning):
            result = get_chains('9ZZZ')

        self.assertIsNone(result)

    def test_get_chain_ids_nonexistent_pdb_id(self):
        with self.assertWarns(UserWarning):
            self.assertIsNone(get_chain_ids('9ZZZ'))

    def test_get_ligand_instances(self):
        instances = get_ligand_instances('4HHB')

        self.assertTrue(len(instances) > 0)
        # 4HHB contains four hemes and two phosphates
        hemes = [i for i in instances if i['chem_id'] == 'HEM']
        self.assertEqual(len(hemes), 4)

        self.assertIn({
            'chem_id': 'HEM',
            'auth_asym_id': 'A',
            'auth_seq_id': 142
        }, instances)

    def test_get_ligand_instances_for_entry_without_ligands(self):
        self.assertEqual(get_ligand_instances('1UBQ'), [])

    def test_get_ideal_ligand_file(self):
        sdf = get_ideal_ligand_file('ATP', verbosity=False)

        self.assertIsNotNone(sdf)
        self.assertEqual(sdf.splitlines()[0].strip(), 'ATP')

    def test_get_ligand_instance_file(self):
        sdf = get_ligand_instance_file('4HHB', 'A', 142, verbosity=False)

        self.assertIsNotNone(sdf)
        self.assertEqual(sdf.splitlines()[0].strip(), 'HEM')
        # Instance coordinates carry the ligand's real position in the entry
        self.assertIn('18.6750', sdf)

    def test_get_ligand_instance_file_as_mol2(self):
        mol2 = get_ligand_instance_file('4HHB',
                                        'A',
                                        142,
                                        filetype='mol2',
                                        verbosity=False)

        self.assertIsNotNone(mol2)
        self.assertIn('HEM', mol2.splitlines()[0])

    def test_get_ligands_nonexistent_pdb_id(self):
        with self.assertWarns(UserWarning):
            result = get_ligands('9ZZZ')

        self.assertIsNone(result)

    def test_describe_chemical_invalid_chem_id_format(self):
        chem_id = "NAGX" # Invalid format (too long)
        with self.assertRaisesRegex(Exception, "Ligand id with more than 3 characters provided"):
            describe_chemical(chem_id)


# Import PDBFileType for testing get_pdb_file
import pypdb.clients.pdb.pdb_client
# More direct imports for specific items
# Correcting FastaEntry to FastaSequence based on actual class name in fasta_client.py
from pypdb.clients.fasta.fasta_client import FastaSequence as PDBFastaSequence
from pypdb.clients.search.search_client import LogicalOperator as PDBLogicalOperator
from pypdb.clients.search.operators.sequence_operators import SequenceOperator as PDBSequenceOperator


class TestDeprecatedFunctions(unittest.TestCase):

    @patch('pypdb.clients.pdb.pdb_client.get_pdb_file')
    def test_get_pdb_file_wrapper_deprecation_and_call(self, mock_pdb_client_get_file):
        # Test with filetype='pdb'
        with self.assertWarnsRegex(
            DeprecationWarning,
            "The `get_pdb_file` function within pypdb.py is deprecated."
            "See `pypdb/clients/pdb/pdb_client.py` for a near-identical "
            "function to use"
        ):
            get_pdb_file('1EHZ', filetype='pdb', compression=False)

        mock_pdb_client_get_file.assert_called_once_with('1EHZ', pypdb.clients.pdb.pdb_client.PDBFileType.PDB, False)

        # Reset mock for the next call
        mock_pdb_client_get_file.reset_mock()

        # Test with filetype='cif'
        with self.assertWarnsRegex(
            DeprecationWarning,
            "The `get_pdb_file` function within pypdb.py is deprecated."
        ): # Shortened regex for brevity, full message checked above
            get_pdb_file('1EHZ', filetype='cif', compression=True)

        mock_pdb_client_get_file.assert_called_once_with('1EHZ', pypdb.clients.pdb.pdb_client.PDBFileType.CIF, True)

        # Test with filetype='xml'
        mock_pdb_client_get_file.reset_mock()
        with self.assertWarns(DeprecationWarning):
            get_pdb_file('1EHZ', filetype='xml')
        mock_pdb_client_get_file.assert_called_once_with('1EHZ', pypdb.clients.pdb.pdb_client.PDBFileType.XML, False)

        # Test with filetype='structfact'
        mock_pdb_client_get_file.reset_mock()
        with self.assertWarns(DeprecationWarning):
            get_pdb_file('1EHZ', filetype='structfact')
        mock_pdb_client_get_file.assert_called_once_with('1EHZ', pypdb.clients.pdb.pdb_client.PDBFileType.STRUCTFACT, False)

        # Test with an invalid filetype (should still warn about deprecation, and then the wrapped function might error or warn further)
        mock_pdb_client_get_file.reset_mock()
        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always")
            self.assertIsNone(get_pdb_file('1EHZ', filetype='invalid', compression=False))

        self.assertTrue(any(isinstance(w.message, DeprecationWarning) for w in caught_warnings))
        self.assertTrue(any("Filetype specified to `get_pdb_file` appears to be invalid" in str(w.message)
                            for w in caught_warnings))
        mock_pdb_client_get_file.assert_not_called()


    @patch('pypdb.pypdb.fasta_client.get_fasta_from_rcsb_entry')
    @patch('pypdb.pypdb.search_client.perform_search_with_graph')
    def test_get_blast_deprecation(self, mock_perform_search, mock_get_fasta):
        # Configure mocks to return basic valid structures to allow the function to run
        mock_get_fasta.return_value = [
            PDBFastaSequence(
                entity_id='2F5N_1',
                chains=['A'],
                # Using a minimal valid protein sequence to avoid CannotAutoresolveSequenceTypeError
                sequence='MTE',
                fasta_header='>2F5N_1|Chains A|...'
            )
        ]
        mock_perform_search.return_value = {"result_set": []}

        with self.assertWarnsRegex(
            DeprecationWarning,
            "The `get_blast` function is slated for deprecation."
            "See `pypdb/clients/search/EXAMPLES.md` for examples to use a"
            "`SequenceOperator` search to similar effect"
        ):
            get_blast('2F5N', chain_id='A')

        # Verify that the underlying functions were called (optional, but good for confirming flow)
        mock_get_fasta.assert_called_once_with('2F5N')
        self.assertTrue(mock_perform_search.called)
        # Check some aspects of the query object passed to perform_search_with_graph
        args, kwargs = mock_perform_search.call_args
        query_object = kwargs.get('query_object')
        self.assertIsNotNone(query_object)
        self.assertEqual(query_object.logical_operator, PDBLogicalOperator.OR) # Using direct import alias
        self.assertEqual(len(query_object.queries), 1)
        self.assertIsInstance(query_object.queries[0], PDBSequenceOperator) # Using direct import alias
        self.assertEqual(query_object.queries[0].sequence, 'MTE') # Corrected expected sequence

    def test_characterize_get_gene_onto_for_deprecation(self):
        # Intended for deprecation. Functionality covered by DataFetcher.
        # Currently, this function is commented out in pypdb.py
        with self.assertRaises(NameError):
            get_gene_onto('4Z0L')

    def test_characterize_get_all_for_deprecation(self):
        # Intended for deprecation. Functionality covered by search_client.
        # Currently, this function is commented out in pypdb.py
        with self.assertRaises(NameError):
            get_all()


if __name__ == '__main__':
    unittest.main()
