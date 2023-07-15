import unittest

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
        self.assertTrue(type(found_pdbs[0]) == str)

        # an error page would be a longer string
        self.assertTrue(len(found_pdbs[0]) < 10)

    def test_pubmed(self):
        found_pdbs = Query(27499440, "PubmedIdQuery").search()
        self.assertTrue(len(found_pdbs) > 0)

    def test_treeentity(self):
        found_pdbs = Query('6239', 'TreeEntityQuery').search()
        self.assertTrue(len(found_pdbs) > 0)

    def test_exptype(self):
        found_pdbs = Query('SOLID-STATE NMR', 'ExpTypeQuery').search()
        self.assertTrue(len(found_pdbs) > 0)

    def test_structure(self):
        found_pdbs = Query('2E8D', 'structure').search()
        self.assertTrue(len(found_pdbs) > 0)

    def test_advancedauthor(self):
        found_pdbs = Query('Perutz, M.F.', 'AdvancedAuthorQuery').search()
        self.assertTrue(len(found_pdbs) > 0)

    def test_organism(self):
        found_pdbs = Query('Dictyostelium', 'OrganismQuery').search()
        self.assertTrue(len(found_pdbs) > 0)

    

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


if __name__ == '__main__':
    unittest.main()
