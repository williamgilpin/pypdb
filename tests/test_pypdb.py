import unittest

## Import from local directory
import sys
sys.path.insert(0, '../pypdb')
from pypdb import *

# TODO(ejwilliams): Write generic logic, to execute `test_*.py` files
# within the pypdb directory (removing need for sys.path hack)


class TestSearchFunctions(unittest.TestCase):

    # Returns True if it successfully connects to the
    # protein data bank
    def test_querying(self):
        found_pdbs = Query('actin network').search()
        self.assertTrue(len(found_pdbs) > 0)
        self.assertTrue(type(found_pdbs[0]) == str)

        # an error page would be a longer string
        self.assertTrue(len(found_pdbs[0]) < 10)

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
