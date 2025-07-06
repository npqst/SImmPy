import unittest

import Bio


class TestSymmetryMates(unittest.TestCase):
    def test_get_symmetry_mates(self):
        import stcrpy

        tcr = stcrpy.fetch_TCRs("6ULR")
        assert len(tcr[0].get_MHC()) == 1
