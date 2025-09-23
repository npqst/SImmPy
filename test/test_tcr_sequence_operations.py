import unittest


class TestTCRSequenceOperations(unittest.TestCase):

    def test_get_germlines(self):
        import stcrpy

        tcr = stcrpy.fetch_TCRs("8gvb")[0]

        germline_info = tcr.get_germline_assignments()

        mhc_info = tcr.get_MHC_allele_assignments()

        print(germline_info)
        print(mhc_info)
