from addons import *

@ctest_labeler("scf;noc1")
def test_scf_coverage():
    ctest_runner(__file__)

