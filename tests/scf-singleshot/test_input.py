from addons import *

@ctest_labeler("scf;dft")
def test_scf_singleshot():
    ctest_runner(__file__)
