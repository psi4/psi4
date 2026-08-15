from addons import *

@ctest_labeler("scf;mom;dft;misc")
def test_mom_roks():
    ctest_runner(__file__)
