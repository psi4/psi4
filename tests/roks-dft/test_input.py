from addons import *

@ctest_labeler("dft;scf;roks")
def test_roks_dft():
    ctest_runner(__file__)
