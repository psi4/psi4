from addons import *

@ctest_labeler("dft;scf")
def test_dft_ghost():
    ctest_runner(__file__)

