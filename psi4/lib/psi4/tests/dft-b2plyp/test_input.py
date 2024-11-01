from addons import *

@ctest_labeler("dft;scf")
def test_dft_b2plyp():
    ctest_runner(__file__)

