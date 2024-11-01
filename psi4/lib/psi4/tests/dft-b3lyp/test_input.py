from addons import *

@ctest_labeler("dft;scf")
def test_dft_b3lyp():
    ctest_runner(__file__)

