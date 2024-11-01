from addons import *

@ctest_labeler("dft;scf;cart")
def test_dft_grad1():
    ctest_runner(__file__)

