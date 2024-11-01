from addons import *

@ctest_labeler("dft;scf;cart")
def test_dft_grad_meta():
    ctest_runner(__file__)

