from addons import *

@ctest_labeler("shorttests;dft;scf;cart")
def test_dft_grad2():
    ctest_runner(__file__)

