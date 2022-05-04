from addons import *

@ctest_labeler("dft;scf")
def test_dft_grad_disk():
    ctest_runner(__file__)

