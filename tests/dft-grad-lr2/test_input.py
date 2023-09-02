from addons import *

@ctest_labeler("dft;scf;findif")
def test_dft_grad_lr2():
    ctest_runner(__file__)

