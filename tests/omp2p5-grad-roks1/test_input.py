from addons import *

@ctest_labeler("omp;roks;gradient")
def test_omp2p5_grad_roks1():
    ctest_runner(__file__)
