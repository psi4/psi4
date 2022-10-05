from addons import *

@ctest_labeler("omp;gradient")
def test_omp2p5_grad1():
    ctest_runner(__file__)

