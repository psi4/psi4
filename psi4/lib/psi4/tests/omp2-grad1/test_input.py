from addons import *

@ctest_labeler("quick;omp;gradient")
def test_omp2_grad1():
    ctest_runner(__file__)

