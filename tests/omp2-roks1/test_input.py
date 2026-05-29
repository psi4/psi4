from addons import *

@ctest_labeler("omp;dft")
def test_omp2_roks1():
    ctest_runner(__file__)
