from addons import *

@ctest_labeler("dft;scf;mp2")
def test_dft_reference():
    ctest_runner(__file__)

