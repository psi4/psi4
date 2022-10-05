from addons import *

@ctest_labeler("dft;scf;quick")
def test_dft_freq_analytic():
    ctest_runner(__file__)

