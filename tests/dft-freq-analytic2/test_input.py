from addons import *

@ctest_labeler("dft;scf;quick;findif")
def test_dft_freq_analytic2():
    ctest_runner(__file__)

