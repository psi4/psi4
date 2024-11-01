from addons import *

@uusing("gcp")
@uusing("dftd3")
@ctest_labeler("dft")
def test_gcp_dft_custom_hybrid():
    ctest_runner(__file__)

