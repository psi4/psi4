from addons import *

@uusing("gcp")
@uusing("dftd3")
@ctest_labeler("gradient")
def test_gcp_hf3c_gradient():
    ctest_runner(__file__)

