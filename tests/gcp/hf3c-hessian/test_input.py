from addons import *

@uusing("gcp")
@uusing("dftd3")
@ctest_labeler("hessian")
def test_gcp_hf3c_hessian():
    ctest_runner(__file__)

