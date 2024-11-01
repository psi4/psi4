from addons import *

@uusing("gcp")
@uusing("dftd3")
@ctest_labeler("hessian;d2ints")
def test_gcp_hf3c_hessian():
    ctest_runner(__file__)

