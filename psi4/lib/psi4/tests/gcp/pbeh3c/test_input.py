from addons import *

@uusing("gcp")
@uusing("dftd3")
@ctest_labeler("cart")
def test_gcp_pbeh3c():
    ctest_runner(__file__)

