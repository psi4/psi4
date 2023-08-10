from addons import *

@uusing("mctc-gcp")
@uusing("s-dftd3")
@ctest_labeler("quick")
def test_gcp_b973c():
    ctest_runner(__file__)

