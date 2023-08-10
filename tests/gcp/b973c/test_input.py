from addons import *

@uusing("gcp")
@uusing("dftd3")
@ctest_labeler("quick")
def test_gcp_b973c():
    ctest_runner(__file__)

