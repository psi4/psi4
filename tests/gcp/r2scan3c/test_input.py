from addons import *

@uusing("mctc-gcp")
@uusing("dftd4_350")
@ctest_labeler("quick")
def test_gcp_r2scan3c():
    ctest_runner(__file__)

