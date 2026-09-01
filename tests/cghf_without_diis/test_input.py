from addons import *

@uusing("einsums")
@ctest_labeler("quick;cghf")
def test_cghf_without_diis():
    ctest_runner(__file__)

