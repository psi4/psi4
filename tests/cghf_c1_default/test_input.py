from addons import *

@uusing("einsums")
@ctest_labeler("smoke;quick;cghf")
def test_cghf_c1_default():
    ctest_runner(__file__)

