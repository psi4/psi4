from addons import *

@uusing("mrcc")
@ctest_labeler("quick;smoke;cc")
def test_mrcc_ccsdt():
    ctest_runner(__file__)

