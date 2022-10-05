from addons import *

@uusing("mrcc")
@ctest_labeler("cc;opt;freq")
def test_mrcc_optfreq():
    ctest_runner(__file__)

