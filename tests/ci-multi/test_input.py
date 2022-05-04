from addons import *

@ctest_labeler("quick;misc;noc1")
def test_ci_multi():
    ctest_runner(__file__)

