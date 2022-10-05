from addons import *

@uusing("chemps2")
@ctest_labeler("quick;cart")
def test_chemps2_caspt2_small():
    ctest_runner(__file__)

