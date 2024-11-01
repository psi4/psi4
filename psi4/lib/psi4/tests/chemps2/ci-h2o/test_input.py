from addons import *

@uusing("chemps2")
@ctest_labeler("cart")
def test_chemps2_ci_h2o():
    ctest_runner(__file__)

