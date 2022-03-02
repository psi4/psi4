from addons import *

@uusing("dftd3")
@ctest_labeler("quick;smoke;cart")
def test_dftd3_energy():
    ctest_runner(__file__)
