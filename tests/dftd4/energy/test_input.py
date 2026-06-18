from addons import *

@uusing("dftd4")
@uusing("dftd3")
@ctest_labeler("quick;smoke;cart")
def test_dftd4_energy():
    ctest_runner(__file__)

