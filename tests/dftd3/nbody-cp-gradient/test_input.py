from addons import *

@uusing("dftd3")
@ctest_labeler("cart;nbody;gradient")
def test_dftd3_nbody_cp_gradient():
    ctest_runner(__file__)
