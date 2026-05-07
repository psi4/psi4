from addons import *

@uusing("dftd3")
@uusing("qcmanybody")
@ctest_labeler("cart;gradient")
def test_dftd3_nbody_cp_gradient():
    ctest_runner(__file__)
