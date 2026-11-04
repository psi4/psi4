from addons import *

@uusing("qcmanybody")
@ctest_labeler("hessian")
def test_nbody_hessian():
    ctest_runner(__file__)
