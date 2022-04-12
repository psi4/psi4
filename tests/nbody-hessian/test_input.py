from addons import *

@ctest_labeler("nbody;hessian")
def test_nbody_hessian():
    ctest_runner(__file__)

