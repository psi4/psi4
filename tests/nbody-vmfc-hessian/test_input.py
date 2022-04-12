from addons import *

@ctest_labeler("nbody;hessian")
def test_nbody_vmfc_hessian():
    ctest_runner(__file__)

