from addons import *

@uusing("qcmanybody")
@ctest_labeler("gradient")
def test_nbody_vmfc_gradient():
    ctest_runner(__file__)
