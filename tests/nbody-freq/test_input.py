from addons import *

@uusing("qcmanybody")
@ctest_labeler("freq;hessian;gradient;opt")
def test_nbody_freq():
    ctest_runner(__file__)
