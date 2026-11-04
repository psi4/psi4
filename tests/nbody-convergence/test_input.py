from addons import *

@uusing("qcmanybody")
@ctest_labeler("gradient")
def test_nbody_convergence():
    ctest_runner(__file__)
