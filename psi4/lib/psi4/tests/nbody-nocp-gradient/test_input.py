from addons import *

@ctest_labeler("nbody;gradient")
def test_nbody_nocp_gradient():
    ctest_runner(__file__)

