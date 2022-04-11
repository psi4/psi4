from addons import *

@ctest_labeler("nbody;gradient")
def test_nbody_convergence():
    ctest_runner(__file__)

