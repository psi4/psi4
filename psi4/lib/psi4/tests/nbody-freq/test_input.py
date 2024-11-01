from addons import *

@ctest_labeler("nbody;freq;hessian;gradient;opt")
def test_nbody_freq():
    ctest_runner(__file__)

