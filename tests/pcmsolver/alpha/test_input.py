from addons import *

@uusing("pcmsolver")
@ctest_labeler("quick;smoke;alpha")
def test_pcmsolver_alpha():
    ctest_runner(__file__)

