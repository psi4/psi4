from addons import *

@uusing("pcmsolver")
@ctest_labeler("quick;smoke;tdscf")
def test_pcmsolver_tdscf():
    ctest_runner(__file__)

