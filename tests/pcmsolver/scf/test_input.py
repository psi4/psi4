from addons import *

@uusing("pcmsolver")
@ctest_labeler("quick;smoke;scf;cart")
def test_pcmsolver_scf():
    ctest_runner(__file__)

