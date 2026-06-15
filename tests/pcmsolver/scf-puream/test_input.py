from addons import *


@uusing("pcmsolver")
@ctest_labeler("quick;scf;cart")
def test_pcmsolver_scf_puream():
    ctest_runner(__file__)
