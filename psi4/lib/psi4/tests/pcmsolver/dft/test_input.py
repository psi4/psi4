from addons import *

@uusing("pcmsolver")
@ctest_labeler("scf;dft;cart")
def test_pcmsolver_dft():
    ctest_runner(__file__)

