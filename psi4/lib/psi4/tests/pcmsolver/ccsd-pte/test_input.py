from addons import *

@uusing("pcmsolver")
@ctest_labeler("quick;cc;pte;cart")
def test_pcmsolver_ccsd_pte():
    ctest_runner(__file__)

