from addons import *

@uusing("libefp")
@ctest_labeler("quick;scf;cart")
def test_libefp_qchem_qmefp_puream_sp():
    ctest_runner(__file__)

