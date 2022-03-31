from addons import *

@uusing("libefp")
@ctest_labeler("quick;smoke;scf;cart")
def test_libefp_qchem_qmefp_sp():
    ctest_runner(__file__)

