from addons import *

@uusing("libefp")
@ctest_labeler("cart")
def test_libefp_qchem_efp_sp():
    ctest_runner(__file__)

