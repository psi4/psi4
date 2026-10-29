from addons import *

@uusing("ecpint")
@ctest_labeler("scf;ecp;cart;smoke;quick")
@orbital_optimizer_combinations
def test_scf_ecp(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
