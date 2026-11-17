from addons import *

@ctest_labeler("quick;sapt")
@orbital_optimizer_combinations
def test_sapt_sf1(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
