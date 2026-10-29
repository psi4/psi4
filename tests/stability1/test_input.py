from addons import *

@ctest_labeler("quick;stability;cart;noc1")
@orbital_optimizer_combinations
def test_stability1(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
