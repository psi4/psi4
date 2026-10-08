from addons import *

@ctest_labeler("opt;cart")
@orbital_optimizer_combinations
def test_opt2(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
