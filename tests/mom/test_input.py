from addons import *

@ctest_labeler("scf;mom;misc")
@orbital_optimizer_combinations
def test_mom(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
