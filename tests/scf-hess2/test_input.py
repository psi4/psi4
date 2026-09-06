from addons import *

@ctest_labeler("quick;scf;freq;cart;d2ints")
@orbital_optimizer_combinations
def test_scf_hess2(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
