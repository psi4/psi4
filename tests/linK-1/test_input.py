from addons import *

@ctest_labeler("quick;scf;direct-scf")
@first_order_optimizer_combinations
def test_linK_1(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
