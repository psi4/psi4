from addons import *

@ctest_labeler("quick;scf;direct-scf")
@first_order_optimizer_combinations
def test_linK_1(oopkg, soopkg):
    # only vary the first-order package as 2nd-order can't do semi-numerical exch (nonsym)
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
