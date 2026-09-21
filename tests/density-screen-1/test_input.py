from addons import *

@ctest_labeler("quick;scf")
@first_order_optimizer_combinations
def test_density_screen_1(oopkg, soopkg):
    # only vary the first-order package as 2nd-order can't work with incfock
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
