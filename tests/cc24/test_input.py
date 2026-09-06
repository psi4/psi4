from addons import *

@ctest_labeler("cc;cart;noc1;eom;findif")
@orbital_optimizer_combinations
def test_cc24(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
