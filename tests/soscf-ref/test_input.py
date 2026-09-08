from addons import *

@ctest_labeler("shorttests;scf")
@second_order_optimizer_combinations
def test_soscf_ref(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
