from addons import *

@ctest_labeler("quick;dft;scf;cart")
@orbital_optimizer_combinations
def test_dft_dens_cut(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
