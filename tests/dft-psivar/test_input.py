from addons import *

@ctest_labeler("dft;scf")
@orbital_optimizer_combinations
def test_dft_psivar(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
