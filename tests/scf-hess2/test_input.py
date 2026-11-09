from addons import *

@ctest_labeler("quick;scf;freq;cart;d2ints")
@orbital_optimizer_combinations
@integral_package_combinations
def test_scf_hess2(oopkg, soopkg, intpkg):
    setenv = (orbital_optimizer_setenv(oopkg, soopkg) or []) + (integral_package_setenv(intpkg) or []) or None

    ctest_runner(__file__, setenv=setenv)

