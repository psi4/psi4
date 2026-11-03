from addons import *

@ctest_labeler("quick;scf;cart")
@integral_package_combinations
def test_scf1(intpkg):
    ctest_runner(__file__, setenv=integral_package_setenv(intpkg))

