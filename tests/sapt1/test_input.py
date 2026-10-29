from addons import *

@ctest_labeler("quick;smoke;sapt;cart")
@integral_package_combinations
def test_sapt1(intpkg):
    ctest_runner(__file__, setenv=integral_package_setenv(intpkg))

