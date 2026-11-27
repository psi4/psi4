from addons import *

@ctest_labeler("quick;smoke;cc")
@integral_package_combinations
def test_cc1(intpkg):
    ctest_runner(__file__, setenv=integral_package_setenv(intpkg))

