from addons import *

@ctest_labeler("quick;df;dfmp2;gradient")
@integral_package_combinations
def test_dfmp2_grad1(intpkg):
    ctest_runner(__file__, setenv=integral_package_setenv(intpkg))

