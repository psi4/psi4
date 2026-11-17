from addons import *

@ctest_labeler("quick;scf;noc1")
@integral_package_combinations
def test_scf_guess(intpkg):
    ctest_runner(__file__, setenv=integral_package_setenv(intpkg))

