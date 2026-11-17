from addons import *

@uusing("einsums")
@ctest_labeler("f12;df;quick")
@integral_package_combinations
def test_dfmp2f12_1_energy(intpkg):
    ctest_runner(__file__, setenv=integral_package_setenv(intpkg))
