from addons import *

@ctest_labeler("scf;extern;findif;quick")
@integral_package_combinations
def test_extern7(intpkg):
    ctest_runner(__file__, setenv=integral_package_setenv(intpkg))

