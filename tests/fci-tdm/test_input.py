from addons import *

@ctest_labeler("ci;properties;noc1")
@internal_and_ooootr_combinations
def test_fci_tdm(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))

