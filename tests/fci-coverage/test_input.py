from addons import *

@ctest_labeler("ci;cart")
def test_fci_coverage():
    ctest_runner(__file__)

