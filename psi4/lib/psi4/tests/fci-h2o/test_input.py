from addons import *

@ctest_labeler("ci;cart")
def test_fci_h2o():
    ctest_runner(__file__)

