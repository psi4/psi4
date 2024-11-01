from addons import *

@ctest_labeler("ci;properties;noc1")
def test_fci_tdm():
    ctest_runner(__file__)

