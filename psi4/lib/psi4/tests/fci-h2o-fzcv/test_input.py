from addons import *

@ctest_labeler("fci;cart;noc1")
def test_fci_h2o_fzcv():
    ctest_runner(__file__)

