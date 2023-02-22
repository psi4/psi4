from addons import *

@ctest_labeler("quick;scf;freq;cart;d2ints")
def test_scf_hess4():
    ctest_runner(__file__)

