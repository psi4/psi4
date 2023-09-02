from addons import *

@ctest_labeler("quick;scf;freq;cart;findif;d2ints")
def test_scf_hess3():
    ctest_runner(__file__)

