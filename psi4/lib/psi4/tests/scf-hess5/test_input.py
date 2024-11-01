from addons import *

@ctest_labeler("quick;scf;freq;cart;findif;d2ints")
def test_scf_hess5():
    ctest_runner(__file__)

