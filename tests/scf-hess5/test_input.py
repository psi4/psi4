from addons import *

@ctest_labeler("quick;scf;freq;cart")
def test_scf_hess5():
    ctest_runner(__file__)

