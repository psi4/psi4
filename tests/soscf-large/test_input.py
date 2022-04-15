from addons import *

@ctest_labeler("quick;scf;cart")
def test_soscf_large():
    ctest_runner(__file__)

