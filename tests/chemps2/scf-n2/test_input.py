from addons import *

@uusing("chemps2")
@ctest_labeler("quick;smoke;cart")
def test_chemps2_scf_n2():
    ctest_runner(__file__)

