from addons import *

@uusing("dftd3")
@ctest_labeler("quick;smoke;sapt;cart;fsapt")
def test_fsapt_d():
    ctest_runner(__file__)

