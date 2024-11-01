from addons import *

@ctest_labeler("quick;sapt;cart;fsapt;extern")
def test_fsapt_ext():
    ctest_runner(__file__)

