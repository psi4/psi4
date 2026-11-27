from addons import *

@ctest_labeler("misc;cart")
def test_pubchem1():
    with requires_pubchem():
        ctest_runner(__file__)

