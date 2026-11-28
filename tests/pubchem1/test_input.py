from addons import *

@ctest_labeler("misc;cart")
def test_pubchem1():
    with xfail_on_pubchem_unavailable():
        ctest_runner(__file__)

