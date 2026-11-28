from addons import *

@ctest_labeler("minitests;cart")
def test_python_pubchem():
    with xfail_on_pubchem_unavailable():
        ctest_runner(__file__)

