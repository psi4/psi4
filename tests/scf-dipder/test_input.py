from addons import *

@ctest_labeler("quick;scf;d2ints")
def test_scf_dipder():
    ctest_runner(__file__)

