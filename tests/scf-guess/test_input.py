from addons import *

@ctest_labeler("quick;scf;noc1")
def test_scf_guess():
    ctest_runner(__file__)

