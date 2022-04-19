from addons import *

@ctest_labeler("shorttests;scf")
def test_soscf_ref():
    ctest_runner(__file__)

