from addons import *

@ctest_labeler("quicktests;scf;dft")
def test_soscf_dft():
    ctest_runner(__file__)

