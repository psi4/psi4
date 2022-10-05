from addons import *

@ctest_labeler("shorttests;scf;dft")
def test_soscf_dft():
    ctest_runner(__file__)

