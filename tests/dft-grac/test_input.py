from addons import *

@ctest_labeler("shorttests;dft;scf")
def test_dft_grac():
    ctest_runner(__file__)

