from addons import *

@ctest_labeler("dft;scf")
def test_dft_custom_mgga():
    ctest_runner(__file__)

