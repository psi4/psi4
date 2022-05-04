from addons import *

@ctest_labeler("dft;scf")
def test_dft_custom_gga():
    ctest_runner(__file__)

