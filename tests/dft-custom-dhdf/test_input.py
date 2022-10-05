from addons import *

@ctest_labeler("dft;scf;nbody")
def test_dft_custom_dhdf():
    ctest_runner(__file__)

