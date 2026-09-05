from addons import *

@uusing("qcmanybody")
@ctest_labeler("dft;scf")
def test_dft_custom_dhdf():
    ctest_runner(__file__)
