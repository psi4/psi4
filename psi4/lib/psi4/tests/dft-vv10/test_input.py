from addons import *

@ctest_labeler("quick;dft;scf;cart")
def test_dft_vv10():
    ctest_runner(__file__)

