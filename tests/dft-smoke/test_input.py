from addons import *

@ctest_labeler("quick;dft;scf")
def test_dft_smoke():
    ctest_runner(__file__)

