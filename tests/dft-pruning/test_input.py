from addons import *

@ctest_labeler("dft;scf")
def test_dft_pruning():
    ctest_runner(__file__)

