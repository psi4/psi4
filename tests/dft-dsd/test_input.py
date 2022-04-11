from addons import *

@ctest_labeler("dft;scf;cart;nbody")
def hide_test_dft_dsd():
    ctest_runner(__file__)

