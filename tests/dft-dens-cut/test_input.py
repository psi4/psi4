from addons import *

@ctest_labeler("quick;dft;scf;cart")
def test_dft_dens_cut():
    ctest_runner(__file__)

