from addons import *

@ctest_labeler("quick;scf;dft")
def test_roks_soscf_incfock():
    ctest_runner(__file__)
