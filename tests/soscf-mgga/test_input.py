from addons import *

@ctest_labeler("scf;dft")
def test_soscf_mgga():
    ctest_runner(__file__)
