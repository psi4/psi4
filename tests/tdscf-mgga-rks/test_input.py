from addons import *

@ctest_labeler("tdscf;dft")
def test_tdscf_mgga_rks():
    ctest_runner(__file__)
