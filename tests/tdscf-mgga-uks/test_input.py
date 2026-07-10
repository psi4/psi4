from addons import *

@ctest_labeler("tdscf;dft")
def test_tdscf_mgga_uks():
    ctest_runner(__file__)
