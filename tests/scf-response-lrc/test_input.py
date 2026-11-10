from addons import *

@ctest_labeler("quick;scf;properties;dft")
def test_scf_response_lrc():
    ctest_runner(__file__)
