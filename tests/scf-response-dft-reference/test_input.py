from addons import *

@ctest_labeler("quick;scf;properties;dft")
def test_scf_response_dft_reference():
    ctest_runner(__file__)
