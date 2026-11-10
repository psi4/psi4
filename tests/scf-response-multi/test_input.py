from addons import *

@ctest_labeler("quick;scf;properties")
def test_scf_response_multi():
    ctest_runner(__file__)
