from addons import *

@ctest_labeler("quick;scf;properties")
def test_scf_response_lindep():
    ctest_runner(__file__)
