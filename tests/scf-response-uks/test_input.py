from addons import *

@ctest_labeler("quick;scf;properties")
def test_scf_response_uks():
    ctest_runner(__file__)
