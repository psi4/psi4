from addons import *

@ctest_labeler("quick;scf;properties")
def test_scf_response_uhf_open_shell():
    ctest_runner(__file__)
