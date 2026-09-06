from addons import *

@uusing("qcmanybody")
@ctest_labeler("quick;scf;direct-scf")
def test_cfmm_4():
    ctest_runner(__file__)
