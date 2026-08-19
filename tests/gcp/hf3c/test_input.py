from addons import *

@uusing("qcmanybody")
@uusing("gcp")
@uusing("dftd3")
@ctest_labeler("quick;smoke")
def test_gcp_hf3c():
    ctest_runner(__file__)

