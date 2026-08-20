from addons import *

@uusing("qcmanybody")
@ctest_labeler("quick;smoke;df;dfmp2;cart")
def test_dfmp2_1():
    ctest_runner(__file__)
