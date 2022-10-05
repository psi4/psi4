from addons import *

@uusing("ecpint")
@ctest_labeler("quick;df;dfmp2;ecp;nbody")
def test_dfmp2_ecp():
    ctest_runner(__file__)

