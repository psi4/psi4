from addons import *

@uusing("adcc")
@ctest_labeler("quick;smoke")
def test_adcc_h2o_adc1_quick():
    ctest_runner(__file__)

