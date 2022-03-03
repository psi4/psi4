from addons import *

@uusing("adcc")
@ctest_labeler("quick")
def test_adcc_h2o_cvs_adc2x_triplets():
    ctest_runner(__file__)

