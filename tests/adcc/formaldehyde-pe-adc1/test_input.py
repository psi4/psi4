from addons import *

@uusing("adcc")
@uusing("cppe")
@ctest_labeler("quick")
def hide_test_adcc_formaldehyde_pe_adc1():
    ctest_runner(__file__)

