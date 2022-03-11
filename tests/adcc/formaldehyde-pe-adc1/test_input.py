from addons import *

@uusing("adcc")
@uusing("cppe")
@ctest_labeler("quick")
def test_adcc_formaldehyde_pe_adc1():
    ctest_runner(__file__, ["fa_6w.pot"])

