from addons import *

@uusing("adcc")
@uusing("cppe")
@ctest_labeler("")
def test_adcc_formaldehyde_pe_adc2():
    ctest_runner(__file__, ["fa_6w.pot"])

