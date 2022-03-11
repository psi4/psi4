from addons import *

@uusing("dftd3")
@ctest_labeler("quick;smoke;cart")
def test_dftd3_energy(tmp_path):
    ctest_runner(__file__)
