from addons import *

@ctest_labeler("quick;scf")
def test_scf_upcast_custom_basis():
    ctest_runner(__file__)

