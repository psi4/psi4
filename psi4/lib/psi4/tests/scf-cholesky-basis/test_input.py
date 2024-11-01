from addons import *

@ctest_labeler("quick;scf;noc1")
def test_scf_cholesky_basis():
    ctest_runner(__file__)

