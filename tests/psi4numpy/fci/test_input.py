from addons import *

@ctest_labeler("ci;cart")
def test_psi4numpy_fci():
    ctest_runner(__file__)

