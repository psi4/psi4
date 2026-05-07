from addons import *

@uusing("qcmanybody")
@ctest_labeler("cart")
def test_nbody_he_cluster():
    ctest_runner(__file__)
