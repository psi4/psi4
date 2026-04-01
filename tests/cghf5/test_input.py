from addons import *

@uusing("einsums")
@ctest_labeler("quick;cghf")
def test_cghf5():
    ctest_runner(__file__, extra_infiles=["density.dat"])

