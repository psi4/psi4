from addons import *

@uusing("einsums")
@ctest_labeler("quick;cghf")
def test_cghf4():
    ctest_runner(__file__, extra_binary_infiles=["density.dat"])

