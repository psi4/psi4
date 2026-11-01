from addons import *

@uusing("einsums")
@ctest_labeler("quick;cghf;scf")
def test_cghf6():
    ctest_runner(__file__, extra_binary_infiles=["density.dat"])

