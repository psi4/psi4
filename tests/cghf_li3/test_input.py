from addons import *

@uusing("einsums")
@ctest_labeler("quick;cghf")
def test_cghf_li3():
    ctest_runner(__file__, extra_infiles=["li3_cocc_guess.dat"])

