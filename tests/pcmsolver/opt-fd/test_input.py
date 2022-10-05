from addons import *

@uusing("pcmsolver")
@ctest_labeler("quick;smoke;scf;opt;cart")
def test_pcmsolver_opt_fd():
    ctest_runner(__file__)

