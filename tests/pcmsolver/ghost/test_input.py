from addons import *

@uusing("qcmanybody")
@uusing("pcmsolver")
@ctest_labeler("quick;smoke;ghosts")
def test_pcmsolver_ghost():
    ctest_runner(__file__)

