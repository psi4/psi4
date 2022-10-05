from addons import *

@ctest_labeler("nbody;extern")
def test_nbody_multi_level():
    ctest_runner(__file__)

