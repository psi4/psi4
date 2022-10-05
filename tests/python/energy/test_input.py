from addons import *

@ctest_labeler("quick;smoke")
def test_python_energy():
    ctest_runner(__file__)

