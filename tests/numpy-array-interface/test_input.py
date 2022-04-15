from addons import *

@ctest_labeler("quick;pywrap;numpy;cart")
def test_numpy_array_interface():
    ctest_runner(__file__)

