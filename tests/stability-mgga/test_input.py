from addons import *

@ctest_labeler("stability;dft")
def test_stability_mgga():
    ctest_runner(__file__)
