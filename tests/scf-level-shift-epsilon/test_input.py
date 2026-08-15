from addons import *

@ctest_labeler("quick;scf;dft")
def test_scf_level_shift_epsilon():
    ctest_runner(__file__)
