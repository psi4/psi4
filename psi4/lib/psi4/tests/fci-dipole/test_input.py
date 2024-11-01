from addons import *

@ctest_labeler("quick;ci;properties;cart")
def test_fci_dipole():
    ctest_runner(__file__)

