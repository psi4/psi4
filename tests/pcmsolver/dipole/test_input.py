from addons import *

@uusing("pcmsolver")
@ctest_labeler("quick;scf;properties;dipole")
def test_pcmsolver_dipole():
    ctest_runner(__file__)

