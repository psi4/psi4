from addons import *

@ctest_labeler("quick;smoke;scf;properties;noc1")
def test_scf_property():
    ctest_runner(__file__, extra_infiles=["grid.dat"])

