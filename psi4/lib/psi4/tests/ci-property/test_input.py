from addons import *

@ctest_labeler("quick;ci;cas;properties;cart;noc1")
def test_ci_property():
    ctest_runner(__file__, extra_infiles=["grid.dat"])

