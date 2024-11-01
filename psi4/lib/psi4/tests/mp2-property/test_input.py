from addons import *

@ctest_labeler("quick;mp2;dfmp2;properties;cart")
def test_mp2_property():
    ctest_runner(__file__, extra_infiles=["grid.dat"])

