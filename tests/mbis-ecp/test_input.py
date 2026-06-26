from addons import *

@uusing("ecpint")
@ctest_labeler("properties;quick;ecp")
def test_mbis_ecp():
    ctest_runner(__file__)
