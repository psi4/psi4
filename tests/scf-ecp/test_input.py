from addons import *

@uusing("ecpint")
@ctest_labeler("scf;ecp;cart;smoke;quick")
def test_scf_ecp():
    ctest_runner(__file__)

