from addons import *

@ctest_labeler("quick;smoke;casscf;noc1")
def test_casscf_sp():
    ctest_runner(__file__)

