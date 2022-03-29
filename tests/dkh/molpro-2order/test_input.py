from addons import *

@uusing("dkh")
@ctest_labeler("quick;smoke;scf")
def test_dkh_molpro_2order():
    ctest_runner(__file__)

