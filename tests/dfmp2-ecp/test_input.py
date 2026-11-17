from addons import *

@uusing("ecpint")
@uusing("qcmanybody")
@ctest_labeler("quick;df;dfmp2;ecp")
@orbital_optimizer_combinations
def test_dfmp2_ecp(oopkg, soopkg):
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
