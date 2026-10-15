from addons import *

@ctest_labeler("quick;scf;soscf")
@orbital_optimizer_combinations
def test_scf_soscf_start(oopkg, soopkg):
    """SOSCF_START=1 reaches the occupation-restart path in the OpenTrustRegion combinations.

    With the default SOSCF_START the first-order iterations settle the aufbau assignment
    before second-order engages and the path is unreachable; this is the case that exercises
    it for real rather than by injection.
    """
    ctest_runner(__file__, setenv=orbital_optimizer_setenv(oopkg, soopkg))
