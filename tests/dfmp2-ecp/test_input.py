import pytest
from addons import *

@uusing("ecpint")
@ctest_labeler("quick;df;dfmp2;ecp;nbody")
@pytest.mark.parametrize("oopkg", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="openorbitaloptimizer", marks=using("ooo")),
])
def test_dfmp2_ecp(oopkg):
    setenv = ["_PSI4_USE_OOPKG"] if oopkg else None

    ctest_runner(__file__, setenv=setenv)
