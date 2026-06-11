import pytest
from addons import *

@uusing("ecpint")
@ctest_labeler("scf;ecp;cart;smoke;quick")
@pytest.mark.parametrize("oopkg", [
    pytest.param(False, id="internal"),
    pytest.param("ooo", id="openorbitaloptimizer", marks=using("ooo")),
    pytest.param("otr", id="opentrustregion", marks=using("otr")),
])
def test_scf_ecp(oopkg):
    setenv = ["_PSI4_USE_OOPKG"] if oopkg == "ooo" else None
    setenv = ["_PSI4_USE_OTRPKG"] if oopkg == "otr" else None
    ctest_runner(__file__, setenv=setenv)
