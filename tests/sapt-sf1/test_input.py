import pytest
from addons import *

@ctest_labeler("quick;sapt")
@pytest.mark.parametrize("oopkg", [
    pytest.param(False, id="internal"),
    pytest.param("ooo", id="openorbitaloptimizer", marks=using("ooo")),
    pytest.param("otr", id="opentrustregion", marks=using("otr")),
])
def test_sapt_sf1(oopkg):
    setenv = {"ooo": ["_PSI4_USE_OOPKG"], "otr": ["_PSI4_USE_OTRPKG"]}.get(oopkg)
    ctest_runner(__file__, setenv=setenv)
