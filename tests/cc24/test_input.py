import pytest
from addons import *

@ctest_labeler("cc;cart;noc1;eom;findif")
@pytest.mark.parametrize("oopkg", [
    pytest.param(False, id="internal"),
    pytest.param("ooo", id="openorbitaloptimizer", marks=using("ooo")),
    pytest.param("otr", id="opentrustregion", marks=using("otr")),
])
def test_cc24(oopkg):
    setenv = {"ooo": ["_PSI4_USE_OOPKG"], "otr": ["_PSI4_USE_OTRPKG"]}.get(oopkg)

    ctest_runner(__file__, setenv=setenv)

