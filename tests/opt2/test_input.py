import pytest
from addons import *

@ctest_labeler("opt;cart")
@pytest.mark.parametrize("oopkg", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="openorbitaloptimizer", marks=using("ooo")),
])
def test_opt2(oopkg):
    setenv = ["_PSI4_USE_OOPKG"] if oopkg else None

    ctest_runner(__file__, setenv=setenv)

