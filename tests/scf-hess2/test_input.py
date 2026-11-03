import pytest
from addons import *

@ctest_labeler("quick;scf;freq;cart;d2ints")
@integral_package_combinations
@pytest.mark.parametrize("oopkg", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="openorbitaloptimizer", marks=using("ooo")),
])
def test_scf_hess2(oopkg, intpkg):
    setenv = (["_PSI4_USE_OOPKG"] if oopkg else []) + (integral_package_setenv(intpkg) or []) or None

    ctest_runner(__file__, setenv=setenv)

