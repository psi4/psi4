import pytest
from addons import *

@ctest_labeler("dft;scf")
@pytest.mark.parametrize("oopkg", [
    pytest.param(False, id="internal"),
    pytest.param("OOO", id="openorbitaloptimizer", marks=using("ooo")),
    pytest.param("OTR", id="opentrustregion", marks=using("otr")),
])
def test_dft_psivar(oopkg):
    setenv = [f"_PSI4_USE_{oopkg}PKG"]
    ctest_runner(__file__, setenv=setenv)
