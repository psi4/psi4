import pytest
from addons import *

@ctest_labeler("nbody;hessian")

@pytest.mark.parametrize("distributed", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="snowflake", marks=using("qcfractal_next")),
])
def test_nbody_vmfc_hessian(distributed):
    setenv = ["_PSI4_USE_QCF"] if distributed else None
    ctest_runner(__file__, setenv=setenv)

