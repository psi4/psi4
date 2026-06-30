import pytest
from addons import *

@uusing("qcmanybody")
@pytest.mark.parametrize("distributed", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="snowflake", marks=using("qcfractal_next")),
])
@ctest_labeler("gradient")
def test_nbody_cp_gradient(distributed):
    setenv = ["_PSI4_USE_QCF"] if distributed else None
    ctest_runner(__file__, setenv=setenv)
