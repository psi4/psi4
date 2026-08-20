import pytest
from addons import *

@uusing("qcmanybody")
@pytest.mark.parametrize("distributed", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="snowflake", marks=using("qcfractal_next")),
])
@ctest_labeler("cart")
def test_nbody_he_4b(distributed):
    setenv = ["_PSI4_USE_QCF"] if distributed else None
    ctest_runner(__file__, setenv=setenv)
