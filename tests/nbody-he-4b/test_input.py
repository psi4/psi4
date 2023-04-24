import pytest
from addons import *

@pytest.mark.parametrize("distributed", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="snowflake", marks=using("qcfractal")),
])
@ctest_labeler("nbody;cart")
def test_nbody_he_4b(distributed):
    setenv = ["_PSI4_USE_QCF"] if distributed else None
    ctest_runner(__file__, setenv=setenv)

