import pytest
from addons import *

@ctest_labeler("findif;quick;smoke")
@pytest.mark.parametrize("distributed", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="snowflake", marks=using("qcfractal")),
])
def test_ddd_deriv(distributed):
    setenv = ["_PSI4_USE_QCF"] if distributed else None

    ctest_runner(__file__, setenv=setenv)

