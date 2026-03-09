import pytest
from addons import *

@ctest_labeler("nbody;extern")
@pytest.mark.parametrize("distributed", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="snowflake", marks=[*using("qcfractal_next"), pytest.mark.medlong]),
])
def test_nbody_multi_level(distributed, monkeypatch):
    setenv = ["_PSI4_USE_QCF"] if distributed else None
    monkeypatch.setenv("QCMANYBODY_EMBEDDING_CHARGES", "1")
    ctest_runner(__file__, setenv=setenv)

