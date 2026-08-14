import pytest
from addons import *

@ctest_labeler("dft;scf")
@pytest.mark.parametrize("cuest", [
    pytest.param(False, id="internal"),
    pytest.param(True,  id="cuest", marks=[*using("cuest"), *using("cuda_cc8")]),
])
def test_dft1(cuest):
    setenv = ["_PSI4_USE_CUEST"] if cuest else None

    ctest_runner(__file__, setenv=setenv)
