from addons import *
import pytest

@uusing("einsums")
@ctest_labeler("quick;cghf")
@pytest.mark.xfail
def test_cghf2():
    ctest_runner(__file__)

