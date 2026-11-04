import pytest
from addons import *

@pytest.mark.xfail(reason="Pubchem is sometimes too busy")
@ctest_labeler("misc;cart")
def test_pubchem1():
    ctest_runner(__file__)

