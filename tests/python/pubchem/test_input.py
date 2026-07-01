import pytest
from addons import *

@pytest.mark.xfail(reason="Pubchem is sometimes too busy")
@ctest_labeler("minitests;cart")
def test_python_pubchem():
    ctest_runner(__file__)

