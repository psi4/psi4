from pathlib import Path
import psi4
from addons import *
import pytest

@uusing("dftd4")
@ctest_labeler("quick;smoke;sapt;cart;fsapt")
def test_fsapt_d4():
    fsaptpy_installed = (Path(psi4.core.get_datadir()) / "fsapt" / "fsapt.py").resolve()

    ctest_runner(__file__, extra_infiles=[
        fsaptpy_installed,
    ])
