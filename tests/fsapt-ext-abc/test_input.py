from pathlib import Path
import psi4
from addons import *

@ctest_labeler("quick;sapt;cart;fsapt;extern")
def test_fsapt_ext_abc():
    fsaptpy_installed = (Path(psi4.core.get_datadir()) / "fsapt" / "fsapt.py").resolve()

    ctest_runner(__file__, [
        fsaptpy_installed,
    ])

