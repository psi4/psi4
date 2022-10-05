from pathlib import Path
import psi4
from addons import *

@ctest_labeler("sapt;cart;fsapt")
def test_fsapt_allterms():
    fsaptpy_installed = (Path(psi4.core.get_datadir()) / "fsapt" / "fsapt.py").resolve()

    ctest_runner(__file__, [
        fsaptpy_installed,
    ])

