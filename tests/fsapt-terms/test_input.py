from addons import *

@ctest_labeler("quick;sapt;cart;fsapt")
def test_fsapt_terms():
    ctest_runner(__file__, [
        "../../psi4/share/psi4/fsapt/fsapt.py",
    ])

