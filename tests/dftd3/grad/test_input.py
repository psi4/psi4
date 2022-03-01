from pathlib import Path
from addons import uusing


# edit any marks
@uusing("dftd3")
def test_dftd3_grad():  # edit test name
    from qcengine.util import execute

    success, output = execute(["psi4", Path(__file__).parent / "input.dat"])

    if not success:
        print(output["stdout"])
    assert success
