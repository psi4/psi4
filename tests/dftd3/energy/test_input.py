from pathlib import Path
from addons import uusing, ctlabels


@uusing("dftd3")
@ctlabels("quick;smoke;cart")
def test_dftd3_energy():
    from qcengine.util import execute

    success, output = execute(["psi4", Path(__file__).parent / "input.dat"])

    if not success:
        print(output["stdout"])
    assert success
