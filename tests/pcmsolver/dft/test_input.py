from pathlib import Path
from addons import uusing


@uusing("pcmsolver")
def test_pcmsolver_dft():
    from qcengine.util import execute

    success, output = execute(["psi4", Path(__file__).parent / "input.dat"])

    if not success:
        print(output["stdout"])
    assert success
