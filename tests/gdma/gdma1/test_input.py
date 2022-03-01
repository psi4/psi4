from pathlib import Path
from addons import uusing


@uusing("gdma")
def test_gdma_gdma1():
    from qcengine.util import execute

    success, output = execute(["psi4", Path(__file__).parent / "input.dat"])

    if not success:
        print(output["stdout"])
    assert success
