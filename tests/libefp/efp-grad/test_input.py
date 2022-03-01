from pathlib import Path
from addons import uusing


@uusing("libefp")
def test_libefp_efp_grad():
    from qcengine.util import execute

    success, output = execute(["psi4", Path(__file__).parent / "input.dat"])

    if not success:
        print(output["stdout"])
    assert success
