import json
import os
import subprocess
import sys
from typing import List

from qcelemental.models import AtomicInput, OptimizationInput

from qcengine import cli, get_molecule, util
from qcengine.testing import using


def run_qcengine_cli(args: List[str], stdin: str = None) -> str:
    """
    Runs qcengine via its CLI

    This method was chosen over the more sophisticated util.execute in order to get pytest-cov to work.
    For the same reason, qcengine is invoked as "python -m qcengine.cli" rather than "qcengine".

    Parameters
    ----------
    args: List[str]
        List of CLI arguments.
    stdin: Optional[str]
        Standard input for the process.

    Returns
    -------
    str
        QCEngine CLI standard output.
    """
    if stdin is not None:
        stdin = stdin.encode("utf-8")

    return subprocess.check_output([sys.executable, "-m", "qcengine"] + args, input=stdin)


def test_no_args():
    """Test for qcengine with no arguments"""
    try:
        run_qcengine_cli([])
    except subprocess.CalledProcessError as e:
        assert e.returncode == 1


def test_info():
    """Test for qcengine info"""
    outputs = []
    for arg in cli.info_choices:
        output = run_qcengine_cli(["info", arg])
        if arg not in {"all", "config"}:  # output of config changes call-to-call depending e.g. on mem available
            outputs.append(output)

    default_output = run_qcengine_cli(["info"])
    for output in outputs:
        assert output in default_output


@using("psi4")
def test_run_psi4(tmp_path):
    """Tests qcengine run with psi4 and JSON input"""

    def check_result(stdout):
        output = json.loads(stdout)
        assert output["provenance"]["creator"].lower() == "psi4"
        assert output["success"] is True

    inp = AtomicInput(molecule=get_molecule("hydrogen"), driver="energy", model={"method": "hf", "basis": "6-31G"})

    args = ["run", "psi4", inp.json()]
    check_result(run_qcengine_cli(args))

    args = ["run", "psi4", os.path.join(tmp_path, "input.json")]
    with util.disk_files({"input.json": inp.json()}, {}, cwd=tmp_path):
        check_result(run_qcengine_cli(args))

    args = ["run", "psi4", "-"]
    check_result(run_qcengine_cli(args, stdin=inp.json()))


@using("geometric")
@using("psi4")
def test_run_procedure(tmp_path):
    """Tests qcengine run-procedure with geometric, psi4, and JSON input"""

    def check_result(stdout):
        output = json.loads(stdout)
        assert output["provenance"]["creator"].lower() == "geometric"
        assert output["success"] is True

    inp = {
        "keywords": {"coordsys": "tric", "maxiter": 100, "program": "psi4"},
        "input_specification": {"driver": "gradient", "model": {"method": "HF", "basis": "sto-3g"}, "keywords": {}},
        "initial_molecule": get_molecule("hydrogen"),
    }
    inp = OptimizationInput(**inp)

    args = ["run-procedure", "geometric", inp.json()]
    check_result(run_qcengine_cli(args))

    args = ["run-procedure", "geometric", os.path.join(tmp_path, "input.json")]
    with util.disk_files({"input.json": inp.json()}, {}, cwd=tmp_path):
        check_result(run_qcengine_cli(args))

    args = ["run-procedure", "geometric", inp.json()]
    check_result(run_qcengine_cli(args, stdin=inp.json()))
