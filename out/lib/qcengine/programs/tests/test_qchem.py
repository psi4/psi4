import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_recursive, compare_values

import qcengine as qcng
from qcengine.testing import qcengine_records, using

qchem_info = qcengine_records("qchem")
qchem_logonly_info = qcengine_records("qchem_logonly")

qchem_forgive = [
    "root.molecule.provenance.version",
    "root.molecule.provenance.routine",
    "root.provenance.version",
    "root.provenance.routine",
    "root.provenance.creator",
    "root.extras",
]


@pytest.mark.parametrize("test_case", qchem_info.list_test_cases())
def test_qchem_output_parser(test_case):

    # Get output file data
    data = qchem_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    outfiles = qcel.util.deserialize(data["outfiles.msgpack"], "msgpack-ext")
    output = qcng.get_program("qchem", check=False).parse_output(outfiles, inp).dict()
    output.pop("provenance", None)

    output_ref = qcel.models.AtomicResult.parse_raw(data["output.json"]).dict()
    output_ref.pop("provenance", None)
    output_ref.pop("extras", None)
    output.pop("extras", None)

    check, message = compare_recursive(output_ref, output, return_message=True)
    assert check, message


@pytest.mark.parametrize("test_case", qchem_info.list_test_cases())
def test_qchem_input_formatter(test_case):

    # Get input file data
    data = qchem_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program("qchem", check=False).build_input(inp, qcng.get_config())
    assert input_file.keys() >= {"commands", "infiles"}


@pytest.mark.parametrize("test_case", qchem_info.list_test_cases())
def test_qchem_input_formatter_template(test_case):

    # Get input file data
    data = qchem_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program("qchem", check=False).build_input(inp, qcng.get_config(), template="Test template")
    assert input_file.keys() >= {"commands", "infiles"}


@using("qchem")
@pytest.mark.parametrize("test_case", qchem_info.list_test_cases())
def test_qchem_executor(test_case):
    # Get input file data
    data = qchem_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    # Run qchem
    result = qcng.compute(inp, "qchem")
    assert result.success is True

    # Get output file data
    output_ref = qcel.models.AtomicResult.parse_raw(data["output.json"])

    atol = 1e-6
    assert compare_recursive(output_ref.return_result, result.return_result, atol=atol)


@using("qchem")
def test_qchem_orientation():

    mol = qcel.models.Molecule.from_data(
        """
        He 0.0  0.7  0.7
        He 0.0 -0.7 -0.7
        """
    )

    # Compare with rotation
    inp = {"molecule": mol, "driver": "gradient", "model": {"method": "HF", "basis": "6-31g"}}
    ret = qcng.compute(inp, "qchem", raise_error=True)
    assert compare_values(np.linalg.norm(ret.return_result, axis=0), [0, 0, 0.00791539])

    # Compare without rotations
    mol_noorient = mol.copy(update={"fix_orientation": True})
    inp = {"molecule": mol_noorient, "driver": "gradient", "model": {"method": "HF", "basis": "6-31g"}}
    ret = qcng.compute(inp, "qchem", raise_error=True)

    assert compare_values(np.linalg.norm(ret.return_result, axis=0), [0, 0.00559696541, 0.00559696541])


@pytest.mark.parametrize("test_case", qchem_logonly_info.list_test_cases())
def test_qchem_logfile_parser(test_case):

    # Get output file data
    data = qchem_logonly_info.get_test_data(test_case)
    outfiles = {"dispatch.out": data["qchem.out"]}
    with pytest.warns(Warning):
        output = qcng.get_program("qchem", check=False).parse_logfile(outfiles).dict()
    output["stdout"] = None

    output_ref = qcel.models.AtomicResult.parse_raw(data["output.json"]).dict()
    for key in list(output["provenance"].keys()):
        if key not in output_ref["provenance"]:
            output["provenance"].pop(key)

    # Modify ref to trim down total data as a molecule is now sparse
    output_ref["molecule"] = {k: v for k, v in output_ref["molecule"].items() if k in output["molecule"]}

    check, message = compare_recursive(output_ref, output, return_message=True, forgive=qchem_forgive)
    assert check, message


@pytest.mark.parametrize("test_case", qchem_info.list_test_cases())
def test_qchem_logfile_parser_qcscr(test_case):

    # Get output file data
    data = qchem_info.get_test_data(test_case)
    outfiles = qcel.util.deserialize(data["outfiles.msgpack"], "msgpack-ext")

    with pytest.warns(Warning):
        output = qcng.get_program("qchem", check=False).parse_logfile(outfiles).dict()
    output["stdout"] = None

    output_ref = qcel.models.AtomicResult.parse_raw(data["output.json"]).dict()
    for key in list(output["provenance"].keys()):
        if key not in output_ref["provenance"]:
            output["provenance"].pop(key)

    output_ref["stdout"] = None

    # Modify ref to trim down total data as a molecule is now sparse
    output_ref["molecule"] = {k: v for k, v in output_ref["molecule"].items() if k in output["molecule"]}

    output_ref["model"]["method"] = output_ref["model"]["method"].lower()
    check, message = compare_recursive(output_ref, output, return_message=True, forgive=qchem_forgive)
    assert check, message
