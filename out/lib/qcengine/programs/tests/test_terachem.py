import pytest
import qcelemental as qcel
from qcelemental.testing import compare_recursive

import qcengine as qcng
from qcengine.testing import qcengine_records, using

# Prep globals
terachem_info = qcengine_records("terachem")


@pytest.mark.parametrize("test_case", terachem_info.list_test_cases())
def test_terachem_output_parser(test_case):
    # Get output file data
    data = terachem_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    output = qcng.get_program("terachem", check=False).parse_output(data, inp).dict()
    output_ref = qcel.models.AtomicResult.parse_raw(data["output.json"]).dict()

    # Forgiving molecule since it is now sparse
    assert compare_recursive(output_ref, output, forgive={"stdout", "provenance", "molecule"})


@pytest.mark.parametrize("test_case", terachem_info.list_test_cases())
def test_terachem_input_formatter(test_case):
    # Get input file data
    data = terachem_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    # TODO add actual comparison of generated input file
    input_file = qcng.get_program("terachem", check=False).build_input(inp, qcng.get_config())
    assert input_file.keys() >= {"commands", "infiles"}


@using("terachem")
@pytest.mark.parametrize("test_case", terachem_info.list_test_cases())
def test_terachem_executor(test_case):
    # Get input file data
    data = terachem_info.get_test_data(test_case)
    inp = qcel.models.AtomicInput.parse_raw(data["input.json"])

    # Run Terachem
    result = qcng.compute(inp, "terachem")
    # result = qcng.get_program('terachem').compute(inp, qcng.get_config())
    assert result.success is True

    # Get output file data
    output_ref = qcel.models.AtomicResult.parse_raw(data["output.json"])

    atol = 1e-6
    if result.driver == "gradient":
        atol = 1e-3
    assert compare_recursive(output_ref.return_result, result.return_result, atol=atol)
