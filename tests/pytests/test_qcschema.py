import json

import numpy as np
import pytest
import qcelemental as qcel

import psi4

from .utils import *

pytestmark = pytest.mark.quick

@pytest.fixture(scope="function")
def result_data_fixture():
    return {
        "molecule": {
            "geometry": [
                0.0, 0.0, -0.1294769411935893,
                0.0, -1.494187339479985, 1.0274465079245698,
                0.0, 1.494187339479985, 1.0274465079245698
            ],
            "symbols": ["O", "H", "H"],
            "connectivity": [(0, 1, 1.0), (0, 2, 1.0)]
        },
        "driver": "energy",
        "model": {
            "method": "MP2",
            "basis": "cc-pVDZ"
        },
        "keywords": {
            "scf_type": "df",
            "mp2_type": "df",
            "scf_properties": ["mayer_indices"]
        }
    }

def test_qcschema_energy(result_data_fixture):
    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture)

    expected_return_result = -76.22831410222477
    assert compare_integers(True, ret.success, "Computation Status")
    assert compare_values(expected_return_result, ret.return_result, 5, "Return Energy")
    assert compare_values(expected_return_result, ret.properties.return_energy, 5, "Properties Return Energy")
    assert compare_values(-122.44529682915068, ret.properties.scf_one_electron_energy, 5, "SCF One-Electron Energy")

    # Check stdout
    assert compare_integers(True, 'Psi4: An Open' in ret.stdout, "Stdout Header")
    assert compare_integers(True, 'beer' in ret.stdout, "Stdout Beer")

    # Check Array data
    assert compare_integers(True, isinstance(ret.extras["qcvars"]["MAYER_INDICES"], np.ndarray),
                            "Extras: Mayer Indices is Array")
    assert compare_integers(True, ret.provenance.routine == "psi4.schema_runner.run_qcschema", "Provenance: Routine")


def test_qcschema_gradient(result_data_fixture):
    result_data_fixture["driver"] = "gradient"
    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture)

    bench_grad = np.array([[0.00000000e+00, 0.00000000e+00, -3.14003710e-02],
                           [0.00000000e+00, -3.02856483e-02, 1.57001855e-02],
                           [0.00000000e+00, 3.02856483e-02, 1.57001855e-02]])

    assert compare_integers(True, ret.success, "Computation Status")
    assert compare_arrays(bench_grad, ret.return_result, 5, "Return Gradient")
    assert compare_values(-76.22831410222477, ret.properties.return_energy, 5, "Properties Return Energy")
    assert compare_values(-122.44529682915068, ret.properties.scf_one_electron_energy, 5, "SCF One-Electron Energy")


def test_qcschema_keyword_error(result_data_fixture):
    result_data_fixture["keywords"] = {"unicorn": "bad option"}
    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture)

    assert compare_integers(False, ret.success, "Computation Status")
    assert compare_integers(True, "unicorn" in ret.error.error_message, "Error Message")


@pytest.mark.parametrize("input_enc, input_fn, output_enc, output_fn, ", [
    ("json", "input.json", "json", None),
    ("msgpack-ext", "input.msgpack", "msgpack-ext", None),
    ("json", "input.json", "msgpack-ext", "output.msgpack"),
    ("msgpack-ext", "input.msgpack", "json", "output.json"),
    ("msgpack-ext", "input.msgpack", "msgpack-ext", "output.something"),
])
def test_qcschema_cli(input_enc, input_fn, output_enc, output_fn, result_data_fixture):

    data = qcel.models.ResultInput(**result_data_fixture)

    if (input_enc == "msgpack-ext") or (output_enc == "msgpack-ext"):
        try:
            import msgpack
        except:
            pytest.skip("Msgpack could not be found, skipping.")

    inputs = {input_fn: data.serialize(input_enc)}

    cmds = ["--qcschema"]
    if output_fn:
        outfiles = [output_fn]
        cmds.extend(["-o", output_fn])
    else:
        outfiles = [input_fn]
        output_fn = input_fn

    as_binary = []
    if input_enc == "msgpack-ext":
        as_binary.append(input_fn)
    if output_enc == "msgpack-ext":
        as_binary.append(output_fn)

    success, ret = run_psi4_cli(inputs, outfiles, cmds, as_binary=as_binary)
    assert compare_integers(True, success, "Computation Status")
    assert compare_integers(True, ret['stdout'] == '', "Empty stdout")

    try:
        parsed = True
        ret = qcel.models.Result.parse_raw(ret["outfiles"][output_fn], encoding=output_enc)
    except Exception as e:
        parsed = False
        print(e)

    assert compare_integers(True, parsed, "Result Model Parsed")
    assert compare_values(-76.22831410207938, ret.return_result, "Return")

def test_qcschema_wavefunction_basis(result_data_fixture):
    result_data_fixture["protocols"] = {"wavefunction": "orbitals_and_eigenvalues"}
    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture)
    wfn = ret.wavefunction

    assert len(wfn.basis.center_data) == 2
    assert len(wfn.basis.atom_map) == 3
    assert wfn.basis.atom_map[1] == wfn.basis.atom_map[2]
    assert wfn.basis.nbf == 24
    assert wfn.restricted

def test_qcschema_wavefunction_scf_orbitals(result_data_fixture):
    result_data_fixture["protocols"] = {"wavefunction": "orbitals_and_eigenvalues"}
    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture)
    wfn = ret.wavefunction

    expected_keys = {'basis', 'restricted', 'scf_orbitals_a', 'scf_eigenvalues_a', 'orbitals_a', 'eigenvalues_a'}
    assert wfn.dict().keys() == expected_keys
