import sys
import json
import pprint

import numpy as np
import pytest
import qcelemental as qcel
from qcengine.testing import schema_versions, from_v2, checkver_and_convert

import psi4

from utils import *

pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

_ispy314 = sys.version_info >= (3, 14)


@pytest.fixture(scope="function")
def result_data_fixture(request):
    atin_v1 = {
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
            "scf_properties": ["mayer_indices"],
        }
    }
    if from_v2(request.node.name):
        mol = atin_v1.pop("molecule")
        atin_v2 = {"molecule": mol, "specification": atin_v1}
        return atin_v2
    else:
        return atin_v1

def test_qcschema_energy(result_data_fixture, schema_versions, request):
    models, retver, _ = schema_versions

    if not (_ispy314 and "as_v2" not in request.node.name):
        # checkver_and_convert does a little QCSchema version casting & checking. not always on b/c requires models not avail for py314
        checkver_and_convert(result_data_fixture, request.node.name, "pre")
    # return_dict=T to permit v1 w/py314
    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture, return_dict=True, return_version=retver)
    if not (_ispy314 and "as_v2" not in request.node.name):
        checkver_and_convert(ret, request.node.name, "post")

    expected_return_result = -76.22831410222477
    assert compare_integers(True, ret["success"], "Computation Status"), pp.pprint(ret)
    assert compare_values(expected_return_result, ret["return_result"], 5, "Return Energy")
    assert compare_values(expected_return_result, ret["properties"]["return_energy"], 5, "Properties Return Energy")
    assert compare_values(-122.44529682915068, ret["properties"]["scf_one_electron_energy"], 5, "SCF One-Electron Energy")

    # Check stdout
    assert compare_integers(True, 'Psi4: An Open' in ret["stdout"], "Stdout Header")
    assert compare_integers(True, 'beer' in ret["stdout"], "Stdout Beer")

    # Check Array data
    mayer = ret["extras"]["qcvars"]["MAYER INDICES"]
    assert compare_integers(True, isinstance(mayer, np.ndarray),
                            f"Extras: Mayer Indices is Array: {type(mayer)}")
    assert compare_integers(True, ret["provenance"]["routine"] == "psi4.schema_runner.run_qcschema", "Provenance: Routine")


def test_qcschema_gradient(result_data_fixture, schema_versions, request):
    _, retver, _ = schema_versions

    if sys.version_info >= (3, 14) and not ("as_v2" in request.node.name or "to_v2" in request.node.name):
        pytest.skip("Py314 and v1.AtomicResult object incompatible.")

    if psi4.core.get_option("scf", "orbital_optimizer_package") != "INTERNAL":
        result_data_fixture["keywords"].update({"e_convergence": 9, "d_convergence": 5e-9})

    if from_v2(request.node.name):
        result_data_fixture["specification"]["driver"] = "gradient"
    else:
        result_data_fixture["driver"] = "gradient"
    if not (_ispy314 and "as_v2" not in request.node.name):
        result_data_fixture = checkver_and_convert(result_data_fixture, request.node.name, "pre")
    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture, return_version=retver)
    if not (_ispy314 and "as_v2" not in request.node.name):
        ret = checkver_and_convert(ret , request.node.name, "post")

    bench_grad = np.array([[0.00000000e+00, 0.00000000e+00, -3.14003710e-02],
                           [0.00000000e+00, -3.02856483e-02, 1.57001855e-02],
                           [0.00000000e+00, 3.02856483e-02, 1.57001855e-02]])

    assert compare_integers(True, ret.success, "Computation Status")
    assert compare_arrays(bench_grad, ret.return_result, 5, "Return Gradient")
    assert compare_values(-76.22831410222477, ret.properties.return_energy, 5, "Properties Return Energy")
    assert compare_values(-122.44529682915068, ret.properties.scf_one_electron_energy, 5, "SCF One-Electron Energy")


def test_qcschema_keyword_error(result_data_fixture, schema_versions, request):
    _, retver, _ = schema_versions

    if from_v2(request.node.name):
        result_data_fixture["specification"]["keywords"] = {"unicorn": "bad option"}
    else:
        result_data_fixture["keywords"] = {"unicorn": "bad option"}

    if not (_ispy314 and "as_v2" not in request.node.name):
        result_data_fixture = checkver_and_convert(result_data_fixture, request.node.name, "pre")
    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture, return_dict=True, return_version=retver)
    if not (_ispy314 and "as_v2" not in request.node.name):
        ret = checkver_and_convert(ret , request.node.name, "post", cast_dict_as="FailedOperation", vercheck=False)

    assert compare_integers(False, ret["success"], "Computation Status")
    assert compare_integers(True, "UNICORN" in ret["error"]["error_message"], "Error Message")


@pytest.mark.parametrize("input_enc, input_fn, output_enc, output_fn, ", [
    ("json", "input.json", "json", None),
    ("msgpack-ext", "input.msgpack", "msgpack-ext", None),
    ("json", "input.json", "msgpack-ext", "output.msgpack"),
    ("msgpack-ext", "input.msgpack", "json", "output.json"),
    ("msgpack-ext", "input.msgpack", "msgpack-ext", "output.something"),
])
def test_qcschema_cli(input_enc, input_fn, output_enc, output_fn, result_data_fixture, schema_versions, request):
    models, retver, models_out = schema_versions

    if _ispy314 and "as_v2" not in request.node.name:
        # this test we can't work around models.v1 with return_dict=T but we want to test the CLI route,
        #   so have to publically admit _v1v2 exists. NO FUTURE GUARANTEE
        if "as_v1" in request.node.name:
            models, models_out = qcel.models._v1v2, qcel.models._v1v2
        elif "to_v1" in request.node.name:
            models_out = qcel.models._v1v2
        elif "to_v2" in request.node.name:
            models = qcel.models._v1v2
        else:
            models, models_out = qcel.models._v1v2, qcel.models._v1v2

    data = models.AtomicInput(**result_data_fixture)

    if (input_enc == "msgpack-ext") or (output_enc == "msgpack-ext"):
        try:
            import msgpack
        except ImportError:
            pytest.skip("Msgpack could not be found, skipping.")

    inputs = {input_fn: data.serialize(input_enc)}

    cmds = ["--qcschema", "--return-version", str(retver)]
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
    pp.pprint(ret)
    assert compare_integers(True, success, "Computation Status"), pp.pprint(ret)
    # command shows up in stdout for Windows
    if not sys.platform.startswith("win"):
        assert compare_integers(True, ret['stdout'] == '', "Empty stdout"), ret['stdout']

    try:
        parsed = True
        ret = models_out.AtomicResult.parse_raw(ret["outfiles"][output_fn], encoding=output_enc)
    except Exception as e:
        parsed = False
        print(e)

    assert compare_integers(True, parsed, "Result Model Parsed")
    assert compare_values(-76.22831410207938, ret.return_result, "Return")


def test_qcschema_wavefunction_basis(result_data_fixture, schema_versions, request):
    models, retver, models_out = schema_versions

    if from_v2(request.node.name):
        result_data_fixture["specification"]["protocols"] = {"wavefunction": "orbitals_and_eigenvalues"}
    else:
        result_data_fixture["protocols"] = {"wavefunction": "orbitals_and_eigenvalues"}

    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture, return_dict=True, return_version=retver)

    assert compare_integers(True, ret["success"], "Computation Status"), pp.pprint(ret)
    wfn = ret["wavefunction"]

    assert len(wfn["basis"]["center_data"]) == 2
    assert len(wfn["basis"]["atom_map"]) == 3
    assert wfn["basis"]["atom_map"][1] == wfn["basis"]["atom_map"][2]
    assert wfn["basis"]["nbf"] == 24
    assert wfn["restricted"]


def test_qcschema_wavefunction_scf_orbitals(result_data_fixture, schema_versions, request):
    models, retver, models_out = schema_versions

    if from_v2(request.node.name):
        result_data_fixture["specification"]["protocols"] = {"wavefunction": "orbitals_and_eigenvalues"}
    else:
        result_data_fixture["protocols"] = {"wavefunction": "orbitals_and_eigenvalues"}

    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture, return_dict=True, return_version=retver)

    assert compare_integers(True, ret["success"], "Computation Status"), pp.pprint(ret)
    wfn = ret["wavefunction"]

    expected_keys = {'basis', 'restricted', 'scf_orbitals_a', 'scf_eigenvalues_a', 'orbitals_a', 'eigenvalues_a'}
    assert wfn.keys() == expected_keys


def test_qcschema_wavefunction_scf_occupations_gs(result_data_fixture, schema_versions, request):
    models, retver, models_out = schema_versions

    if from_v2(request.node.name):
        result_data_fixture["specification"]["protocols"] = {"wavefunction": "all"}
        result_data_fixture["specification"]["keywords"]["docc"] = [3, 0, 1, 1] # ground state
    else:
        result_data_fixture["protocols"] = {"wavefunction": "all"}
        result_data_fixture["keywords"]["docc"] = [3, 0, 1, 1] # ground state

    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture, return_dict=True, return_version=retver)

    assert compare_integers(True, ret["success"], "Computation Status"), pp.pprint(ret)
    wfn = ret["wavefunction"]

    # correctness check

    ref_occupations_a = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,])
    ref_occupied_energies_a = np.array([-20.55785069,  -1.31618596,  -0.6770761,  -0.55872283,  -0.49037545])

    assert compare_arrays(ref_occupations_a, wfn["scf_occupations_a"], 6, "Orbital Occupations")
    assert compare_arrays(ref_occupied_energies_a, wfn["scf_eigenvalues_a"][wfn["scf_occupations_a"] == 1], 6, "Occupied Orbital Energies")

    # consistency check

    Ca = wfn["scf_orbitals_a"]
    Fa = wfn["scf_fock_a"]
    Da = wfn["scf_density_a"]
    ea = wfn["scf_eigenvalues_a"]
    na = wfn["scf_occupations_a"]
    Ca_occ = Ca[:, na == 1]

    assert compare_arrays(ea, (Ca.T @ Fa @ Ca).diagonal(), 10, "Orbital Consistency")
    assert compare_arrays(Da, Ca_occ @ Ca_occ.T, 10, "Occupied Orbital Consistency")


def test_qcschema_wavefunction_scf_occupations_es(result_data_fixture, schema_versions, request):
    models, retver, models_out = schema_versions

    if from_v2(request.node.name):
        result_data_fixture["specification"]["protocols"] = {"wavefunction": "all"}
        result_data_fixture["specification"]["keywords"]["docc"] = [2, 1, 1, 1] # excited state
    else:
        result_data_fixture["protocols"] = {"wavefunction": "all"}
        result_data_fixture["keywords"]["docc"] = [2, 1, 1, 1] # excited state

    ret = psi4.schema_wrapper.run_qcschema(result_data_fixture, return_dict=True, return_version=retver)

    assert compare_integers(True, ret["success"], "Computation Status"), pp.pprint(ret)
    wfn = ret["wavefunction"]

    # correctness check

    ref_occupations_a = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,])
    ref_occupied_energies_a = np.array([-20.86710716,  -1.38875044,  -0.77442913,  -0.6598582,    1.10374473])

    assert compare_arrays(ref_occupations_a, wfn["scf_occupations_a"], 6, "Orbital Occupations")
    assert compare_arrays(ref_occupied_energies_a, wfn["scf_eigenvalues_a"][wfn["scf_occupations_a"] == 1], 6, "Occupied Orbital Energies")

    # consistency check

    Ca = wfn["scf_orbitals_a"]
    Fa = wfn["scf_fock_a"]
    Da = wfn["scf_density_a"]
    ea = wfn["scf_eigenvalues_a"]
    na = wfn["scf_occupations_a"]
    Ca_occ = Ca[:, na == 1]

    assert compare_arrays(ea, (Ca.T @ Fa @ Ca).diagonal(), 10, "Orbital Consistency")
    assert compare_arrays(Da, Ca_occ @ Ca_occ.T, 10, "Occupied Orbital Consistency")

