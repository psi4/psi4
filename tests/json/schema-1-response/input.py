#! test QC_JSON Schema for response properties

import pytest
import psi4
import numpy as np
import json
import os
from distutils import dir_util


@pytest.fixture
def datadir(tmpdir, request):
    """
    from: https://stackoverflow.com/a/29631801
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    """
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, str(tmpdir))

    return tmpdir


def test_response(datadir):
    # Generate JSON input
    json_data = {
        "schema_name": "qcschema_input",
        "schema_version": 1,
        "molecule": {
            "geometry": [
                0.0, 0.0, -0.1294769411935893, 0.0, -1.494187339479985, 1.0274465079245698, 0.0, 1.494187339479985,
                1.0274465079245698
            ],
            "symbols": ["O", "H", "H"],
            "fix_com":
            True,
            "fix_orientation":
            True
        },
        "driver": "properties",
        "model": {
            "method": "CC2",
            "basis": "6-31G",
        },
        "keywords": {
            "scf_type": "df",
            "mp2_type": "df",
            "e_convergence": 9,
            "omega": [355, 439, 'nm'],
            "gauge": "velocity",
            "function_kwargs": {
                "properties": ["dipole", "polarizability", "rotation", "roa_tensor"]
            }
        }
    }

    # Load expected output (dipole & quadrupole in au)
    reference_file = datadir.join(f"output.json.ref")
    with open(reference_file) as f:
        expected_response = json.load(f)

    # Convert lists to np arrays
    arrays = [
        'CC2 DIPOLE', 'CC2 QUADRUPOLE', 'CC2 DIPOLE POLARIZABILITY TENSOR @ 355NM',
        'CC2 DIPOLE POLARIZABILITY TENSOR @ 439NM', 'CC2 OPTICAL ROTATION TENSOR (MVG) @ 355NM',
        'CC2 OPTICAL ROTATION TENSOR (MVG) @ 439NM', 'CC2 OPTICAL ROTATION TENSOR (VEL) @ 0NM',
        'CC2 OPTICAL ROTATION TENSOR (VEL) @ 355NM', 'CC2 OPTICAL ROTATION TENSOR (VEL) @ 439NM',
        'CC2 QUADRUPOLE POLARIZABILITY TENSOR COMPONENT 0 @ 355NM',
        'CC2 QUADRUPOLE POLARIZABILITY TENSOR COMPONENT 0 @ 439NM',
        'CC2 QUADRUPOLE POLARIZABILITY TENSOR COMPONENT 1 @ 355NM',
        'CC2 QUADRUPOLE POLARIZABILITY TENSOR COMPONENT 1 @ 439NM',
        'CC2 QUADRUPOLE POLARIZABILITY TENSOR COMPONENT 2 @ 355NM',
        'CC2 QUADRUPOLE POLARIZABILITY TENSOR COMPONENT 2 @ 439NM'
    ]
    for a in arrays:
        expected_response[a] = np.asarray(expected_response[a])

    json_ret = psi4.json_wrapper.run_qcschema(json_data).dict()

    psi4.compare_integers(True, json_ret["success"], "JSON Success")  #TEST
    psi4.compare_strings("qcschema_output", json_ret["schema_name"], "Schema Name")  #TEST
    for k in expected_response.keys():  #TEST
        psi4.compare_values(expected_response[k], json_ret["extras"]["qcvars"][k], 5, "Result: " + k.upper())  #TEST
