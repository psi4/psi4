import pytest

import numpy as np
import psi4
import json

from .utils import *

pytestmark = pytest.mark.quick

data_blob = {
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
} # yapf: disable

def test_qcschema_energy():
    inp = data_blob.copy()
    ret = psi4.schema_wrapper.run_qcschema(inp)

    expected_return_result = -76.22831410222477
    assert compare_integers(True, ret.success, "Computation Status")
    assert compare_values(expected_return_result, ret.return_result, 5, "Return Energy")
    assert compare_values(expected_return_result, ret.properties.return_energy, 5, "Properties Return Energy")
    assert compare_values(-122.44529682915068, ret.properties.scf_one_electron_energy, 5, "SCF One-Electron Energy")

    # Check Array data
    assert compare_integers(True, isinstance(ret.extras["qcvars"]["MAYER_INDICES"], np.ndarray),
                            "Extras: Mayer Indices is Array")
    assert compare_integers(True, ret.provenance.routine == "psi4.schema_runner.run_qcschema", "Provenance: Routine")


def test_qcschema_gradient():
    inp = data_blob.copy()
    inp["driver"] = "gradient"
    ret = psi4.schema_wrapper.run_qcschema(inp)

    bench_grad = np.array([[0.00000000e+00, 0.00000000e+00, -3.14003710e-02],
                           [0.00000000e+00, -3.02856483e-02, 1.57001855e-02],
                           [0.00000000e+00, 3.02856483e-02, 1.57001855e-02]])

    assert compare_integers(True, ret.success, "Computation Status")
    assert compare_arrays(bench_grad, ret.return_result, 5, "Return Gradient")
    assert compare_values(-76.22831410222477, ret.properties.return_energy, 5, "Properties Return Energy")
    assert compare_values(-122.44529682915068, ret.properties.scf_one_electron_energy, 5, "SCF One-Electron Energy")


def test_qcschema_keyword_error():
    inp = data_blob.copy()
    inp["keywords"] = {"unicorn": "bad option"}
    ret = psi4.schema_wrapper.run_qcschema(inp)

    assert compare_integers(False, ret.success, "Computation Status")
    assert compare_integers(True, "unicorn" in ret.error.error_message, "Error Message")
