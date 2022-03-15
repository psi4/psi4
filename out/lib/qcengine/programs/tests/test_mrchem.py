"""Tests for MRChem functionality"""
import numpy as np
import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using


@pytest.fixture
def h2o():
    return qcel.models.Molecule(
        geometry=[[0, 0, -0.1250], [-1.4375, 0, 1.0250], [1.4375, 0, 1.0250]],
        symbols=["O", "H", "H"],
        connectivity=[[0, 1, 1], [0, 2, 1]],
    )


@using("mrchem")
def test_energy(h2o):
    mr_kws = {
        "world_prec": 1.0e-3,
        "world_size": 6,
        "world_unit": "bohr",
    }

    inp = qcel.models.AtomicInput(
        molecule=h2o,
        driver="energy",
        model={
            "method": "BLYP",
        },
        keywords=mr_kws,
    )

    res = qcng.compute(inp, "mrchem", raise_error=True, return_dict=True)

    # Make sure the calculation completed successfully
    assert compare_values(-76.4546307, res["return_result"], atol=1e-3)
    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # Make sure the properties parsed correctly
    assert compare_values(-76.4546307, res["properties"]["return_energy"], atol=1e-3)
    assert res["properties"]["calcinfo_natom"] == 3
    assert res["properties"]["calcinfo_nalpha"] == 5
    assert res["properties"]["calcinfo_nbeta"] == 5
    assert res["properties"]["calcinfo_nmo"] == 10
    assert compare_values([-3.766420e-07, 0.0, 0.720473], res["properties"]["scf_dipole_moment"], atol=1e-3)


@using("mrchem")
def test_dipole(h2o):
    mr_kws = {
        "world_prec": 1.0e-3,
        "world_size": 6,
        "world_unit": "bohr",
    }

    inp = qcel.models.AtomicInput(
        molecule=h2o,
        driver="properties",
        model={
            "method": "BLYP",
        },
        keywords=mr_kws,
    )

    res = qcng.compute(inp, "mrchem", raise_error=True, return_dict=True)

    # Make sure the calculation completed successfully
    assert compare_values([-3.766420e-07, 0.0, 0.720473], res["return_result"]["dipole_moment"]["dip-1"], atol=1e-3)
    assert res["driver"] == "properties"
    assert "provenance" in res
    assert res["success"] is True

    # Make sure the properties parsed correctly
    assert compare_values(-76.4546307, res["properties"]["return_energy"], atol=1e-3)
    assert res["properties"]["calcinfo_natom"] == 3
    assert res["properties"]["calcinfo_nalpha"] == 5
    assert res["properties"]["calcinfo_nbeta"] == 5
    assert res["properties"]["calcinfo_nmo"] == 10
    assert compare_values([-3.766420e-07, 0.0, 0.720473], res["properties"]["scf_dipole_moment"], atol=1e-3)
