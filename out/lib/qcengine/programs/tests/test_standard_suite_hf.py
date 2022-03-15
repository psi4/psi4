import pytest
import qcelemental as qcel
from qcelemental.testing import compare_values

import qcengine as qcng
from qcengine.testing import using


@pytest.fixture
def h2o():
    smol = """
 # R=0.958 A=104.5
 H                  0.000000000000     1.431430901356     0.984293362719
 O                  0.000000000000     0.000000000000    -0.124038860300
 H                  0.000000000000    -1.431430901356     0.984293362719
 units au
"""
    return qcel.models.Molecule.from_data(smol)


@pytest.fixture
def nh2():
    smol = """
 # R=1.008 #A=105.0
 0 2
 N   0.000000000000000   0.000000000000000  -0.145912918634892
 H   0.000000000000000  -1.511214298139000   1.013682596946108
 H   0.000000000000000   1.511214298139000   1.013682596946108
 units au
 symmetry c1
"""
    return qcel.models.Molecule.from_data(smol)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param("cfour", "aug-pvdz", {"scf_conv": 12}, marks=using("cfour")),
        pytest.param("cfour", "aug-pvdz", {}, marks=using("cfour")),
        pytest.param(
            "qcore",
            "aug-cc-pVDZ",
            {"coulomb_method": "direct_4idx", "exchange_method": "direct_4idx"},
            marks=using("qcore"),
        ),
        pytest.param("gamess", "accd", {"contrl__ispher": 1}, marks=using("gamess")),
        pytest.param("molpro", "aug-cc-pvdz", {}, marks=using("molpro")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True}, marks=using("nwchem")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "qc_module": "tce"}, marks=using("nwchem")),
        pytest.param("psi4", "aug-cc-pvdz", {"scf_type": "direct"}, marks=using("psi4")),
        pytest.param("qchem", "aug-cc-pvdz", {}, marks=using("qchem")),
        pytest.param("turbomole", "aug-cc-pVDZ", {}, marks=using("turbomole")),
        pytest.param("terachem_pbs", "aug-cc-pvdz", {}, marks=using("terachem_pbs")),
    ],
)
def test_sp_hf_rhf(program, basis, keywords, h2o):
    """cfour/sp-rhf-hf/input.dat
    #! single point HF/adz on water

    """
    resi = {"molecule": h2o, "driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords}

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    scf_tot = -76.0413815332

    atol = 1.0e-6
    assert compare_values(scf_tot, res["return_result"], atol=atol)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param(
            "cfour",
            "aug-pvdz",
            {"reference": "uhf", "occupation": [[3, 1, 1, 0], [3, 0, 1, 0]], "scf_conv": 12},
            marks=using("cfour"),
        ),
        pytest.param("cfour", "aug-pvdz", {"reference": "uhf"}, marks=using("cfour")),
        pytest.param(
            "qcore",
            "aug-cc-pVDZ",
            {"ansatz": "u", "coulomb_method": "direct_4idx", "exchange_method": "direct_4idx"},
            marks=using("qcore"),
        ),
        pytest.param("gamess", "accd", {"contrl__ispher": 1, "contrl__scftyp": "uhf"}, marks=using("gamess")),
        pytest.param("molpro", "aug-cc-pvdz", {"reference": "unrestricted"}, marks=using("molpro")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "scf__uhf": True}, marks=using("nwchem")),
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"basis__spherical": True, "qc_module": "tce", "scf__uhf": True},
            marks=using("nwchem"),
        ),
        pytest.param("psi4", "aug-cc-pvdz", {"reference": "uhf", "scf_type": "direct"}, marks=using("psi4")),
        pytest.param("qchem", "aug-cc-pvdz", {}, marks=using("qchem")),
        pytest.param("turbomole", "aug-cc-pVDZ", {}, marks=using("turbomole")),
    ],
)
def test_sp_hf_uhf(program, basis, keywords, nh2):
    resi = {"molecule": nh2, "driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords}

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    scf_tot = -55.57513805253009

    atol = 1.0e-6
    assert compare_values(scf_tot, res["return_result"], atol=atol)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param(
            "cfour",
            "aug-pvdz",
            {"reference": "rohf", "occupation": [[3, 1, 1, 0], [3, 0, 1, 0]], "scf_conv": 12},
            marks=using("cfour"),
        ),
        pytest.param("cfour", "aug-pvdz", {"reference": "rohf"}, marks=using("cfour")),
        pytest.param("gamess", "accd", {"contrl__ispher": 1, "contrl__scftyp": "rohf"}, marks=using("gamess")),
        pytest.param("molpro", "aug-cc-pvdz", {}, marks=using("molpro")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "scf__rohf": True}, marks=using("nwchem")),
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"basis__spherical": True, "qc_module": "tce", "scf__rohf": True},
            marks=using("nwchem"),
        ),
        pytest.param("psi4", "aug-cc-pvdz", {"reference": "rohf", "scf_type": "direct"}, marks=using("psi4")),
        pytest.param("qchem", "aug-cc-pvdz", {"UNRESTRICTED": False}, marks=using("qchem")),
    ],
)
def test_sp_hf_rohf(program, basis, keywords, nh2):
    resi = {"molecule": nh2, "driver": "energy", "model": {"method": "hf", "basis": basis}, "keywords": keywords}

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    scf_tot = -55.570724348574

    atol = 1.0e-6
    assert compare_values(scf_tot, res["return_result"], atol=atol)
