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
"""
    return qcel.models.Molecule.from_data(smol)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param("cfour", "aug-pvdz", {"scf_conv": 12, "cc_conv": 12}, marks=using("cfour")),
        pytest.param("cfour", "aug-pvdz", {}, marks=using("cfour")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True}, marks=using("nwchem")),
        pytest.param("nwchem", "aug-cc-pvdz", {"basis__spherical": True, "qc_module": "tce"}, marks=using("nwchem")),
        pytest.param("psi4", "aug-cc-pvdz", {}, marks=using("psi4")),
        pytest.param("gamess", "accd", {"ccinp__ncore": 0, "contrl__ispher": 1}, marks=using("gamess")),
    ],
)
def test_sp_ccsd_t_rhf_full(program, basis, keywords, h2o):
    """cfour/sp-rhf-ccsd/input.dat
    #! single point CCSD(T)/adz on water

    """
    resi = {"molecule": h2o, "driver": "energy", "model": {"method": "ccsd(t)", "basis": basis}, "keywords": keywords}

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    ccsd_t_tot = -76.276030676767

    atol = 1.0e-6
    assert compare_values(ccsd_t_tot, res["return_result"], atol=atol)


@pytest.mark.parametrize(
    "program,basis,keywords,errmsg",
    [
        pytest.param(
            "nwchem",
            "aug-cc-pvdz",
            {"ccsd__freeze": 1, "scf__uhf": True},
            "ccsd: nopen is not zero",
            marks=using("nwchem"),
        ),
        pytest.param(
            "gamess",
            "accd",
            {"contrl__scftyp": "uhf"},
            "CCTYP IS PROGRAMMED ONLY FOR SCFTYP=RHF OR ROHF",
            marks=using("gamess"),
        ),
    ],
)
def test_sp_ccsd_t_uhf_fc_error(program, basis, keywords, nh2, errmsg):
    resi = {"molecule": nh2, "driver": "energy", "model": {"method": "ccsd(t)", "basis": basis}, "keywords": keywords}

    with pytest.raises(qcng.exceptions.InputError) as e:
        qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert errmsg in str(e.value)


@pytest.mark.parametrize(
    "program,basis,keywords",
    [
        pytest.param(
            "cfour",
            "aug-pvdz",
            {"reference": "rohf", "occupation": [[3, 1, 1, 0], [3, 0, 1, 0]], "scf_conv": 12, "cc_conv": 12},
            marks=using("cfour"),
        ),
        pytest.param("cfour", "aug-pvdz", {"reference": "rohf"}, marks=using("cfour")),
        # pytest.param('nwchem', 'aug-cc-pvdz', {'basis__spherical': True, 'qc_module': 'tce', 'scf__rohf': True}, marks=using("nwchem")),
        pytest.param("psi4", "aug-cc-pvdz", {"reference": "rohf"}, marks=using("psi4")),
    ],
)
def test_sp_ccsd_t_rohf_full(program, basis, keywords, nh2):
    resi = {"molecule": nh2, "driver": "energy", "model": {"method": "ccsd(t)", "basis": basis}, "keywords": keywords}

    res = qcng.compute(resi, program, raise_error=True, return_dict=True)

    assert res["driver"] == "energy"
    assert "provenance" in res
    assert res["success"] is True

    # aug-cc-pvdz
    ccsd_t_tot = -55.752861467462

    atol = 1.0e-6
    assert compare_values(ccsd_t_tot, res["return_result"], atol=atol)
