import pytest

from utils import compare_values

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]

@pytest.fixture
def mols():
    return {
        "h2o" : psi4.geometry("""
0 1
O    0.000000000000     0.000000000000    -0.124038860300
H    0.000000000000    -1.431430901356     0.984293362719
H    0.000000000000     1.431430901356     0.984293362719
units au
"""),
        "nh2" : psi4.geometry("""
0 2
N    0.000000000000000   0.000000000000000  -0.145912918634892
H    0.000000000000000  -1.511214298139000   1.013682596946108
H    0.000000000000000   1.511214298139000   1.013682596946108
units au
"""),
        "h2o_nap1" : psi4.geometry("""
0 1
O    0.000000000000     0.000000000000    -0.124038860300
H    0.000000000000    -1.431430901356     0.984293362719
H    0.000000000000     1.431430901356     0.984293362719
--
1 1
Na   0.000000000000     0.000000000000    -4.124038860300
units au
""")
        }


@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      "ref" : -76.026780223322},
                      id="h2o (rhf)"),
        pytest.param({"method" : "b3lyp",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      "ref" : -76.420402720419},
                      id="h2o (rks)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "uhf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      "ref" : -55.566890252551},
                      id="nh2 (uhf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rohf"},
                      "molecule" : "nh2",
                      "bsse_type" : None,
                      "ref" : -55.562689948780},
                      id="nh2 (rohf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o_nap1",
                      "bsse_type" : "CP",
                      "ref" :  -0.040121884077},
                      id="h2o/na+ (rhf ie)"),
    ],
)
def test_dfjcosk(inp, mols):
    """Test the DFJCOSK JK object via SCF calculations"""

    molecule = mols[inp["molecule"]]
    psi4.set_options({"scf_type" : "cosx", "basis": "cc-pvdz"})
    psi4.set_options(inp["options"])

    # does the DFJCOSK SCF energy match a pre-computed reference?
    energy_dfjcosk = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"])
    assert compare_values(inp["ref"], energy_dfjcosk, atol=1e-6)

    # is the DFJCOSK SCF energy reasonably close to a conventional SCF?
    psi4.set_options({"scf_type" : "pk"})
    energy_pk = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"])
    assert compare_values(energy_pk, energy_dfjcosk, atol=1e-4)
