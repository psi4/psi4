import pytest

from utils import compare_values

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]

@pytest.fixture
def mols():
    return {
        "h2o" : psi4.geometry("""
        0 1
        O
        H 1 0.96
        H 1 0.96 2 104.5
        symmetry c1
        no_reorient
        no_com
        """),
        "benzene" : psi4.geometry("""
        0 1
        C    -1.0478252   -1.4216736    0.0000000
        C    -1.4545034   -0.8554459    1.2062048
        C    -1.4545034   -0.8554459   -1.2062048
        C    -2.2667970    0.2771610    1.2069539
        C    -2.6714781    0.8450211    0.0000000
        C    -2.2667970    0.2771610   -1.2069539
        H    -1.1338534   -1.2920593   -2.1423150
        H    -2.5824943    0.7163066   -2.1437977
        H    -3.3030422    1.7232700    0.0000000
        H    -2.5824943    0.7163066    2.1437977
        H    -1.1338534   -1.2920593    2.1423150
        H    -0.4060253   -2.2919049    0.0000000
        symmetry c1
        no_reorient
        no_com
        """),
        "benzene+" : psi4.geometry("""
        1 2
        C    -1.0478252   -1.4216736    0.0000000
        C    -1.4545034   -0.8554459    1.2062048
        C    -1.4545034   -0.8554459   -1.2062048
        C    -2.2667970    0.2771610    1.2069539
        C    -2.6714781    0.8450211    0.0000000
        C    -2.2667970    0.2771610   -1.2069539
        H    -1.1338534   -1.2920593   -2.1423150
        H    -2.5824943    0.7163066   -2.1437977
        H    -3.3030422    1.7232700    0.0000000
        H    -2.5824943    0.7163066    2.1437977
        H    -1.1338534   -1.2920593    2.1423150
        H    -0.4060253   -2.2919049    0.0000000
        symmetry c1
        no_reorient
        no_com
        """)
    }

@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : "h2o",
                      "bsse_type" : None,
                      "ref" : -76.0104095301084},
                      id="h2o (rhf)"),
        pytest.param({"method" : "b3lyp",
                      "options": {"reference" : "rhf"},
                      "molecule" : "benzene",
                      "bsse_type" : None,
                      "ref" : -232.24852940340557},
                      id="benzene (rks)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "uhf"},
                      "molecule" : "benzene+",
                      "bsse_type" : None,
                      "ref" : -230.41729673026288},
                      id="benzene+ (uhf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rohf"},
                      "molecule" : "benzene+",
                      "bsse_type" : None,
                      "ref" : -230.40287183406298},
                      id="benzene+ (rohf)"),
    ],
)

def test_dfjcosk(inp, mols):
    """Test the DFJ+LinK Composite JK algorithm via SCF calculations"""

    molecule = mols[inp["molecule"]]
    psi4.set_options({"scf_type" : "linK", "basis": "6-31G*"})
    psi4.set_options(inp["options"])

    # does the DFJLinK SCF energy match a pre-computed reference?
    energy_dfjlink = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"])
    assert compare_values(inp["ref"], energy_dfjlink, atol=1e-6)

    # is the DFJLinK SCF energy reasonably close to a conventional SCF?
    psi4.set_options({"scf_type" : "pk"})
    energy_pk = psi4.energy(inp["method"], molecule=molecule, bsse_type=inp["bsse_type"])
    assert compare_values(energy_pk, energy_dfjlink, atol=1e-4)