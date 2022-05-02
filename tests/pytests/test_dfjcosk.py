import pytest

from utils import compare_values

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]


h2o = psi4.geometry("""
0 1
O    0.000000000000     0.000000000000    -0.124038860300
H    0.000000000000    -1.431430901356     0.984293362719
H    0.000000000000     1.431430901356     0.984293362719
units au
""")

nh2 = psi4.geometry("""
0 2
N    0.000000000000000   0.000000000000000  -0.145912918634892
H    0.000000000000000  -1.511214298139000   1.013682596946108
H    0.000000000000000   1.511214298139000   1.013682596946108
units au
""")

@pytest.mark.parametrize(
    "inp",
    [
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rhf"},
                      "molecule" : h2o,
                      "ref" : -76.0267902379935},
                      id="h2o (rhf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "uhf"},
                      "molecule" : nh2,
                      "ref" : -55.56691163612347},
                      id="nh2 (uhf)"),
        pytest.param({"method" : "hf",
                      "options": {"reference" : "rohf"},
                      "molecule" : nh2,
                      "ref" : -55.56271068759112},
                      id="nh2 (rohf)"),
    ],

)
def test_dfjcosk(inp):
    """Formerly known as dfmp2-4"""

    #psi4.set_options({"e_convergence" : 1e-6})
    psi4.set_options({"scf_type" : "cosk", "basis": "cc-pvdz"})
    psi4.set_options(inp["options"])

    energy, wfn = psi4.energy(inp["method"], molecule = inp["molecule"], return_wfn=True)
    assert compare_values(inp["ref"], energy, atol=1e-6)

