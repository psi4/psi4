import pytest

import psi4

from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]

@pytest.mark.parametrize("inp", [
    pytest.param({'options': {"guess": "core"}}, id="core"),
    pytest.param({'options': {"guess": "gwh"}}, id="gwh"),
    pytest.param({'options': {"guess": "huckel"}}, id="huckel"),
    pytest.param({'options': {"guess": "modhuckel"}}, id="modhuckel"),
    pytest.param({'late_options': {"guess": "read"}}, id="read"),
    pytest.param({'options': {"guess": "sad"}}, id="sad"),
    pytest.param({'options': {"guess": "sadno"}}, id="sadno"),
    pytest.param({'options': {"guess": "sap"}}, id="sap"),
    ]
)
def test_guess_mix_for_broken_symmetry(inp):

    refENuc  =  0.17639240356
    refSCF   = -0.82648407827446 
    refBSSCF = -0.99872135103903

    h2 = psi4.geometry("""
        0 1
        H
        H 1 3.0
        symmetry c1
        """)

    psi4.set_options({"reference": "uhf", "e_convergence": 12, "basis": "cc-pvdz"})
    psi4.set_options(inp.get("options", {}))

    thisSCF = psi4.energy("scf")
    psi4.set_options(inp.get("late_options", {}))
    psi4.set_options({"guess_mix": True})
    thisBSSCF = psi4.energy("scf")

    assert compare_values(refENuc, h2.nuclear_repulsion_energy(), 10, "Nuclear repulsion energy")
    assert compare_values(refSCF, thisSCF, 10, "Reference energy")
    assert compare_values(refBSSCF, thisBSSCF, 10, "Reference broken-symmetry energy")

