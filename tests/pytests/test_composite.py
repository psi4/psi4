import pytest
import psi4
from utils import *
from addons import uusing

pytestmark = [pytest.mark.psi, pytest.mark.api]


@uusing("mrcc")
@pytest.mark.cbs
def test_allen_focal_point():
    be = psi4.geometry("Be")
    
    psi4.energy("allen_focal_point", scf_basis="cc-pV[TQ5]Z", corl_basis="cc-pV[Q5]Z", delta_basis="cc-pV[Q5]Z", delta2_basis="cc-pV[Q5]Z")
    psi4.core.print_variables()
    
    assert compare_values(-14.57305004, psi4.variable("CBS REFERENCE ENERGY"), 7, "scf tq5")
    #ADD VAR assert compare_values(-0.06148737, psi4.variable("CBS CORL ENERGY"), 7, "corl mp2 q5")
    assert compare_values(-0.01767880, psi4.variable("CBS DELTA1 TOTAL ENERGY"), 7, "delta1 ccsd q5")
    assert compare_values(-0.00063008, psi4.variable("CBS DELTA2 TOTAL ENERGY"), 7, "delta2 (t) q5")
    assert compare_values(-0.00001090, psi4.variable("CBS DELTA3 TOTAL ENERGY"), 7, "delta3 t 3")
    assert compare_values(-0.00000205, psi4.variable("CBS DELTA4 TOTAL ENERGY"), 7, "delta4 (q) 2")
    assert compare_values(-14.65285924, psi4.variable("CBS TOTAL ENERGY"), 7, "cbs")
    assert compare_values(-14.65285924, psi4.variable("CURRENT ENERGY"), 7, "current")
    assert compare(5, psi4.variable("CBS NUMBER"), "cbs no")

    # 5 tasks: (t)/tz, (t)/qz, (t)/5z, t/tz, t(q)/dz
    # lowering the basis sets keeps the test within packaged L2's AM bounds

@pytest.mark.cbs
def test_cbs_with_managed_conditions():

    n = psi4.geometry("""
            0 4
            N 0.00 0.00 0.00
        """)

    psi4.set_options({
        "scf_type": "direct",
        "reference": "rohf",
        "r_convergence": 6,
        "d_convergence": 7,
        "e_convergence": 8,
        "freeze_core": True,
    })

    e_cbs = psi4.energy('cbs',
        scf_basis='aug-cc-pV[TQ5]Z',
        corl_wfn='mp2',
        corl_basis='aug-cc-pV[TQ]Z',
        delta_wfn='ccsd(t)',
        delta_basis='aug-cc-pV[DT]Z')

    assert compare_values(-54.5325027481897, e_cbs, 5, "cbs")

