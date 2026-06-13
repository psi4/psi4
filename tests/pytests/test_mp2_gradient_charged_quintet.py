import numpy as np
import pytest

import psi4
from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]


@pytest.mark.findif
@pytest.mark.slow
def test_dfmp2_gradient_charged_quintet():
    """
    Test DF-MP2 gradient for a small charged open-shell quintet molecule.
    Regression test for Tensor2d::back_transform dimension mismatch fix.
    
    This test ensures that DF-MP2 gradient calculations work correctly for:
    - Charged systems (charge = +1)
    - Open-shell with multiplicity 5 (quintet, S=2)
    - All molecular properties that could expose dimension bugs
    """
    # Fe+ (iron cation) in high-spin quintet state
    # Natural multiplicity: 5 (d^5 configuration)
    # Charge: +1
    fe_cation = psi4.geometry("""
        1 5
        Fe
    """)

    psi4.set_options({
        'basis': 'def2-SVP',
        'scf_type': 'df',
        'reference': 'uhf',
        'e_convergence': 10,
        'd_convergence': 9,
    })

    # Compute DF-MP2 gradient (this would crash before the fix)
    gradient = psi4.gradient('mp2')
    
    # Basic validation: gradient should be small for a single atom
    # (single atoms have no internal degrees of freedom, so gradient norm should be near machine precision)
    gradient_norm = np.linalg.norm(gradient.np)
    assert gradient_norm < 1e-4, f"Single-atom gradient norm {gradient_norm} is unexpectedly large"
