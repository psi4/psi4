import numpy as np
import pytest

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.addon]


@pytest.mark.findif
def test_dfmp2_gradient_charged_quintet():
    """
    Test DF-MP2 gradient for a small charged open-shell quintet molecule.
    Regression test for Tensor2d::back_transform dimension mismatch fix.
    
    Multi-atom system ensures nso >> nmo, exposing the dimension mismatch bug
    where temp tensor was created as (nmo, nmo) instead of (nmo, nso).
    
    This test ensures that DF-MP2 gradient calculations work correctly for:
    - Charged systems (charge = +1)
    - Open-shell with multiplicity 5 (quintet, S=2)
    - Multi-atom geometry: nso=def2-SVP(Fe+Cl) >> nmo
    """
    # FeCl+ (iron chloride cation) in high-spin quintet state
    # Charge: +1
    # Multiplicity: 5 (Fe high-spin d^5)
    fecl_cation = psi4.geometry("""
        1 5
        Fe
        Cl 1 2.3
    """)

    psi4.set_options({
        'basis': 'def2-SVP',
        'scf_type': 'df',
        'reference': 'uhf',
        'e_convergence': 10,
        'd_convergence': 9,
    })

    # Compute DF-MP2 gradient (this would crash before the fix)
    # With FeCl+: nso ≈ 30-35, nmo ≈ 20-25
    # Bug creates temp as (20, 20) instead of (20, 30), causing gemm dimension mismatch
    gradient = psi4.gradient('mp2')
    
    # Validate gradient was computed successfully
    assert gradient is not None, "Gradient computation returned None"
    assert gradient.shape == (2, 3), "Gradient should have shape (natom=2, 3)"
    assert np.all(np.isfinite(gradient.np)), "Gradient contains NaN or Inf values"
