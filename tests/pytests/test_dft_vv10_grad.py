import numpy as np
import pytest

import psi4

from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.scf, pytest.mark.dft]


@pytest.mark.findif
@pytest.mark.slow
@pytest.mark.parametrize("func,ref_key", [
    pytest.param("BLYP-NL", "rks", id="BLYP-NL RKS"),
    pytest.param("BLYP-NL", "uks", id="BLYP-NL UKS"),
])
def test_vv10_gradient(func, ref_key):
    """Compare analytic vs finite-difference gradients for VV10 functionals."""

    psi4.set_num_threads(1)

    mol = psi4.geometry("""
        0 1
        N   -0.0034118    3.5353926    0.0000000
        C    0.0751963    2.3707040    0.0000000
        H    0.1476295    1.3052847    0.0000000
    """)

    psi4.set_options({
        "basis": "cc-pvdz",
        "reference": ref_key,
        "scf_type": "df",
        "dft_radial_points": 99,
        "dft_spherical_points": 590,
        "dft_vv10_spherical_points": 194,
        "dft_vv10_radial_points": 50,
        "e_convergence": 9,
        "d_convergence": 9,
        "points": 5,
    })

    analytic_gradient = psi4.gradient(func, dertype=1)
    psi4.core.clean()
    findif_gradient = psi4.gradient(func, dertype=0)

    assert compare_values(findif_gradient, analytic_gradient, 4,
                          f"VV10 {func} {ref_key.upper()} analytic vs. findif gradient")
    psi4.core.clean()
