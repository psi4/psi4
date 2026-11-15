import pytest

import psi4

from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.scf, pytest.mark.dft]

_HCN = """
    0 1
    N   -0.0034118    3.5353926    0.0000000
    C    0.0751963    2.3707040    0.0000000
    H    0.1476295    1.3052847    0.0000000
"""

_NH2 = """
    0 2
    N
    H 1 1.03
    H 1 1.03 2 103.0
"""


@pytest.mark.findif
@pytest.mark.slow
@pytest.mark.parametrize("func,ref_key,geom,nl_grid,digits", [
    pytest.param("BLYP-NL", "rks", _HCN, (50, 194), 4, id="BLYP-NL RKS"),
    pytest.param("BLYP-NL", "uks", _NH2, (50, 194), 4, id="BLYP-NL UKS doublet"),
    # Meta-GGA host: regression for stale V_TAU_A in the VV10 gradient path.
    pytest.param("wB97M-V", "rks", _HCN, (50, 194), 4, id="wB97M-V RKS"),
])
def test_vv10_gradient(func, ref_key, geom, nl_grid, digits):
    psi4.set_num_threads(1)

    mol = psi4.geometry(geom)

    psi4.set_options({
        "basis": "cc-pvdz",
        "reference": ref_key,
        "scf_type": "df",
        "dft_radial_points": 99,
        "dft_spherical_points": 590,
        "dft_vv10_radial_points": nl_grid[0],
        "dft_vv10_spherical_points": nl_grid[1],
        "e_convergence": 9,
        "d_convergence": 9,
        "points": 5,
    })

    analytic_gradient = psi4.gradient(func, dertype=1)
    psi4.core.clean()
    findif_gradient = psi4.gradient(func, dertype=0)

    assert compare_values(findif_gradient, analytic_gradient, digits,
                          f"VV10 {func} {ref_key.upper()} analytic vs. findif gradient")
    psi4.core.clean()
