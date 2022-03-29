import pytest
from utils import *
from addons import uusing

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]

@uusing("dftd3")
def test_dftd3_dft_grad_lr3():
    """modified VV10-less B97 functional gradient wB97X-V -> wB97X-D3BJ"""

    # stored data from finite differences
    FD_wb97x_d3 = psi4.core.Matrix.from_list([
       [  0.03637802044642,    0.06718963272193,    0.00000000000000],
       [  0.04955519892514,   -0.06340333481039,    0.00000000000000],
       [ -0.07009043821383,   -0.00834477190196,    0.00000000000000],
       [  0.02732425404378,   -0.05883094637658,    0.00000000000000],
       [ -0.02158351760075,    0.03169471018350,    0.05342791683461],
       [ -0.02158351760075,    0.03169471018350,   -0.05342791683461]])

    psi4.geometry("""
    0 1
    O         -1.65542       -0.12330        0.00000
    O          1.24621        0.10269        0.00000
    H         -0.70409        0.03193        0.00000
    H         -2.03867        0.75372        0.00000
    H          1.57599       -0.38252       -0.75856
    H          1.57599       -0.38252        0.75856
    """)

    psi4.set_options({
        'scf_type': 'pk',
        'basis': 'minix',
        'dft_radial_points': 99,
        'dft_spherical_points': 302,
        'e_convergence': 8,
    })

    analytic = psi4.gradient('wB97X-D3BJ', dertype=1)
    assert compare_matrices(analytic, FD_wb97x_d3, 5, "wB97X-D3BJ Analytic vs Store")

