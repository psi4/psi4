import numpy as np
import pytest

import psi4

from utils import compare_values

pytestmark = [pytest.mark.psi, pytest.mark.api]

# Reference data generated from Psi's dfmp2 module
data = {
    "df-mp2 ae": np.array([[ 0, 0,  9.62190509e-03],
         [0,  5.49835030e-03, -4.81095255e-03],
         [0, -5.49835030e-03, -4.81095255e-03]]),
    "df-mp2 fc": np.array([[ 0, 0,  1.02432654e-02],
         [0,  5.88581965e-03, -5.12163268e-03],
         [0, -5.88581965e-03, -5.12163268e-03]]),
    "df-mp2 fv": np.array([[ 0, 0,  1.22918833e-02],
         [0,  5.52107556e-03, -6.14594166e-03],
         [0, -5.52107556e-03, -6.14594166e-03]]),
    "df-mp2 fc/fv": np.array([[ 0, 0,  1.25984187e-02],
         [0,  5.71563223e-03, -6.29920936e-03],
         [0, -5.71563223e-03, -6.29920936e-03]]),
    "df-dct": np.array([[0, 0, 0.008477558394],
        [0,  0.005148825942, -0.004238779197],
        [0, -0.005148825942, -0.004238779197]]),
    "df-cc2": np.array([[0, 0, 0.011903811700],
        [0,  0.006730035450, -0.005951905850],
        [0, -0.006730035450, -0.005951905850]])
    }

@pytest.mark.findif
@pytest.mark.slow
# TODO: That "true" needs to be a string is silly. Convert it to a boolean when you can do that without incurring a NaN energy.
@pytest.mark.parametrize("inp", [
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df'}, 'ref': data["df-mp2 ae"]}, id='df-mp2 ae'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df', 'freeze_core': 'true'}, 'ref': data["df-mp2 fc"]}, id='df-mp2 fc'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df', 'num_frozen_uocc': 4}, 'ref': data["df-mp2 fv"]}, id='df-mp2 fv'),
    pytest.param({'name': 'mp2', 'options': {'mp2_type': 'df', 'freeze_core': 'true', 'num_frozen_uocc': 4}, 'ref': data["df-mp2 fc/fv"]}, id='df-omp2 fc/fv'),
    pytest.param({'name': 'dct', 'options': {'dct_type': 'df'}, 'ref': data["df-dct"]}, id='df-rdct'),
    pytest.param({'name': 'dct', 'options': {'dct_type': 'df', 'reference': 'uhf'}, 'ref': data["df-dct"]}, id='df-udct'),
    pytest.param({'name': 'cc2', 'options': {'cc_type': 'df', 'qc_module': 'ccenergy'}, 'ref': data["df-cc2"]}, id='df-cc2')
    ]
)
def test_gradient(inp):
    h2o = psi4.geometry("""
        O
        H 1 0.958
        H 1 0.958 2 104.5
    """)

    psi4.set_options({'basis': 'aug-cc-pvdz', 'points': 5})
    psi4.set_options(inp['options'])
    
    analytic_gradient = psi4.gradient(inp['name'], dertype=1)
    print(analytic_gradient)
    findif_gradient = psi4.gradient(inp['name'], dertype=0)
    reference_gradient = inp["ref"]

    assert compare_values(findif_gradient, analytic_gradient, 5, "analytic vs. findif gradient")
    assert compare_values(reference_gradient, analytic_gradient.np, 5, "analytic vs. reference gradient")

def test_gradient_ref():
    h2o = psi4.geometry("""
        O
        H 1 0.958
        H 1 0.958 2 104.5
    """)

    psi4.set_options({"basis": "cc-pVDZ"})
    mp2_wfn = psi4.energy("mp2", return_wfn=True)[1]
    with pytest.raises(TypeError):
        psi4.gradient("scf", ref_wfn=mp2_wfn)
    scf_wfn = psi4.energy("scf", return_wfn=True)[1]
    psi4.gradient("scf", ref_wfn=scf_wfn)


# Distorted, non-planar, no symmetry, so all twelve gradient components are independent and nonzero
_vv10_molecules = {
    "nh3": """0 1
        N  0.000  0.000  0.120
        H  0.000  0.980 -0.270
        H  0.870 -0.450 -0.330
        H -0.820 -0.520 -0.250
        symmetry c1
        no_reorient
        no_com""",
    "ch3": """0 2
        C  0.000  0.000  0.060
        H  0.000  1.120 -0.090
        H  0.950 -0.540 -0.160
        H -0.930 -0.560 -0.110
        symmetry c1
        no_reorient
        no_com""",
}


@pytest.mark.findif
@pytest.mark.slow
@pytest.mark.parametrize("mol,reference,functional,grid,atol", [
    pytest.param("nh3", "rks", "wb97x-v", {}, 3.e-5, id="nh3-rks-wb97x-v"),
    pytest.param("ch3", "uks", "wb97x-v", {}, 3.e-5, id="ch3-uks-wb97x-v"),
    # B97M-V's meta-GGA is grid sensitive: on the default local grid analytic and findif differ by 1.3e-3 here
    pytest.param("nh3", "rks", "b97m-v", {"dft_radial_points": 99, "dft_spherical_points": 590}, 2.e-4,
                 id="nh3-rks-b97m-v", marks=pytest.mark.medlong),
])
def test_vv10_gradient(mol, reference, functional, grid, atol):
    """Analytic gradients with VV10 nonlocal correlation against 3-point finite differences of energies. As for the
    local XC terms, the analytic gradient holds the grid fixed, so the tolerances are its error on these
    grids (observed: 9.6e-6, 3.5e-6, 6.6e-5), which plain PBE shares (6.5e-5 for this NH3)."""
    psi4.geometry(_vv10_molecules[mol])
    psi4.set_options({"basis": "def2-svp", "scf_type": "pk", "reference": reference, "e_convergence": 11,
                      "d_convergence": 9, "maxiter": 300, **grid})

    analytic = psi4.gradient(functional, dertype=1).np
    findif = psi4.gradient(functional, dertype=0).np

    assert compare_values(findif, analytic, atol=atol, label=f"{functional} analytic vs. findif gradient")


def test_vv10_postscf_gradient_by_findif():
    """A VV10 energy added after the SCF (DFT_VV10_POSTSCF) is not stationary in the orbitals, so its gradient
    is taken by finite differences, whereas the self-consistent one is analytic."""
    from psi4.driver.driver_util import negotiate_derivative_type

    psi4.geometry(_vv10_molecules["nh3"])
    assert negotiate_derivative_type("gradient", "wb97x-v", None) == (1, 1)
    psi4.set_options({"dft_vv10_postscf": True})
    assert negotiate_derivative_type("gradient", "wb97x-v", None) == (1, 0)
    assert negotiate_derivative_type("gradient", "pbe", None) == (1, 1)
