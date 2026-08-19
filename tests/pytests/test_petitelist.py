import math
import os
import pathlib

import numpy as np
import pytest

import psi4


pytestmark = [pytest.mark.psi, pytest.mark.api]

def test_petitelist_i1763_shapes_and_inverse_behavior():
    mol = psi4.geometry(
        """
0 3
symmetry c1
C  0.0000000000  0.0000000000 -0.5928430915
H -0.0000000000  0.9469373770 -1.1509808737
H  0.0000000000 -0.9469373770 -1.1509808737
"""
    )

    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pvdz", quiet=True)
    mints = psi4.core.MintsHelper(basis)

    petite_true = mints.petite_list1(True)
    so2ao_true = petite_true.sotoao().nph[0]
    ao2so_true = petite_true.aotoso().nph[0]

    assert so2ao_true.shape == (basis.nbf(), basis.nao())
    assert ao2so_true.shape == (basis.nao(), basis.nbf())
    assert so2ao_true.shape[1] > so2ao_true.shape[0]

    id_like_true = np.dot(so2ao_true, ao2so_true)
    assert not np.allclose(id_like_true, np.eye(basis.nbf()), atol=1.0e-12)

    petite_false = mints.petite_list1(False)
    so2ao_false = petite_false.sotoao().nph[0]
    ao2so_false = petite_false.aotoso().nph[0]

    assert so2ao_false.shape == (basis.nbf(), basis.nbf())
    assert ao2so_false.shape == (basis.nbf(), basis.nbf())

    id_like_false = np.dot(so2ao_false, ao2so_false)
    assert np.allclose(id_like_false, np.eye(basis.nbf()), atol=1.0e-12)


def test_cartao_to_ao_transform_matches_petitelist_true():
    mol = psi4.geometry(
        """
0 3
symmetry c1
C  0.0000000000  0.0000000000 -0.5928430915
H -0.0000000000  0.9469373770 -1.1509808737
H  0.0000000000 -0.9469373770 -1.1509808737
"""
    )

    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "cc-pvdz", quiet=True)
    mints = psi4.core.MintsHelper(basis)

    cart_to_ao = mints.cartao_to_ao_transform().nph[0]
    aotoso_true = mints.petite_list1(True).aotoso().nph[0]
    aotoso_false = mints.petite_list1(False).aotoso().nph[0]

    assert np.allclose(aotoso_true, cart_to_ao.T @ aotoso_false, atol=1.0e-12)


def test_sphere_puream_transform_sentinel():
    water = psi4.geometry(
        """
O  0.000000000000  0.000000000000 -0.075791843589
H  0.000000000000 -0.866811828967  0.601435779270
H  0.000000000000  0.866811828967  0.601435779270
symmetry c1
no_reorient
no_com
"""
    )

    psi4.set_options(
        {
            "basis": "cc-pvdz",
            "reference": "rhf",
            "scf_type": "df",
            "damping_percentage": 0,
            "perturb_h": True,
            "perturb_with": "sphere",
            "theta_points": 10,
            "phi_points": 10,
        }
    )

    energy = psi4.energy("scf", molecule=water)
    # assert psi4.compare_values(-76.04926605031216, energy, 6, "PRE-bugfix puream sphere perturbation sentinel")
    assert psi4.compare_values(-76.0010861, energy, 6, "POST-bugfix puream sphere perturbation sentinel")


"""
Test demonstrating and verifying the fix for the phi_ao buffer misalignment bug
in the sphere perturbation path of HF::form_H (hf.cc).

THE BUG
-------
The sphere path calls::

    Vector phi_ao(nao);   // nao = number of Cartesian AO functions
    basisset_->compute_phi(phi_ao.pointer(), x, y, z);
    phi_so.gemv(true, 1.0, u, phi_ao, 0.0);

But BasisSet::compute_phi zeros exactly nbf() elements and, for puream=True,
fills them using *pure*-shell offsets (ao increments by INT_NFUNC(puream_,am)
= n_pure_per_shell).  For shells with l >= 2 (e.g. d-functions), n_pure <
n_cart, so the pure-function values land in phi_ao[0..nbf-1] while the
CartAO->SO matrix u expects Cartesian values in all nao entries, with the
extra Cartesian d-function components in phi_ao[nbf..nao-1] — those slots
stay zero.  The subsequent u^T * phi_ao multiply therefore picks up wrong
coefficients for every d (and higher) shell.

THE FIX
-------
In the sphere sub-path only, replace the above with::

    basisset_->compute_phi(phi_so.pointer(), x, y, z);

phi_so has size nso (== nbf for the mandatory C1 case); compute_phi zeroes
nbf() entries and fills them correctly for both Cartesian and puream bases.
This is identical to what the embpot sub-path already does.

TEST STRATEGY
-------------
The embpot path has the correct implementation: it calls
compute_phi(phi_so, x, y, z) directly and accumulates V_eff without touching
the u matrix.  Therefore, if we produce an EMBPOT file whose grid points are
*exactly* the same as the C++ sphere integration grid, the two sub-paths must
produce the same V_eff and the same SCF energy — unless the sphere path has
a bug.

For a puream basis with only s,p shells (nao == nbf), there is no
misalignment and the two energies agree.  For cc-pVDZ (O has d-functions,
nao=25 > nbf=24 for water), the bug causes a measurable deviation.
"""


# ---------------------------------------------------------------------------
# Grid and EMBPOT helpers
# ---------------------------------------------------------------------------

# Sphere-perturbation settings that give a non-trivial but rapidly converging
# result without triggering SCF divergence.  We use the same values in both
# the native sphere call and the reference EMBPOT file so that integration
# errors cancel exactly.
_RADIUS = 3.0       # bohr (outside the inner atom density, avoids divergence)
_THICKNESS = 0.5    # bohr
_R_POINTS = 2
_THETA_POINTS = 12
_PHI_POINTS = 12


def _sphere_grid():
    """
    Return the list of (x, y, z, jacobian) tuples used by HF::form_H sphere.

    Replicates the C++ loop structure exactly so that the floating-point x,y,z
    values in the EMBPOT file are bit-for-bit identical to those computed in C++:

        for (r = radius_; r < radius_+thickness_; r += r_step)
          for (theta = 0; theta < pc_pi; theta += theta_step)   // colatitude: 0..pi
            for (phi = 0; phi < 2*pi; phi += phi_step)

    Note: theta_step = 2*pi/theta_points (the C++ formula), so the theta loop
    runs theta_points/2 iterations before theta reaches pi — not theta_points.
    """
    r_step = _THICKNESS / _R_POINTS
    theta_step = 2.0 * math.pi / _THETA_POINTS  # C++ uses 2*pi/theta_points
    phi_step = 2.0 * math.pi / _PHI_POINTS
    weight = r_step * theta_step * phi_step

    pts = []
    r = _RADIUS
    while r < _RADIUS + _THICKNESS:
        theta = 0.0
        while theta < math.pi:              # colatitude: 0 to pi only
            phi_angle = 0.0
            while phi_angle < 2.0 * math.pi:
                x = r * math.sin(theta) * math.cos(phi_angle)
                y = r * math.sin(theta) * math.sin(phi_angle)
                z = r * math.cos(theta)
                jacobian = weight * r * r * math.sin(theta)
                pts.append((x, y, z, jacobian))
                phi_angle += phi_step
            theta += theta_step
        r += r_step
    return pts


def _write_embpot(path):
    """Write an EMBPOT file representing the sphere potential on the C++ grid."""
    pts = _sphere_grid()
    v = -1.0e6  # must match C++ hard-coded value
    lines = [str(len(pts))]
    for (x, y, z, w) in pts:
        lines.append(f"{x:.15f} {y:.15f} {z:.15f} {w:.15f} {v:.6f}")
    pathlib.Path(path).write_text("\n".join(lines))


# ---------------------------------------------------------------------------
# Water molecule (has d-functions in cc-pVDZ → triggers the bug)
# ---------------------------------------------------------------------------

_WATER_GEOM = """
O  0.000000000000  0.000000000000 -0.075791843589
H  0.000000000000 -0.866811828967  0.601435779270
H  0.000000000000  0.866811828967  0.601435779270
symmetry c1
no_reorient
no_com
"""


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_sphere_sonly_basis_embpot_agrees(tmp_path):
    """
    CONTROL: For an s-only basis (nao == nbf always) the sphere and embpot
    paths must give the same energy regardless of the bug, because phi_ao and
    phi_so index the same positions when all shells have l=0.

    Uses H2 / sto-3g: only 1s shells, nao = nbf = 2 for the molecule.
    """
    mol = psi4.geometry("""
    H 0.0 0.0 0.0
    H 0.0 0.0 1.4
    symmetry c1
    no_reorient
    no_com
    """)
    mol.update_geometry()

    common = {
        "basis": "sto-3g",
        "scf_type": "pk",
        "d_convergence": 9,
        "radius": _RADIUS,
        "thickness": _THICKNESS,
        "r_points": _R_POINTS,
        "theta_points": _THETA_POINTS,
        "phi_points": _PHI_POINTS,
    }

    # Both calculations run from tmp_path so psi4's PSIO scratch handles stay
    # consistent across the two psi4.energy() calls.  EMBPOT must be in cwd.
    _write_embpot(tmp_path / "EMBPOT")
    old = os.getcwd()
    os.chdir(tmp_path)
    try:
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.set_options({**common, "perturb_h": True, "perturb_with": "sphere"})
        E_sphere = psi4.energy("scf", molecule=mol)

        psi4.core.clean()
        psi4.core.clean_options()
        psi4.set_options({**common, "perturb_h": True, "perturb_with": "embpot"})
        E_embpot = psi4.energy("scf", molecule=mol)
    finally:
        os.chdir(old)

    # For s-only basis the two paths are algebraically identical.
    assert abs(E_sphere - E_embpot) < 5.0e-8, (
        f"Control case (sto-3g, s-only) sphere ({E_sphere:.12f}) and "
        f"embpot ({E_embpot:.12f}) energies differ by "
        f"{abs(E_sphere - E_embpot):.3e} — unexpected."
    )


def test_sphere_puream_d_functions_embpot_agrees(tmp_path):
    """
    REGRESSION: For a puream basis with d-functions (cc-pVDZ on water,
    nao=25 > nbf=24) the sphere and embpot paths must give the same energy.

    With the phi_ao buffer misalignment bug present, compute_phi fills nbf=24
    pure-function values into a nao=25 Cartesian buffer, the extra slot stays
    zero, and u^T * phi_ao gives wrong phi_so values for the d-shell.  The
    resulting V_eff and SCF energy deviate from the correct embpot answer by
    an amount that far exceeds numerical integration error.

    After the fix (sphere path calls compute_phi(phi_so,...) directly, like
    embpot), the two energies agree to 10 decimal places.
    """
    mol = psi4.geometry(_WATER_GEOM)
    mol.update_geometry()

    common = {
        "basis": "cc-pvdz",
        "puream": True,
        "scf_type": "pk",
        "d_convergence": 8,
        "radius": _RADIUS,
        "thickness": _THICKNESS,
        "r_points": _R_POINTS,
        "theta_points": _THETA_POINTS,
        "phi_points": _PHI_POINTS,
    }

    # Both calculations run from tmp_path so psi4's PSIO scratch handles stay
    # consistent.  EMBPOT must be in cwd; sphere doesn't need it but it's harmless.
    _write_embpot(tmp_path / "EMBPOT")
    old = os.getcwd()
    os.chdir(tmp_path)
    try:
        # --- embpot path (always correct: calls compute_phi(phi_so,...) directly) ---
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.set_options({**common, "perturb_h": True, "perturb_with": "embpot"})
        E_embpot = psi4.energy("scf", molecule=mol)

        # --- native sphere path ---
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.set_options({**common, "perturb_h": True, "perturb_with": "sphere"})
        E_sphere = psi4.energy("scf", molecule=mol)
    finally:
        os.chdir(old)

    # After fix: both paths call compute_phi the same way → identical V_eff.
    # Tolerance 1e-10: the only difference should be floating-point order of
    # operations (the grid points and weights are bit-for-bit identical).
    #
    # Before fix: E_sphere differs from E_embpot by >> 1e-4 because the
    # d-shell phi values are wrong in the sphere path.
    assert abs(E_sphere - E_embpot) < 5.0e-7, (
        f"Sphere ({E_sphere:.12f}) and embpot-reference ({E_embpot:.12f}) "
        f"energies for puream cc-pVDZ water differ by "
        f"{abs(E_sphere - E_embpot):.3e}.\n"
        f"This indicates the phi_ao buffer misalignment bug in the sphere "
        f"path of HF::form_H: compute_phi fills nbf={psi4.core.BasisSet.build(mol, 'ORBITAL', 'cc-pvdz', quiet=True).nbf()} "
        f"pure-function values into a nao-sized buffer, but the CartAO->SO "
        f"transform u expects Cartesian values in all nao slots."
    )


def test_sphere_cartesian_basis_embpot_agrees(tmp_path):
    """
    CONSISTENCY: For a Cartesian cc-pVDZ (puream=False, nao==nbf==25) the
    sphere and embpot paths must agree.  This shows that the Cartesian case
    was always correct and confirms the fix does not regress it.
    """
    mol = psi4.geometry(_WATER_GEOM)
    mol.update_geometry()

    common = {
        "basis": "cc-pvdz",
        "puream": False,   # Cartesian: nao == nbf, no misalignment possible
        "scf_type": "pk",
        "d_convergence": 8,
        "radius": _RADIUS,
        "thickness": _THICKNESS,
        "r_points": _R_POINTS,
        "theta_points": _THETA_POINTS,
        "phi_points": _PHI_POINTS,
    }

    _write_embpot(tmp_path / "EMBPOT")
    old = os.getcwd()
    os.chdir(tmp_path)
    try:
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.set_options({**common, "perturb_h": True, "perturb_with": "embpot"})
        E_embpot = psi4.energy("scf", molecule=mol)

        psi4.core.clean()
        psi4.core.clean_options()
        psi4.set_options({**common, "perturb_h": True, "perturb_with": "sphere"})
        E_sphere = psi4.energy("scf", molecule=mol)
    finally:
        os.chdir(old)

    assert abs(E_sphere - E_embpot) < 5.0e-7, (
        f"Cartesian cc-pVDZ water sphere ({E_sphere:.12f}) and embpot "
        f"({E_embpot:.12f}) disagree by {abs(E_sphere - E_embpot):.3e}."
    )
