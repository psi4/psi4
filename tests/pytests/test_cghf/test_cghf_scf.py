import pytest
import numpy as np

import psi4
from addons import using

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.smoke, pytest.mark.cghf, *using("einsums")]

REFERENCE_ENERGY = -1.1287000935604175

def test_cghf_basic_scf():
    """Minimal CGHF full SCF example: CGHF with no spin-orbit coupling should
    reproduce RHF for closed-shell molecules."""
    mol = psi4.geometry("""
    0 1
        H
        H 1 0.74
    symmetry c1
    """)
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "diis": False,
        "scf_initial_accelerator": "none"
    })
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", psi4.core.get_global_option("BASIS"))
    e = psi4.energy("scf", molecule=mol)
    assert type(e) == float
    assert e == pytest.approx(REFERENCE_ENERGY, abs=1e-6)


def test_cghf_sad_guess():
    """CGHF with GUESS=SAD should converge to the same closed-shell energy as RHF."""
    mol = psi4.geometry("""
    0 1
        H
        H 1 0.74
    symmetry c1
    """)
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "sad",
        "scf_type": "direct",
        "df_scf_guess": False,
        "diis": False,
        "scf_initial_accelerator": "none",
    })
    e = psi4.energy("scf", molecule=mol)
    assert type(e) == float
    assert e == pytest.approx(REFERENCE_ENERGY, abs=1e-6)


def test_cghf_read_guess():
    """Same-basis CGHF READ: checkpoint from a converged SCF restarts in one iteration."""
    mol = psi4.geometry("""
    0 1
        H
        H 1 0.74
    symmetry c1
    """)
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "diis": False,
        "scf_initial_accelerator": "none",
    })
    e0, wfn = psi4.energy("scf", molecule=mol, return_wfn=True, write_orbitals="pytest_cghf_mos")
    assert e0 == pytest.approx(REFERENCE_ENERGY, abs=1e-6)

    psi4.core.clean()
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "scf_type": "direct",
        "df_scf_guess": False,
        "diis": False,
        "scf_initial_accelerator": "none",
        "maxiter": 1,
    })
    e1 = psi4.energy("scf", molecule=mol, restart_file="pytest_cghf_mos")
    assert e1 == pytest.approx(REFERENCE_ENERGY, abs=1e-6)


def test_cghf_guess_C_kwarg():
    """CGHF direct guess_C kwarg: a plain NumPy array of occupied spinor MOs,
    handed straight to energy(), reproduces a converged calculation in one
    iteration -- no checkpoint file or manually-built Wavefunction needed."""
    mol = psi4.geometry("""
    0 1
        H
        H 1 0.74
    symmetry c1
    """)
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "diis": False,
        "scf_initial_accelerator": "none",
    })
    e0, wfn0 = psi4.energy("scf", molecule=mol, return_wfn=True)
    assert e0 == pytest.approx(REFERENCE_ENERGY, abs=1e-6)

    nocc = wfn0.nelec()
    C_occ = wfn0.C().to_array()[:, :nocc]  # plain NumPy array, not a ComplexMatrix

    psi4.core.clean()
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "scf_type": "direct",
        "df_scf_guess": False,
        "diis": False,
        "scf_initial_accelerator": "none",
        "maxiter": 1,
    })
    e1 = psi4.energy("scf", molecule=mol, guess_C=C_occ)
    assert e1 == pytest.approx(REFERENCE_ENERGY, abs=1e-6)


def test_cghf_read_guess_basis_projection_nyi():
    """Cross-basis CGHF READ should raise until spinor basis_projection exists."""
    mol = psi4.geometry("""
    0 1
        H
        H 1 0.74
    symmetry c1
    """)
    psi4.set_options({
        "basis": "sto-3g",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "diis": False,
        "scf_initial_accelerator": "none",
    })
    psi4.energy("scf", molecule=mol, write_orbitals="pytest_cghf_sto3g")

    psi4.core.clean()
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "scf_type": "direct",
        "df_scf_guess": False,
    })
    with pytest.raises(NotImplementedError, match="basis projection"):
        psi4.energy("scf", molecule=mol, restart_file="pytest_cghf_sto3g")
