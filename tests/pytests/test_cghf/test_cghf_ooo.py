import pytest

import psi4
from addons import using

pytestmark = [
    pytest.mark.psi,
    pytest.mark.api,
    pytest.mark.smoke,
    pytest.mark.cghf,
    *using("einsums"),
    *using("ooo"),
]

# CGHF without SOC should match RHF for closed-shell H2 (same ref as test_cghf_scf.py)
REFERENCE_ENERGY = -1.1287000935604175


def test_cghf_openorbital_scf():
    """CGHF + OpenOrbitalOptimizer should reproduce the internal CGHF/RHF energy."""
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
        "orbital_optimizer_package": "openorbitaloptimizer",
        "e_convergence": 1e-8,
        "d_convergence": 1e-8,
    })
    e = psi4.energy("scf", molecule=mol)
    assert e == pytest.approx(REFERENCE_ENERGY, abs=1e-6)


def test_cghf_openorbital_matches_internal():
    """OOO and INTERNAL CGHF energies should agree for a small closed-shell case."""
    mol = psi4.geometry("""
    0 1
        H
        H 1 0.74
    symmetry c1
    """)
    common = {
        "basis": "sto-3g",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "e_convergence": 1e-10,
        "d_convergence": 1e-10,
    }

    psi4.core.clean()
    psi4.set_options({**common, "orbital_optimizer_package": "internal", "diis": True})
    e_internal = psi4.energy("scf", molecule=mol)

    psi4.core.clean()
    psi4.set_options({**common, "orbital_optimizer_package": "openorbitaloptimizer"})
    e_ooo = psi4.energy("scf", molecule=mol)

    assert e_ooo == pytest.approx(e_internal, abs=1e-8)
