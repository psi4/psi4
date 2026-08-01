import pytest
import numpy as np

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.smoke, pytest.mark.cghf]

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
