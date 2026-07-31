import pytest
import numpy as np

import psi4
psi4.core.be_quiet()

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.smoke, pytest.mark.cghf]

def test_cghf_basic_scf():
    """Minimal CGHF full SCF example."""
    mol = psi4.geometry("""
    0 1
        O 0 0 0
        O 0 0 1.5
    symmetry c1
    """)
    psi4.set_options({
        "basis": "cc-pVDZ",
        "reference": "cghf",
        "guess": "core",
        "print": 0,
        "scf_type": "direct",
        "diis": False,
        "scf_initial_accelerator": "none",
        "df_scf_guess": False
    })
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", psi4.core.get_global_option("BASIS"))
    e = psi4.energy("scf", molecule=mol, scf_do_properties=False)
    assert type(e) == float
    # TODO: get reference energy


