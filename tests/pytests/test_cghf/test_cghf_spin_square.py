"""
This test covers the spin_square method. It runs calculations for closed- and
open-shell molecules to determine if spin_square is being computed accurately.

core.CGHF.spin_square(self) returns a tuple of (⟨S²⟩, 2S+1), similar to PySCF.
"""
import pytest

import psi4
from addons import using
import numpy as np

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.cghf, *using("einsums")]

@pytest.fixture(autouse=True)
def clean_psi():
    """Provide a clean workspace for each test."""
    # pre-test cleanup
    psi4.core.clean()
    psi4.core.clean_options()

    yield

    # post-test cleanup
    psi4.core.clean()
    psi4.core.clean_options()


def test_cghf_spin_square_helium():
    """He should be a pure singlet: (0, 1)."""
    mol = psi4.geometry("""
    0 1
      He
    symmetry c1
    """)
    psi4.set_options({
        "basis": "sto-3g",
        "reference": "cghf",
        "scf_type": "direct",
        "df_scf_guess": False,
    })
    _, wfn = psi4.energy("scf", molecule=mol, return_wfn=True)

    ss, multiplicity = wfn.spin_square()

    assert ss == pytest.approx(0.0, abs=1e-8)
    assert multiplicity == pytest.approx(1.0, abs=1e-8)


def test_cghf_spin_square_nitric_oxide():
    """NO should be about a doublet: (.75, 2), but it's contaminated."""
    mol = psi4.geometry("""
    0 2
      N
      O 1 1.165
    symmetry c1
    """)
    psi4.set_options({
        "basis": "6-31G",
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
        "e_convergence": 5e-7
    })
    calc_e, wfn = psi4.energy("scf", molecule=mol, return_wfn=True)

    assert calc_e == pytest.approx(-129.174262,abs=1e-6)

    ss, multiplicity = wfn.spin_square()

    assert ss == pytest.approx(0.932, abs=1e-2)
    assert multiplicity == pytest.approx(2.174, abs=1e-2)


def test_cghf_spin_square_nitrogen():
    """N should be about a quartet: (3.75, 4).
    Starting with a UHF guess until stability analysis is implemented.
    """
    mol = psi4.geometry("""
    0 4
      N
    symmetry c1
    """)
    psi4.set_options({
        "basis": "6-31G",
        "reference": "uhf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
    })
    _, wfn_uhf = psi4.energy("scf", molecule=mol, return_wfn=True)

    Ca = wfn_uhf.Ca().to_array()
    Cb = wfn_uhf.Cb().to_array()
    nso = wfn_uhf.nso()
    nalpha = wfn_uhf.nalpha()
    nbeta = wfn_uhf.nbeta()
    nocc = nalpha + nbeta

    C_guess = np.zeros((2 * nso, nocc), dtype=np.complex128)
    C_guess[:nso, :nalpha] = Ca[:, :nalpha]
    C_guess[nso:, nalpha:] = Cb[:, :nbeta]

    # Save as a ComplexWavefunction checkpoint (bare wfn + set_C is sufficient;
    # from_file / restart_file only needs C and nelec, both available here).
    psi4.core.clean()
    save_wfn = psi4.core.ComplexWavefunction(mol,
        psi4.core.BasisSet.build(mol, "ORBITAL", "6-31G"))
    save_wfn.set_C(psi4.core.ComplexMatrix.from_array(C_guess, name="C"))
    checkpoint = "cghf_nitrogen_uhf_guess"
    save_wfn.to_file(checkpoint)

    psi4.core.clean()
    psi4.set_options({
        "basis": "6-31G",
        "reference": "cghf",
        "scf_type": "direct",
        "df_scf_guess": False,
        "e_convergence": 1e-8,
    })
    calc_e, wfn = psi4.energy("scf", molecule=mol, return_wfn=True, restart_file=checkpoint)

    assert calc_e == pytest.approx(-54.38500771,abs=2e-8)

    ss, multiplicity = wfn.spin_square()

    assert ss == pytest.approx(3.75, abs=1e-2)
    assert multiplicity == pytest.approx(4, abs=1e-2)
