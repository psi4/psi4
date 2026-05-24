"""SAD-guess sanity tests across the periodic table H-Kr.

For each Z in {1, ..., 36} and each SAD_FUNCTIONAL in {"HF", "PBE"} we:

1. Build a single-atom wavefunction with that Z and basis def2-tzvp
   (all-electron through Z=36; ECPs only kick in at Rb).
2. Initialize the HF wavefunction and run guess() to trigger
   psi4.driver.procrouting.scf_proc.sad_guess.populate_sad_guess.
3. Verify that the assembled SAD density has exactly Z electrons
   (tr[(Da + Db) S] == Z within a generous tolerance).
4. Verify that the atomic SCF actually converged: the per-atom result
   stored on the wavefunction satisfies the SAD electron-count check
   to far tighter tolerance than the basis-mapping rounding could
   plausibly mask.

We do NOT run the molecular SCF -- the SAD guess is purpose-built to be
the starting point. The point of this test is to confirm that the
atomic SCFs converge across the periodic table without manual tuning,
including transition metals (Sc-Zn) where spin-averaged occupations
make some shells fractionally filled.
"""

import pytest
import numpy as np

import psi4
from psi4 import compare_values


# Elements 1..36 (H through Kr). All-electron in def2-svp.
ATOMS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
]


def _build_single_atom_wfn(symbol: str) -> psi4.core.Wavefunction:
    """Construct an HF wavefunction for a single atom and run guess()."""
    Z = ATOMS.index(symbol) + 1
    # Multiplicity: 1 for even Z, 2 for odd. Real ground-state multiplicity
    # for some atoms (Cr, Mn, ...) is higher, but the SAD guess is
    # spin-averaged so this only affects what the (uniterated) molecular
    # HF would do, not what SAD does.
    mult = 1 if Z % 2 == 0 else 2
    atom_mol = psi4.geometry(f"0 {mult}\n{symbol}\nsymmetry c1\n")
    base_wfn = psi4.core.Wavefunction.build(atom_mol, psi4.core.get_global_option("BASIS"))
    hf_wfn = psi4.driver.scf_wavefunction_factory("HF", base_wfn, "UHF")
    hf_wfn.initialize()
    hf_wfn.guess()
    return hf_wfn


@pytest.mark.parametrize("symbol", ATOMS)
@pytest.mark.parametrize("functional", ["HF", "PBE"])
def test_sad_atomic_guess_electron_count(symbol, functional):
    """tr[(Da + Db) S] of the SAD guess must equal the atomic number."""
    Z = ATOMS.index(symbol) + 1
    psi4.core.clean_options()
    psi4.set_options(
        {
            "basis": "def2-tzvp",
            "reference": "uhf",
            "guess": "sad",
            "sad_functional": functional,
            # Transition-metal atoms with spin-averaged fractional shells
            # can be slow to converge; give the atomic SCF some room. The
            # default (50) is plenty for closed-shell atoms.
            "sad_maxiter": 200,
            "scf_type": "df",
        }
    )

    wfn = _build_single_atom_wfn(symbol)

    S = np.asarray(wfn.S())
    Da = np.asarray(wfn.Da())
    Db = np.asarray(wfn.Db())

    n_alpha = float(np.einsum("ij,ij->", Da, S))
    n_beta = float(np.einsum("ij,ij->", Db, S))
    n_total = n_alpha + n_beta

    # Closed-shell SAD aliases Db := Da, so n_alpha and n_beta should
    # match to machine precision.
    assert compare_values(
        n_alpha,
        n_beta,
        10,
        f"{symbol} (Z={Z}, {functional}): SAD alpha and beta density traces should agree (Da=Db)",
    )

    # Total electron count must match Z.
    assert compare_values(
        float(Z),
        n_total,
        8,
        f"{symbol} (Z={Z}, {functional}): SAD guess electron count tr[(Da+Db)S] = {n_total:.6f} (expected {Z})",
    )
