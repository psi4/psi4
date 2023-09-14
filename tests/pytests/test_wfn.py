import pytest

import numpy as np
import psi4

from utils import compare_arrays

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

def test_fock_subset_mo():
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    rhf_e, wfn = psi4.energy('SCF/cc-pVDZ', molecule=h2o, return_wfn=True)

    F_diagonals = []
    for h in wfn.epsilon_a().nph:
        F_diagonals.append(np.diag(h))
    F_expected = psi4.core.Matrix.from_array(F_diagonals)
    assert psi4.compare_matrices(F_expected, wfn.Fa_subset("MO"), 8, "Alpha Fock Matrix")
    assert psi4.compare_matrices(F_expected, wfn.Fb_subset("MO"), 8, "Beta Fock Matrix")
