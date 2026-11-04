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


def test_module_roles():
    """set_module/module are the "qc_module" entry of the role map, and other roles ride along."""
    h2o = psi4.geometry("""
        O
        H 1 1.0
        H 1 1.0 2 101.5
    """)

    _, wfn = psi4.energy('SCF/cc-pVDZ', molecule=h2o, return_wfn=True)

    # set_module/module are a view on one role, so the two agree and neither loses the other
    assert wfn.module() == "scf"
    assert wfn.module_roles()["qc_module"] == "scf"
    assert wfn.module_role("qc_module") == "scf"

    # an unfilled role is empty rather than a lookup error
    assert wfn.module_role("no_such_role") == ""

    wfn.set_module_role("widget_maker", "acme")
    assert wfn.module_role("widget_maker") == "acme"
    assert wfn.module() == "scf"

    wfn.set_module("dfmp2")
    assert wfn.module() == "dfmp2"
    assert wfn.module_role("widget_maker") == "acme"

    # roles survive a file round trip, which carries them flattened under a prefix
    reloaded = psi4.core.Wavefunction.from_file(wfn.to_file())
    assert reloaded.module() == "dfmp2"
    assert reloaded.module_role("widget_maker") == "acme"
