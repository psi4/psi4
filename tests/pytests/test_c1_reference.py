import numpy as np
import pytest

import psi4
from psi4.driver.procrouting import proc_util


pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def test_prepare_c1_reference_preserves_scf_quantities():
    """A post-SCF symmetry cast preserves the AO reference quantities."""

    h2o = psi4.geometry(
        """
        0 1
        O
        H 1 0.957
        H 1 0.957 2 104.5
        """
    )
    psi4.set_options({"basis": "sto-3g", "scf_type": "pk", "reference": "rhf"})

    _, sym_wfn = psi4.energy("scf", molecule=h2o, return_wfn=True)
    assert sym_wfn.nirrep() > 1

    expected_c = np.asarray(sym_wfn.Ca_subset("AO", "ALL")).copy()
    expected_f = np.asarray(sym_wfn.Fa_subset("AO")).copy()
    expected_eps = np.asarray(sym_wfn.epsilon_a_subset("AO", "ALL")).copy()
    expected_s = np.asarray(psi4.core.MintsHelper(sym_wfn.basisset()).ao_overlap()).copy()

    c1_wfn = proc_util.prepare_c1_reference(sym_wfn)

    # Conversion is non-destructive, while an already-C1 wavefunction is a no-op.
    assert sym_wfn.nirrep() > 1
    assert c1_wfn.nirrep() == 1
    assert c1_wfn.molecule().schoenflies_symbol() == "c1"
    assert proc_util.prepare_c1_reference(c1_wfn) is c1_wfn

    np.testing.assert_allclose(np.asarray(c1_wfn.Ca()), expected_c, atol=1.0e-12)
    np.testing.assert_allclose(np.asarray(c1_wfn.S()), expected_s, atol=1.0e-12)
    np.testing.assert_allclose(np.asarray(c1_wfn.Fa()), expected_f, atol=1.0e-12)
    np.testing.assert_allclose(np.asarray(c1_wfn.epsilon_a()), expected_eps, atol=1.0e-12)
    np.testing.assert_allclose(c1_wfn.energy(), sym_wfn.energy(), atol=1.0e-12)

    # The transformed canonical orbitals still solve FC = SC epsilon in C1.
    c = np.asarray(c1_wfn.Ca())
    s = np.asarray(c1_wfn.S())
    f = np.asarray(c1_wfn.Fa())
    eps = np.asarray(c1_wfn.epsilon_a())
    np.testing.assert_allclose(f @ c, (s @ c) * eps[np.newaxis, :], atol=1.0e-10)
