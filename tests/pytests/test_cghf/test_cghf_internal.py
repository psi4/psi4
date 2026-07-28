import pytest
import numpy as np

import psi4
from psi4.driver.procrouting.dft import build_superfunctional

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.cghf]


def _spin_block(A):
    """Embed a real SO matrix into GHF spin-blocked form: diag(A, A)."""
    n = A.shape[0]
    out = np.zeros((2 * n, 2 * n), dtype=np.complex128)
    out[:n, :n] = A
    out[n:, n:] = A
    return out


def _build_cghf(mol, basis_name="sto-3g"):
    """Minimal CGHF instance on a C1 molecule."""
    psi4.set_options({"basis": basis_name, "reference": "cghf"})
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", basis_name)
    ref_wfn = psi4.core.ComplexWavefunction(mol, basis)
    superfunc, _ = build_superfunctional("hf", False)
    return psi4.core.CGHF(ref_wfn, superfunc)


def test_cghf_core_hamiltonian():
    """form_H builds spin-blocked T_, V_, H_ matching real SO kinetic/potential."""
    mol = psi4.geometry(
        """
        0 1
        H
        H 1 0.74
        symmetry c1
        """
    )
    wfn = _build_cghf(mol)
    wfn.form_H()

    mints = wfn.mintshelper()
    T_ref = _spin_block(np.asarray(mints.so_kinetic()))
    V_ref = _spin_block(np.asarray(mints.so_potential()))
    H_ref = T_ref + V_ref

    T = wfn.get_T().to_array()
    V = wfn.get_V().to_array()
    H = wfn.get_H().to_array()

    assert T.shape == T_ref.shape
    assert V.shape == V_ref.shape
    assert H.shape == H_ref.shape
    np.testing.assert_allclose(T, T_ref, atol=1e-12)
    np.testing.assert_allclose(V, V_ref, atol=1e-12)
    np.testing.assert_allclose(H, H_ref, atol=1e-12)
