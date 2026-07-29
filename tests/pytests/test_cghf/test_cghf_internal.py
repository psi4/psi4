import pytest
import numpy as np
from scipy.linalg.lapack import dpstrf

import psi4
from psi4.driver.procrouting.dft import build_superfunctional

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.cghf]

_ORTHOG_METHODS = ["SYMMETRIC", "CANONICAL", "PARTIALCHOLESKY", "AUTO"]


def _h2_c1():
    return psi4.geometry(
        """
        0 1
        H
        H 1 0.74
        symmetry c1
        """
    )


def _spin_block(A):
    """Embed a real SO matrix into GHF spin-blocked form: diag(A, A). Supports rectangular A."""
    A = np.asarray(A)
    n0, n1 = A.shape
    out = np.zeros((2 * n0, 2 * n1), dtype=np.complex128)
    out[:n0, :n1] = A
    out[n0:, n1:] = A
    return out


def _build_cghf(mol, basis_name="sto-3g", **scf_options):
    """Minimal CGHF instance on a C1 molecule."""
    opts = {"basis": basis_name, "reference": "cghf"}
    opts.update(scf_options)
    psi4.set_options(opts)
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", basis_name)
    ref_wfn = psi4.core.ComplexWavefunction(mol, basis)
    superfunc, _ = build_superfunctional("hf", False)
    return psi4.core.CGHF(ref_wfn, superfunc)


def _build_rhf(mol, basis_name="sto-3g", **scf_options):
    """Minimal RHF instance used for nmo / orthog cross-checks."""
    opts = {"basis": basis_name, "reference": "rhf"}
    opts.update(scf_options)
    psi4.set_options(opts)
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", basis_name)
    ref_wfn = psi4.core.Wavefunction(mol, basis)
    superfunc, _ = build_superfunctional("hf", True)
    return psi4.core.RHF(ref_wfn, superfunc)


def _reference_orthog_X(S, method, s_tol=1.0e-7, chol_tol=1.0e-8):
    """Mirror BasisSetOrthogonalization (normalize → orthog → unroll) for a C1 overlap."""
    S = np.asarray(S, dtype=np.float64)
    nbf = S.shape[0]
    norm = 1.0 / np.sqrt(np.diag(S))
    Sn = (norm[:, None] * S) * norm[None, :]
    evals, evecs = np.linalg.eigh(Sn)
    min_S = float(evals.min())
    rcond = float(evals.min() / evals.max()) if evals.max() != 0.0 else 0.0

    m = method.upper()
    if m == "AUTO":
        if rcond <= np.finfo(float).eps or not np.isfinite(rcond):
            m = "PARTIALCHOLESKY"
        elif min_S < s_tol:
            m = "CANONICAL"
        else:
            m = "SYMMETRIC"

    if m == "SYMMETRIC":
        Us = evecs * (1.0 / np.sqrt(evals))
        Xn = Us @ evecs.T
    elif m == "CANONICAL":
        keep = evals > s_tol
        Xn = evecs[:, keep] * (1.0 / np.sqrt(evals[keep]))
    elif m == "PARTIALCHOLESKY":
        # Match BasisSetOrthogonalization::compute_partial_cholesky_orthog
        od = np.sum(np.abs(Sn), axis=1) - np.abs(np.diag(Sn))
        order = np.argsort(od, kind="stable")
        Stmp = Sn[np.ix_(order, order)]
        _, piv, rank, info = dpstrf(Stmp, tol=chol_tol, lower=True)
        if info < 0:
            raise RuntimeError(f"dpstrf failed with info={info}")
        pivots = order[(piv[:rank] - 1).astype(int)]
        Ssub = Sn[np.ix_(pivots, pivots)]
        ev, U = np.linalg.eigh(Ssub)
        keep = ev > s_tol
        Xsub = U[:, keep] * (1.0 / np.sqrt(ev[keep]))
        Xn = np.zeros((nbf, Xsub.shape[1]), dtype=np.float64)
        Xn[pivots, :] = Xsub
    else:
        raise ValueError(f"Unrecognized orthogonalization method: {method}")

    return norm[:, None] * Xn


def _assert_form_Shalf(wfn, method, basis_name, s_tol=1.0e-7, chol_tol=1.0e-8):
    """Check spin-blocked S_/X_ from form_Shalf against SO overlap + orthog reference."""
    wfn.form_Shalf()

    S_real = np.asarray(wfn.mintshelper().so_overlap())
    S = wfn.S().to_array()
    X = wfn.get_X().to_array()
    S_ref = _spin_block(S_real)
    X_ref_real = _reference_orthog_X(S_real, method, s_tol=s_tol, chol_tol=chol_tol)
    X_ref = _spin_block(X_ref_real)

    np.testing.assert_allclose(S, S_ref, atol=1e-12)
    assert X.shape[0] == S.shape[0]

    # X^H S X = I in the orthogonal basis
    XtSX = X.conj().T @ S @ X
    np.testing.assert_allclose(XtSX, np.eye(X.shape[1]), atol=1e-10)

    # Same independent-orbital count as real RHF (doubled by spin blocking)
    rhf = _build_rhf(
        wfn.molecule(),
        basis_name,
        s_orthogonalization=method,
        s_tolerance=s_tol,
        s_cholesky_tolerance=chol_tol,
    )
    rhf.form_Shalf()
    assert X.shape[1] == 2 * rhf.nmo()

    # SYMMETRIC X is unique. CANONICAL/AUTO: compare phase-invariant XX^H.
    # PARTIALCHOLESKY pivots can differ across LAPACK builds, so skip value match.
    m = method.upper()
    if m == "SYMMETRIC":
        assert X.shape == X_ref.shape
        np.testing.assert_allclose(X, X_ref, atol=1e-10)
    elif m in ("CANONICAL", "AUTO"):
        assert X.shape == X_ref.shape
        np.testing.assert_allclose(X @ X.conj().T, X_ref @ X_ref.conj().T, atol=1e-8)

def test_cghf_core_hamiltonian():
    """form_H builds spin-blocked T_, V_, H_ matching real SO kinetic/potential."""
    wfn = _build_cghf(_h2_c1())
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


@pytest.mark.parametrize("method", _ORTHOG_METHODS)
def test_cghf_form_Shalf(method):
    """form_Shalf builds spin-blocked S_ and X_ for each S_ORTHOGONALIZATION option."""
    wfn = _build_cghf(_h2_c1(), s_orthogonalization=method)
    _assert_form_Shalf(wfn, method, "sto-3g")


def test_cghf_form_Shalf_cc_pvtz():
    """One larger-basis form_Shalf check (cc-pVTZ, default AUTO orthogonalization)."""
    wfn = _build_cghf(_h2_c1(), basis_name="cc-pVTZ", s_orthogonalization="AUTO")
    _assert_form_Shalf(wfn, "AUTO", "cc-pVTZ", s_tol=1e-3, chol_tol=1e-3)

