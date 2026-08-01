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
    opts = {
        "basis": basis_name,
        "reference": "cghf",
        "guess": "core",
        "scf_type": "direct",
        "df_scf_guess": False,
    }
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


def _reference_form_C(F, X):
    """Mirror CGHF::form_C: C = X @ eigvecs(X^H F X)."""
    Fp = X.conj().T @ F @ X
    _, evecs = np.linalg.eigh(Fp)
    return X @ evecs


def _phase_align(C, C_ref):
    """Flip each column phase so <C_ref[:,k]|C[:,k]> is real and positive."""
    C = np.array(C, copy=True, dtype=np.complex128)
    for k in range(C.shape[1]):
        ov = np.vdot(C_ref[:, k], C[:, k])
        if abs(ov) > 1e-14:
            C[:, k] *= np.conj(ov) / abs(ov)
    return C


def test_cghf_form_C():
    """form_C diagonalizes F in the X-orthogonal basis; C^H S C = I."""
    wfn = _build_cghf(_h2_c1(), s_orthogonalization="SYMMETRIC")
    wfn.form_H()
    wfn.form_Shalf()

    # Core-Hamiltonian guess: F <- H
    F_view = wfn.get_F().to_array(copy=False)
    F_view[:] = wfn.get_H().to_array()

    wfn.form_C()

    F = wfn.get_F().to_array()
    X = wfn.get_X().to_array()
    S = wfn.S().to_array()
    C = wfn.get_C().to_array()
    C_ref = _reference_form_C(F, X)

    assert C.shape == C_ref.shape
    np.testing.assert_allclose(C.conj().T @ S @ C, np.eye(C.shape[1]), atol=1e-10)

    # MO Fock matrix should be diagonal (eigenvalues of X^H F X)
    F_mo = C.conj().T @ F @ C
    np.testing.assert_allclose(F_mo, np.diag(np.diag(F_mo)), atol=1e-10)

    C_aligned = _phase_align(C, C_ref)
    np.testing.assert_allclose(C_aligned, C_ref, atol=1e-10)

    with pytest.raises(RuntimeError, match="Level shifting"):
        wfn.form_C(shift=0.1)


def test_cghf_compute_initial_E():
    """compute_initial_E returns Enuc + Re(Tr(H D))."""
    wfn = _build_cghf(_h2_c1())
    wfn.form_H()

    H = wfn.get_H().to_array()
    n = H.shape[0]
    rng = np.random.default_rng(0)
    # Hermitian density so Tr(H D) is real for Hermitian H
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    D_ref = A @ A.conj().T

    D_view = wfn.get_D().to_array(copy=False)
    D_view[:] = D_ref

    enuc = wfn.molecule().nuclear_repulsion_energy()
    E_ref = enuc + np.real(np.trace(H @ D_ref))
    E = wfn.compute_initial_E()

    assert E == pytest.approx(E_ref, abs=1e-12)


def _core_guess_cghf(basis_name="sto-3g", **scf_options):
    """CGHF with core-Hamiltonian orbitals (form_H --> form_Shalf --> form_C)."""
    wfn = _build_cghf(_h2_c1(), basis_name=basis_name, s_orthogonalization="SYMMETRIC", **scf_options)
    wfn.form_H()
    wfn.form_Shalf()
    F_view = wfn.get_F().to_array(copy=False)
    F_view[:] = wfn.get_H().to_array()
    wfn.form_C()
    return wfn


def _jk_reference_ao(basis, D):
    """J_pq = sum_rs D_rs (pq|rs), K_pr = sum_qs D_qs (pq|rs) via ao_eri."""
    I = np.asarray(psi4.core.MintsHelper(basis).ao_eri())
    J = np.einsum("pqrs,rs->pq", I, D, optimize=True)
    K = np.einsum("pqrs,qs->pr", I, D, optimize=True)
    return J, K


def test_cghf_form_D():
    """form_D builds D = C_occ @ C_occ^H with Tr(D S) = nelectron."""
    wfn = _core_guess_cghf()
    C = wfn.get_C().to_array()
    S = wfn.S().to_array()
    nocc = int(sum(wfn.nelecpi()))

    wfn.form_D()
    D = wfn.get_D().to_array()
    D_ref = C[:, :nocc] @ C[:, :nocc].conj().T

    assert D.shape == (C.shape[0], C.shape[0])
    np.testing.assert_allclose(D, D_ref, atol=1e-12)
    np.testing.assert_allclose(D, D.conj().T, atol=1e-12)  # Hermitian
    assert np.trace(D @ S).real == pytest.approx(nocc, abs=1e-10)
    # Virtuals must not contribute
    assert not np.allclose(D, C @ C.conj().T)


def _ghf_jk_reference(basis, D_full):
    """Reference spin-blocked J/K for a GHF density, built from plain spatial AO integrals."""
    nbf = basis.nbf()
    D_aa = D_full[:nbf, :nbf]
    D_bb = D_full[nbf:, nbf:]
    D_ab = D_full[:nbf, nbf:]
    D_ba = D_full[nbf:, :nbf]

    J_tot, _ = _jk_reference_ao(basis, D_aa + D_bb)
    _, K_aa = _jk_reference_ao(basis, D_aa)
    _, K_bb = _jk_reference_ao(basis, D_bb)
    _, K_ab = _jk_reference_ao(basis, D_ab)
    _, K_ba = _jk_reference_ao(basis, D_ba)

    J = np.zeros_like(D_full)
    K = np.zeros_like(D_full)
    J[:nbf, :nbf] = J_tot
    J[nbf:, nbf:] = J_tot
    K[:nbf, :nbf] = K_aa
    K[nbf:, nbf:] = K_bb
    K[:nbf, nbf:] = K_ab
    K[nbf:, :nbf] = K_ba
    return J, K


def _full_guess_cghf(basis_name="sto-3g", **scf_options):
    """CGHF after core guess + one real ComplexDirectJK build: H_, S_/X_, C_, D_, J_/K_/G_ all populated."""
    wfn = _core_guess_cghf(basis_name=basis_name, **scf_options)
    wfn.form_D()

    jk = psi4.core.ComplexJK.build_JK(wfn.basisset(), wfn.basisset())
    jk.initialize()
    wfn.set_jk(jk)
    wfn.form_G()
    return wfn


def test_cghf_compute_E():
    wfn = _full_guess_cghf(scf_type="direct", screening="NONE")

    H = wfn.get_H().to_array()
    T = wfn.get_T().to_array()
    D = wfn.get_D().to_array()
    J = wfn.get_J().to_array()
    K = wfn.get_K().to_array()
    enuc = wfn.molecule().nuclear_repulsion_energy()

    one_electron_E = np.trace(H @ D).real
    kinetic_E = np.trace(T @ D).real
    two_electron_E = 0.5 * (np.trace(J @ D).real - np.trace(K @ D).real)
    E_ref = enuc + one_electron_E + two_electron_E

    E = wfn.compute_E()

    assert E == pytest.approx(E_ref, abs=1e-10)
    assert wfn.get_energies("Nuclear") == pytest.approx(enuc, abs=1e-12)
    assert wfn.get_energies("Kinetic") == pytest.approx(kinetic_E, abs=1e-10)
    assert wfn.get_energies("One-Electron") == pytest.approx(one_electron_E, abs=1e-10)
    assert wfn.get_energies("Two-Electron") == pytest.approx(two_electron_E, abs=1e-10)


