import pytest
import numpy as np

import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.cghf]


def _random_complex(shape, rng):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def _build_complex_jk():
    """Minimal ComplexDirectJK via ComplexJK.build_JK (SCF_TYPE=DIRECT)."""
    mol = psi4.geometry(
        """
        0 1
        H
        H 1 0.74
        symmetry c1
        """
    )
    psi4.set_options({"SCF_TYPE": "DIRECT", "SCREENING": "NONE", "DF_SCF_GUESS": False})
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "sto-3g")
    return psi4.core.ComplexJK.build_JK(basis, basis)


def _run_compute_D(C_left, C_right=None):
    """Queue one or more C pairs and form D. C_* may be a matrix or a sequence of matrices."""
    if C_right is None:
        C_right = C_left

    lefts = C_left if isinstance(C_left, (list, tuple)) else [C_left]
    rights = C_right if isinstance(C_right, (list, tuple)) else [C_right]
    assert len(lefts) == len(rights)

    jk = _build_complex_jk()
    jk.C_clear()
    for Cl, Cr in zip(lefts, rights):
        jk.C_left_add(Cl)
        jk.C_right_add(Cr)
    jk.compute_D()
    assert len(jk.D()) == len(lefts)
    return jk.D()


def test_complexjk_compute_D_c1_rectangular():
    """C1 non-square C: D = C @ C.conj().T."""
    rng = np.random.default_rng(0)
    C_arr = _random_complex((5, 2), rng)
    C = psi4.core.ComplexMatrix.from_array(C_arr, name="C")

    D = _run_compute_D(C)[0]
    assert D.num_blocks() == 1
    np.testing.assert_allclose(D.to_array(), C_arr @ C_arr.conj().T)


def test_complexjk_compute_D_c1_asymmetric():
    """C1 non-square C_left != C_right: D = C_left @ C_right.conj().T."""
    rng = np.random.default_rng(2)
    Cl_arr = _random_complex((5, 2), rng)
    Cr_arr = _random_complex((5, 2), rng)
    Cl = psi4.core.ComplexMatrix.from_array(Cl_arr, name="Cl")
    Cr = psi4.core.ComplexMatrix.from_array(Cr_arr, name="Cr")

    D = _run_compute_D(Cl, Cr)[0]
    assert D.num_blocks() == 1
    np.testing.assert_allclose(D.to_array(), Cl_arr @ Cr_arr.conj().T)
    # Sanity: not the same as the symmetric products
    assert not np.allclose(D.to_array(), Cl_arr @ Cl_arr.conj().T)
    assert not np.allclose(D.to_array(), Cr_arr @ Cr_arr.conj().T)


def test_complexjk_compute_D_c1_rectangular_D():
    """C_left/C_right with different AO counts yield rectangular D = Cl @ Cr.H."""
    rng = np.random.default_rng(4)
    Cl_arr = _random_complex((5, 2), rng)  # 5 AOs × 2 occ
    Cr_arr = _random_complex((3, 2), rng)  # 3 AOs × 2 occ
    Cl = psi4.core.ComplexMatrix.from_array(Cl_arr, name="Cl")
    Cr = psi4.core.ComplexMatrix.from_array(Cr_arr, name="Cr")

    D = _run_compute_D(Cl, Cr)[0]
    assert D.num_blocks() == 1
    D_arr = D.to_array()
    assert D_arr.shape == (5, 3)
    np.testing.assert_allclose(D_arr, Cl_arr @ Cr_arr.conj().T)


def test_complexjk_compute_D_multi_irrep_rectangular():
    """Multi-irrep non-square C: D_h = C_h @ C_h.conj().T per diagonal tile."""
    rng = np.random.default_rng(1)
    row_sizes = [4, 2, 3]
    col_sizes = [1, 2, 1]
    C = psi4.core.ComplexMatrix("C", row_sizes, col_sizes)

    refs = []
    for h, view in enumerate(C.array_interface()):
        C_h = _random_complex(view.shape, rng)
        view[:] = C_h
        refs.append(C_h @ C_h.conj().T)

    D = _run_compute_D(C)[0]
    assert D.num_blocks() == len(row_sizes)
    D_blocks = D.to_array()
    assert isinstance(D_blocks, list)
    assert len(D_blocks) == len(row_sizes)
    for got, ref, nso in zip(D_blocks, refs, row_sizes):
        assert got.shape == (nso, nso)
        np.testing.assert_allclose(got, ref)


def test_complexjk_compute_D_multiple_C_pairs():
    """Multiple queued C_left/C_right pairs each produce their own D = Cl @ Cr.H."""
    rng = np.random.default_rng(3)
    shapes = [(4, 1), (5, 2), (3, 2)]
    lefts, rights, refs = [], [], []
    for i, shape in enumerate(shapes):
        Cl_arr = _random_complex(shape, rng)
        Cr_arr = _random_complex(shape, rng)
        lefts.append(psi4.core.ComplexMatrix.from_array(Cl_arr, name=f"Cl{i}"))
        rights.append(psi4.core.ComplexMatrix.from_array(Cr_arr, name=f"Cr{i}"))
        refs.append(Cl_arr @ Cr_arr.conj().T)

    Ds = _run_compute_D(lefts, rights)
    assert len(Ds) == len(shapes)
    for D, ref in zip(Ds, refs):
        assert D.num_blocks() == 1
        np.testing.assert_allclose(D.to_array(), ref)


def _h2_sto3g():
    """C1 H2 / STO-3G molecule and orbital basis."""
    mol = psi4.geometry(
        """
        0 1
        H
        H 1 0.74
        symmetry c1
        """
    )
    psi4.set_options({"SCF_TYPE": "DIRECT", "SCREENING": "NONE", "DF_SCF_GUESS": False})
    basis = psi4.core.BasisSet.build(mol, "ORBITAL", "sto-3g")
    return mol, basis


def _jk_reference(basis, D):
    """J_pq = sum_rs D_rs (pq|rs), K_pr = sum_qs D_qs (pq|rs) via MintsHelper.ao_eri."""
    mints = psi4.core.MintsHelper(basis)
    I = np.asarray(mints.ao_eri())
    J = np.einsum("pqrs,rs->pq", I, D, optimize=True)
    K = np.einsum("pqrs,qs->pr", I, D, optimize=True)
    return J, K


def test_complexdirectjk_matches_einsum():
    """ComplexDirectJK J/K match explicit einsum contractions of ao_eri with complex D."""
    _, basis = _h2_sto3g()
    nbf = basis.nbf()
    rng = np.random.default_rng(7)
    C_arr = _random_complex((nbf, 1), rng)
    D_ref = C_arr @ C_arr.conj().T
    J_ref, K_ref = _jk_reference(basis, D_ref)

    jk = psi4.core.ComplexJK.build_JK(basis, basis)
    jk.initialize()
    jk.C_clear()
    jk.C_add(psi4.core.ComplexMatrix.from_array(C_arr, name="C"))
    jk.compute()
    jk.finalize()

    np.testing.assert_allclose(jk.D()[0].to_array(), D_ref)
    np.testing.assert_allclose(jk.J()[0].to_array(), J_ref, atol=1e-10)
    np.testing.assert_allclose(jk.K()[0].to_array(), K_ref, atol=1e-10)


def test_complexdirectjk_asymmetric_matches_einsum():
    """ComplexDirectJK with C_left != C_right still matches einsum (general complex D)."""
    _, basis = _h2_sto3g()
    nbf = basis.nbf()
    rng = np.random.default_rng(11)
    Cl_arr = _random_complex((nbf, 1), rng)
    Cr_arr = _random_complex((nbf, 1), rng)
    D_ref = Cl_arr @ Cr_arr.conj().T
    J_ref, K_ref = _jk_reference(basis, D_ref)

    jk = psi4.core.ComplexJK.build_JK(basis, basis)
    jk.initialize()
    jk.C_clear()
    jk.C_left_add(psi4.core.ComplexMatrix.from_array(Cl_arr, name="Cl"))
    jk.C_right_add(psi4.core.ComplexMatrix.from_array(Cr_arr, name="Cr"))
    jk.compute()
    jk.finalize()

    np.testing.assert_allclose(jk.D()[0].to_array(), D_ref)
    np.testing.assert_allclose(jk.J()[0].to_array(), J_ref, atol=1e-10)
    np.testing.assert_allclose(jk.K()[0].to_array(), K_ref, atol=1e-10)


def test_complexdirectjk_real_D_matches_jk():
    """Completely real D: ComplexDirectJK agrees with real Direct JK."""
    _, basis = _h2_sto3g()
    nbf = basis.nbf()
    rng = np.random.default_rng(13)
    C_arr = rng.normal(size=(nbf, 1))

    cjk = psi4.core.ComplexJK.build_JK(basis, basis)
    cjk.initialize()
    cjk.C_clear()
    cjk.C_add(psi4.core.ComplexMatrix.from_array(C_arr.astype(np.complex128), name="C"))
    cjk.compute()
    cjk.finalize()

    rjk = psi4.core.JK.build_JK(basis, basis)
    rjk.initialize()
    rjk.C_clear()
    rjk.C_add(psi4.core.Matrix.from_array(C_arr, name="C"))
    rjk.compute()
    rjk.finalize()

    np.testing.assert_allclose(cjk.D()[0].to_array().real, np.asarray(rjk.D()[0]), atol=1e-10)
    np.testing.assert_allclose(cjk.J()[0].to_array().real, np.asarray(rjk.J()[0]), atol=1e-10)
    np.testing.assert_allclose(cjk.K()[0].to_array().real, np.asarray(rjk.K()[0]), atol=1e-10)
    np.testing.assert_allclose(cjk.J()[0].to_array().imag, 0.0, atol=1e-10)
    np.testing.assert_allclose(cjk.K()[0].to_array().imag, 0.0, atol=1e-10)
