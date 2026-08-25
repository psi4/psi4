import pytest
import numpy as np

import psi4
from addons import using
from psi4.driver.p4util.exceptions import ValidationError

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, *using("einsums")]


def test_complexmatrix_constructor_zeroed():
    """ComplexMatrix(name, block_sizes) builds a zeroed multi-block tensor."""
    mat = psi4.core.ComplexMatrix("zeros", [2, 3, 1])
    assert mat.nirrep() == 3

    views = mat.array_interface()
    assert len(views) == 3
    assert [v.shape for v in views] == [(2, 2), (3, 3), (1, 1)]
    assert all(v.dtype == np.complex128 for v in views)

    blocks = mat.to_array()
    assert isinstance(blocks, list)
    for blk in blocks:
        np.testing.assert_allclose(blk, np.zeros_like(blk))


def test_complexmatrix_constructor_rectangular():
    """ComplexMatrix(name, row_sizes, col_sizes) builds rectangular diagonal tiles."""
    mat = psi4.core.ComplexMatrix("Cocc", [4, 2], [1, 3])
    assert mat.nirrep() == 2
    views = mat.array_interface()
    assert [v.shape for v in views] == [(4, 1), (2, 3)]
    views[0][:] = [[1 + 1j], [2], [3], [4]]
    views[1][:] = np.arange(6, dtype=np.complex128).reshape(2, 3)
    blocks = mat.to_array()
    np.testing.assert_allclose(blocks[0], views[0])
    np.testing.assert_allclose(blocks[1], views[1])


def test_complexmatrix_array_interface_write_through():
    """array_interface views share memory with the ComplexMatrix."""
    mat = psi4.core.ComplexMatrix("view", [2])
    view = mat.array_interface()[0]
    view[:] = [[1 + 2j, 3], [4, 5 - 1j]]
    np.testing.assert_allclose(mat.to_array(), view)

    mat.to_array(copy=False)[0, 1] = 9 + 9j
    assert view[0, 1] == pytest.approx(9 + 9j)


def test_complexmatrix_to_from_array_c1():
    """C1 ComplexMatrix NumPy round-trip via to_array / from_array."""
    n = 3
    ref = np.arange(n * n, dtype=np.complex128).reshape(n, n)
    ref += 1.0j * np.flip(ref)

    mat = psi4.core.ComplexMatrix.from_array(ref, name="Test")
    assert mat.nirrep() == 1
    out = mat.to_array()
    assert isinstance(out, np.ndarray)
    np.testing.assert_allclose(out, ref)

    # Real input is cast to complex128
    real = np.eye(2)
    cast = psi4.core.ComplexMatrix.from_array(real)
    assert cast.to_array().dtype == np.complex128
    np.testing.assert_allclose(cast.to_array(), real)

    # Rectangular C1 tiles are supported (TiledTensor)
    rect = np.arange(6, dtype=np.complex128).reshape(3, 2)
    rect_mat = psi4.core.ComplexMatrix.from_array(rect, name="Cocc")
    np.testing.assert_allclose(rect_mat.to_array(), rect)
    assert rect_mat.array_interface()[0].shape == (3, 2)

    with pytest.raises(ValidationError, match="2-dimensional"):
        psi4.core.ComplexMatrix.from_array(np.ones(3, dtype=np.complex128))


def test_complexmatrix_to_array_copy_flag():
    """to_array(copy=True) isolates; copy=False aliases C1 and multi-block tensors."""
    mat = psi4.core.ComplexMatrix.from_array(np.arange(4, dtype=np.complex128).reshape(2, 2))

    copied = mat.to_array(copy=True)
    copied[0, 0] = -1
    assert mat.to_array()[0, 0] == pytest.approx(0)

    aliased = mat.to_array(copy=False)
    aliased[0, 0] = 42
    assert mat.to_array()[0, 0] == pytest.approx(42)

    multi = psi4.core.ComplexMatrix("multi", [2, 1])
    multi.array_interface()[0][:] = np.eye(2, dtype=np.complex128)
    multi.array_interface()[1][:] = [[3 + 1j]]
    blocks = multi.to_array(copy=False)
    assert isinstance(blocks, list)
    blocks[1][0, 0] = 8
    assert multi.to_array()[1][0, 0] == pytest.approx(8)


def test_complexmatrix_to_array_multiblock():
    """Multi-block ComplexMatrix.to_array returns a list of per-block arrays."""
    sizes = [2, 3, 1]
    mat = psi4.core.ComplexMatrix("mb", sizes)
    refs = []
    for h, n in enumerate(sizes):
        ref = np.arange(n * n, dtype=np.complex128).reshape(n, n) + 1j * (h + 1)
        mat.array_interface()[h][:] = ref
        refs.append(ref)

    blocks = mat.to_array()
    assert isinstance(blocks, list)
    assert len(blocks) == len(sizes)
    for got, ref in zip(blocks, refs):
        np.testing.assert_allclose(got, ref)


# ---------------------------------------------------------------------------
#   GEMM operations
# ---------------------------------------------------------------------------


def test_doublet_conjugation():
    """Verify that doublet(A, Id, True, False) computes A.conjT()."""
    rng = np.random.default_rng(19)
    n = 4

    A_np = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    I_np = np.eye(n, dtype=np.complex128)

    A_cm = psi4.core.ComplexMatrix.from_array(A_np, name="A")
    I_cm = psi4.core.ComplexMatrix.from_array(I_np, name="I")

    result = psi4.core.doublet(A_cm, I_cm, True, False).to_array()

    A_H = A_np.conj().T

    assert np.allclose(result, A_H, atol=1e-14)


def test_doublet_atypical_shapes():
    """Test linalg::doublet edge cases where a tile is zero or not allocated."""
    A_cm = psi4.core.ComplexMatrix("A", [1, 4, 3], [4, 0, 2])
    B_cm = psi4.core.ComplexMatrix("B", [1, 4, 3], [3, 0, 5])

    result = psi4.core.doublet(A_cm, B_cm, True, False).to_array()

    assert result[0].shape == (4, 3)
    assert result[1].shape == (0, 0)
    assert result[2].shape == (2, 5)


# ---------------------------------------------------------------------------
#  diagonalize
# ---------------------------------------------------------------------------


def _inv_sqrt_hermitian(M):
    """Inverse matrix square root of a Hermitian positive-definite matrix M.
    By the spectral theorem, this matrix exists and is unique.
    Given  M := U @ diag(s) @ U^H,   M^{-1/2} = U @ diag(1/sqrt(s)) @ U^H"""
    s, U = np.linalg.eigh(M)
    assert np.all(s > 0), "M must be positive-definite"
    return U @ np.diag(1.0 / np.sqrt(s)) @ U.conj().T


def test_complexmatrix_diagonalize_hermitian():
    """diagonalize(F, X) diagonalizes X^H @ F @ X.

    The function returns (evals, U^H) where U^H has eigenvectors as columns
    and satisfies U^H @ (X^H @ F @ X) @ U = diag(evals).

    S is kept real-symmetric (like the physical overlap matrix).
    """
    rng = np.random.default_rng(42)
    n = 4

    # Hermitian F
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    F_np = A + A.conj().T

    # Real-symmetric positive-definite S, then orthogonalizer X = S^{-1/2}
    B = rng.normal(size=(n, n))
    S_np = B @ B.T
    X_np = _inv_sqrt_hermitian(S_np)

    F_cm = psi4.core.ComplexMatrix.from_array(F_np, name="F")
    X_cm = psi4.core.ComplexMatrix.from_array(X_np, name="X")

    evals_vec, evecs_cm = psi4.core.diagonalize(F_cm, X_cm)

    evals_np = np.asarray(evals_vec.array_interface()[0])
    evecs_np = evecs_cm.to_array()  # U^H (eigenvectors as columns)

    # Forth = X^H @ F @ X; verify  U^H @ Forth @ U = diag(evals)
    Forth_np = X_np @ F_np @ X_np  # X is real-symmetric, so X^H = X^T = X
    residual = evecs_np.conj().T @ Forth_np @ evecs_np - np.diag(evals_np)
    np.testing.assert_allclose(residual, 0.0, atol=1e-10)


def test_complexmatrix_diagonalize_multi_block():
    """diagonalize works across multiple irrep blocks."""
    rng = np.random.default_rng(99)
    sizes = [2, 3]

    F_cm = psi4.core.ComplexMatrix("F", sizes)
    X_cm = psi4.core.ComplexMatrix("X", sizes)

    for h, n in enumerate(sizes):
        # Hermitian F per block
        A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        F_blk = A + A.conj().T
        F_cm.array_interface()[h][:] = F_blk

        # Real-symmetric positive-definite S -> orthogonalizer X = S^{-1/2}
        B = rng.normal(size=(n, n))
        S_blk = B @ B.T
        X_cm.array_interface()[h][:] = _inv_sqrt_hermitian(S_blk)

    evals_vec, evecs_cm = psi4.core.diagonalize(F_cm, X_cm)

    evals_arrs = evals_vec.array_interface()
    evecs_blocks = evecs_cm.to_array()  # list of U_h^H per block
    for h in range(len(sizes)):
        evals_np_blk = np.asarray(evals_arrs[h])
        evecs_blk = evecs_blocks[h]  # U^H
        F_blk = F_cm.to_array()[h]
        X_blk = X_cm.to_array()[h]  # real-symmetric, so X^H = X

        Forth_blk = X_blk @ F_blk @ X_blk
        residual = evecs_blk.conj().T @ Forth_blk @ evecs_blk - np.diag(evals_np_blk)
        np.testing.assert_allclose(residual, 0.0, atol=1e-10)


# ---------------------------------------------------------------------------
#  axpy
# ---------------------------------------------------------------------------


def test_complexmatrix_axpy():
    """axpy: in-place self += alpha * other with complex alpha."""
    rng = np.random.default_rng(42)
    n = 4
    A_ref = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    B_ref = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))

    mat = psi4.core.ComplexMatrix.from_array(A_ref.copy(), name="A")
    other = psi4.core.ComplexMatrix.from_array(B_ref, name="B")

    alpha = 2.0 + 3.0j
    mat.axpy(alpha, other)
    np.testing.assert_allclose(mat.to_array(), A_ref + alpha * B_ref, atol=1e-14)


def test_complexmatrix_axpy_multi_block():
    """axpy works across multiple diagonal tiles."""
    sizes = [2, 3]
    mat = psi4.core.ComplexMatrix("A", sizes)
    other = psi4.core.ComplexMatrix("B", sizes)

    ref_blocks = []
    for h, n in enumerate(sizes):
        mat.array_interface()[h][:] = np.full((n, n), (h + 1) + 1j * (h + 2), dtype=np.complex128)
        other.array_interface()[h][:] = np.eye(n, dtype=np.complex128) * (0.5 + 0.5j)
        ref_blocks.append(mat.array_interface()[h].copy() + (1.0 - 2.0j) * other.array_interface()[h])

    mat.axpy(1.0 - 2.0j, other)
    for h, ref in enumerate(ref_blocks):
        np.testing.assert_allclose(mat.to_array()[h], ref, atol=1e-14)


# ---------------------------------------------------------------------------
#  vector_dot
# ---------------------------------------------------------------------------


def _hermitian_dot(A, B):
    """Re(Tr(A^H B)), matching ComplexMatrix.vector_dot."""
    return np.vdot(A, B).real


def test_complexmatrix_vector_dot():
    """vector_dot returns Re(Tr(self^H other)) over diagonal tiles."""
    rng = np.random.default_rng(7)
    n = 5
    A_ref = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    B_ref = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))

    A = psi4.core.ComplexMatrix.from_array(A_ref, name="A")
    B = psi4.core.ComplexMatrix.from_array(B_ref, name="B")

    result = A.vector_dot(B)
    expected = _hermitian_dot(A_ref, B_ref)
    assert result == pytest.approx(expected, abs=1e-12)


def test_complexmatrix_vector_dot_multi_block():
    """vector_dot sums over multiple diagonal tiles."""
    sizes = [2, 2]
    A = psi4.core.ComplexMatrix("A", sizes)
    B = psi4.core.ComplexMatrix("B", sizes)

    rng = np.random.default_rng(13)
    expected = 0.0
    for h, n in enumerate(sizes):
        Ablk = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        Bblk = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        A.array_interface()[h][:] = Ablk
        B.array_interface()[h][:] = Bblk
        expected += _hermitian_dot(Ablk, Bblk)

    result = A.vector_dot(B)
    assert result == pytest.approx(expected, abs=1e-12)


# ---------------------------------------------------------------------------
#  product_trace
# ---------------------------------------------------------------------------


def test_complexmatrix_product_trace():
    """product_trace returns np.einsum('ij,ji->', self, other)."""
    rng = np.random.default_rng(99)
    n = 4
    A_ref = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    B_ref = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))

    A = psi4.core.ComplexMatrix.from_array(A_ref, name="A")
    B = psi4.core.ComplexMatrix.from_array(B_ref, name="B")

    result = A.product_trace(B)
    expected = np.einsum("ij,ji->", A_ref, B_ref)
    np.testing.assert_allclose(result, expected, atol=1e-12)


def test_complexmatrix_product_trace_multi_block():
    """product_trace with multi-block tensors."""
    sizes = [3, 2]
    A = psi4.core.ComplexMatrix("A", sizes)
    B = psi4.core.ComplexMatrix("B", sizes)

    rng = np.random.default_rng(101)
    expected = 0.0 + 0.0j
    for h, n in enumerate(sizes):
        Ablk = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        Bblk = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        A.array_interface()[h][:] = Ablk
        B.array_interface()[h][:] = Bblk
        expected += np.einsum("ij,ji->", Ablk, Bblk)

    result = A.product_trace(B)
    np.testing.assert_allclose(result, expected, atol=1e-12)


# ---------------------------------------------------------------------------
#  save / load
# ---------------------------------------------------------------------------


def test_complexmatrix_save_load(tmp_path):
    """save/load round-trip via PSIO."""
    import os
    old_cwd = os.getcwd()
    os.chdir(tmp_path)

    try:
        rng = np.random.default_rng(55)
        sizes = [3, 2]
        mat = psi4.core.ComplexMatrix("saveload", sizes)
        refs = []
        for h, n in enumerate(sizes):
            blk = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            mat.array_interface()[h][:] = blk
            refs.append(blk)

        psio = psi4.core.IO.shared_object()
        fileno = 150
        psio.open(fileno, 0)  # PSIO_OPEN_NEW

        mat.save(psio, fileno)

        # Load into a fresh ComplexMatrix with the same grid
        loaded = psi4.core.ComplexMatrix("saveload", sizes)
        loaded.load(psio, fileno)

        psio.close(fileno, 0)  # delete

        for h, ref in enumerate(refs):
            np.testing.assert_allclose(loaded.to_array()[h], ref, atol=1e-14)
    finally:
        os.chdir(old_cwd)


# ---------------------------------------------------------------------------
#  conjT
# ---------------------------------------------------------------------------


def test_complexmatrix_conjT_single_block():
    """conjT returns the conjugate transpose of a single-block matrix."""
    rng = np.random.default_rng(17)
    n = 4
    ref = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))

    mat = psi4.core.ComplexMatrix.from_array(ref, name="A")
    ct = mat.conjT()

    np.testing.assert_allclose(ct.to_array(), ref.conj().T, atol=1e-14)
    assert ct.name == "A^H"


def test_complexmatrix_conjT_multi_block():
    """conjT works on multi-block tensors, transposing tile indices."""
    # Rectangular blocks to verify the transpose
    mat = psi4.core.ComplexMatrix("rect", [3, 2], [2, 4])
    views = mat.array_interface()
    assert [v.shape for v in views] == [(3, 2), (2, 4)]

    rng = np.random.default_rng(23)
    refs = []
    for h, v in enumerate(views):
        blk = rng.normal(size=v.shape) + 1j * rng.normal(size=v.shape)
        v[:] = blk
        refs.append(blk)

    ct = mat.conjT()
    ct_views = ct.array_interface()

    # Transposed shapes
    assert ct_views[0].shape == (2, 3)
    assert ct_views[1].shape == (4, 2)

    for h, ref in enumerate(refs):
        np.testing.assert_allclose(ct_views[h], ref.conj().T, atol=1e-14)

    assert ct.name == "rect^H"


def test_complexmatrix_conjT_roundtrip():
    """(A^H)^H == A."""
    rng = np.random.default_rng(31)
    n = 3
    ref = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))

    mat = psi4.core.ComplexMatrix.from_array(ref, name="A")
    roundtrip = mat.conjT().conjT()

    np.testing.assert_allclose(roundtrip.to_array(), ref, atol=1e-14)
    assert roundtrip.name == "A^H^H"
