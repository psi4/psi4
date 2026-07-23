import pytest
import numpy as np

import psi4
from psi4.core import Dimension, Matrix
from psi4.driver.p4util.exceptions import ValidationError

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick, pytest.mark.cghf]


def test_complexmatrix_constructor_zeroed():
    """ComplexMatrix(name, block_sizes) builds a zeroed multi-block tensor."""
    mat = psi4.core.ComplexMatrix("zeros", [2, 3, 1])
    assert mat.num_blocks() == 3

    views = mat.array_interface()
    assert len(views) == 3
    assert [v.shape for v in views] == [(2, 2), (3, 3), (1, 1)]
    assert all(v.dtype == np.complex128 for v in views)

    blocks = mat.to_array()
    assert isinstance(blocks, list)
    for blk in blocks:
        np.testing.assert_allclose(blk, np.zeros_like(blk))


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
    assert mat.num_blocks() == 1
    out = mat.to_array()
    assert isinstance(out, np.ndarray)
    np.testing.assert_allclose(out, ref)

    # Real input is cast to complex128
    real = np.eye(2)
    cast = psi4.core.ComplexMatrix.from_array(real)
    assert cast.to_array().dtype == np.complex128
    np.testing.assert_allclose(cast.to_array(), real)

    with pytest.raises(ValidationError, match="square"):
        psi4.core.ComplexMatrix.from_array(np.ones((2, 3), dtype=np.complex128))

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


def test_share_matrix_to_einsums_c1():
    """share_matrix_to_einsums aliases a C1 Matrix as a BlockTensorView (no copy)."""
    ref = np.arange(9, dtype=float).reshape(3, 3)
    M = Matrix.from_array(ref)
    bt = psi4.core.share_matrix_to_einsums(M)

    assert isinstance(bt, psi4.core.BlockTensorView)
    assert bt.num_blocks() == 1
    assert not hasattr(bt, "from_array")
    np.testing.assert_allclose(bt.to_array(), ref)

    # View -> Matrix
    bt.to_array(copy=False)[1, 1] = 42.0
    assert M.np[1, 1] == pytest.approx(42.0)

    # Matrix -> View
    M.np[0, 0] = -5.0
    assert bt.to_array()[0, 0] == pytest.approx(-5.0)


def test_share_matrix_to_einsums_multi_irrep():
    """Multi-irrep Matrix (including empty blocks) shares correctly with BlockTensorView."""
    rdim = Dimension([2, 1, 3, 0])
    cdim = Dimension([2, 1, 3, 0])
    M = Matrix("sym", rdim, cdim)
    rng = np.random.default_rng(0)
    for h in range(M.nirrep()):
        if M.rows(h) and M.cols(h):
            M.nph[h][:] = rng.normal(size=(M.rows(h), M.cols(h)))

    bt = psi4.core.share_matrix_to_einsums(M)
    assert bt.num_blocks() == M.nirrep()

    mat_blocks = M.to_array()
    bt_blocks = bt.to_array()
    assert isinstance(bt_blocks, list)
    assert len(bt_blocks) == 4
    for h, (a, b) in enumerate(zip(mat_blocks, bt_blocks)):
        assert a.shape == b.shape
        np.testing.assert_allclose(a, b)

    # Empty irrep is a (0, 0) view
    assert bt.array_interface()[3].shape == (0, 0)

    # Write-through on a non-empty block
    bt.to_array(copy=False)[2][0, 0] = 99.0
    assert M.nph[2][0, 0] == pytest.approx(99.0)


def test_blocktensorview_keeps_matrix_alive():
    """Returned BlockTensorView keeps the source Matrix alive after the Python ref drops."""
    bt = psi4.core.share_matrix_to_einsums(Matrix.from_array(np.eye(2) * 3.0))
    np.testing.assert_allclose(bt.to_array(), np.eye(2) * 3.0)
    bt.to_array(copy=False)[0, 1] = 1.5
    assert bt.to_array()[0, 1] == pytest.approx(1.5)
