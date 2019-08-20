"""
Tests for Tensor constructors
"""

import numpy as np
import pytest

from psi4.core import Dimension
from psi4.linalg import Matrix_, Vector_

from .utils import compare_values

pytestmark = pytest.mark.quick


def check_dense_vec(v, exp_d, exp_label=None):
    assert v.dim == exp_d
    if exp_label is not None:
        assert v.label == exp_label
    assert v.nirrep == 1


def check_block_vec(v, exp_nirrep, exp_d, exp_label=None):
    assert v.nirrep == exp_nirrep
    assert v.dimpi == exp_d
    for irr in range(v.nirrep):
        assert v[irr].shape[0] == exp_d[irr]

    if exp_label is not None:
        assert v.label == exp_label


def check_dense_mat(m, exp_r, exp_c, exp_label=None):
    assert m.rows() == exp_r
    assert m.cols() == exp_c
    if exp_label is not None:
        assert m.label == exp_label
    assert m.symmetry == 0
    assert m.nirrep == 1


def check_block_sparse_mat(m, exp_nirrep, exp_rdim, exp_cdim, exp_label=None, exp_sym=0):
    assert m.symmetry == exp_sym
    assert m.nirrep == exp_nirrep
    assert m.rowspi == exp_rdim
    assert m.colspi == exp_cdim
    for row_irr in range(m.nirrep):
        col_irr = row_irr ^ m.symmetry
        r, c = m[row_irr].shape
        assert r == exp_rdim[row_irr]
        assert c == exp_cdim[col_irr]

    if exp_label is not None:
        assert m.label == exp_label


@pytest.mark.parametrize(
    "dtype", [pytest.param(dtype, id=f"{np.dtype(dtype).name}") for dtype in [np.float32, np.float64, np.complex128]])
def test_vector_ctors(dtype):
    int_d = 10

    # unlabeled 1-irrep vector
    v1 = Vector_(dim=int_d, dtype=dtype)
    check_dense_vec(v1, int_d)

    # labeled 1-irrep vector
    v2 = Vector_(label='v2', dim=int_d, dtype=dtype)
    check_dense_vec(v2, int_d, 'v2')

    dim = Dimension([3, 2, 1, 0])
    # unlabeled 4-irrep vector, from blocks and shapes
    v3 = Vector_(dimpi=dim, dtype=dtype)
    check_block_vec(v3, dim.n(), dim)

    # labeled 4-irrep vector
    v4 = Vector_(label='v4', dimpi=dim, dtype=dtype)
    check_block_vec(v4, dim.n(), dim, 'v4')


@pytest.mark.parametrize(
    "dtype", [pytest.param(dtype, id=f"{np.dtype(dtype).name}") for dtype in [np.float32, np.float64, np.complex128]])
def test_matrix_ctors(dtype):
    int_row = 10
    int_col = 20

    # unlabeled 1-irrep matrix
    m1 = Matrix_(rows=int_row, cols=int_col, dtype=dtype)
    check_dense_mat(m1, int_row, int_col)

    # labeled 1-irrep matrix
    m2 = Matrix_(label='m2', rows=int_row, cols=int_col, dtype=dtype)
    check_dense_mat(m2, int_row, int_col, 'm2')

    dim_row = Dimension([3, 2, 1, 4])
    dim_col = Dimension([4, 2, 0, 2])
    # labeled, 4-irrep matrix
    m3 = Matrix_(label='m3', rowspi=dim_row, colspi=dim_col, dtype=dtype)
    check_block_sparse_mat(m3, dim_row.n(), dim_row, dim_col, 'm3')

    # labeled, 4-irrep, symmetry-assigned matrix
    m4 = Matrix_(label='m4', rowspi=dim_row, colspi=dim_col, symmetry=2, dtype=dtype)
    check_block_sparse_mat(m4, dim_row.n(), dim_row, dim_col, 'm4', 2)


@pytest.mark.parametrize(
    "dtype", [pytest.param(dtype, id=f"{np.dtype(dtype).name}") for dtype in [np.float32, np.float64, np.complex128]])
def test_numpy_interop(dtype):
    int_row = 10
    int_col = 20

    # unlabeled 1-irrep matrix
    m1 = Matrix_(rows=int_row, cols=int_col, dtype=dtype)
    np_m1 = np.arange(int_row * int_col, dtype=dtype).reshape(int_row, int_col)
    m1[0] = np.arange(int_row * int_col, dtype=dtype).reshape(int_row, int_col)
    assert compare_values(np_m1, m1[0], "Assign NumPy array to unlabeled, 1-irrep matrix")

    dim_row = Dimension([3, 2, 1, 4])
    dim_col = Dimension([4, 2, 0, 2])
    # labeled, 4-irrep matrix
    m3 = Matrix_(label='m3', rowspi=dim_row, colspi=dim_col, dtype=dtype)
    # Get block, then modify it
    block = m3[3]
    block[:] = np.arange(m3[3].size, dtype=dtype).reshape(*m3[3].shape)
    assert compare_values(
        np.arange(m3[3].size, dtype=dtype).reshape(*m3[3].shape), m3[3],
        "Assign NumPy array to block extracted from labeled, 4-irrep matrix")

    # Assign to block
    m3[1] = np.arange(m3[1].size, dtype=dtype).reshape(*m3[1].shape)
    assert compare_values(
        np.arange(m3[1].size, dtype=dtype).reshape(*m3[1].shape), m3[1],
        "Assign NumPy array to a block in labeled, 4-irrep matrix")
