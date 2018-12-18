"""
Tests for NewVector and NewMatrix classes.
"""

import numpy as np
import pytest
from psi4.core import Dimension, NewMatrix_D, NewVector_D
from utils import compare_arrays


def check_dense_vec(v, exp_d, exp_label=None):
    assert v.dim == exp_d
    if exp_label is not None:
        assert v.label == exp_label
    assert v.nirrep == 1


def check_block_vec(v, exp_nirrep, exp_d, exp_label=None):
    assert v.nirrep == exp_nirrep
    assert v.dimpi == exp_d
    for irr in range(v.nirrep):
        assert v.nph[irr][0] == exp_d[irr]

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
        r, c = m.nph[row_irr]
        assert r == exp_rdim[row_irr]
        assert c == exp_cdim[col_irr]

    if exp_label is not None:
        assert m.label == exp_label


def test_vector_ctors():
    int_d = 10

    # unlabeled 1-irrep vector
    v1 = NewVector_D(int_d)
    check_dense_vec(v1, int_d)

    # labeled 1-irrep vector
    v2 = NewVector_D('v2', int_d)
    check_dense_vec(v2, int_d, 'v2')

    dim = Dimension([3, 2, 1, 0])
    # unlabeled 4-irrep vector, from blocks and shapes
    v3 = NewVector_D(dim)
    check_block_vec(v3, dim.n(), dim)

    # labeled 4-irrep vector
    v4 = NewVector_D('v4', dim)
    check_block_vec(v4, dim.n(), dim, 'v4')


def test_matrix_ctors():
    int_row = 10
    int_col = 20

    # unlabeled 1-irrep matrix
    m1 = NewMatrix_D(int_row, int_col)
    check_dense_mat(m1, int_row, int_col)

    # labeled 1-irrep matrix
    m2 = NewMatrix_D('m2', int_row, int_col)
    check_dense_mat(m2, int_row, int_col, 'm2')

    dim_row = Dimension([3, 2, 1, 4])
    dim_col = Dimension([4, 2, 0, 2])
    # labeled, 4-irrep matrix
    m3 = NewMatrix_D('m3', dim_row, dim_col)
    check_block_sparse_mat(m3, dim_row.n(), dim_row, dim_col, 'm3')

    # labeled, 4-irrep, symmetry-assigned matrix
    m4 = NewMatrix_D('m4', dim_row, dim_col, 2)
    check_block_sparse_mat(m4, dim_row.n(), dim_row, dim_col, 'm4', 2)
