"""
Tests for Vector and Matrix classes.
"""

import itertools

import numpy as np
import pytest
from psi4.core import Dimension, Matrix_D, Vector_D, doublet
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


def test_vector_ctors():
    int_d = 10

    # unlabeled 1-irrep vector
    v1 = Vector_D(int_d)
    check_dense_vec(v1, int_d)

    # labeled 1-irrep vector
    v2 = Vector_D('v2', int_d)
    check_dense_vec(v2, int_d, 'v2')

    dim = Dimension([3, 2, 1, 0])
    # unlabeled 4-irrep vector, from blocks and shapes
    v3 = Vector_D(dim)
    check_block_vec(v3, dim.n(), dim)

    # labeled 4-irrep vector
    v4 = Vector_D('v4', dim)
    check_block_vec(v4, dim.n(), dim, 'v4')


def test_matrix_ctors():
    int_row = 10
    int_col = 20

    # unlabeled 1-irrep matrix
    m1 = Matrix_D(int_row, int_col)
    check_dense_mat(m1, int_row, int_col)

    # labeled 1-irrep matrix
    m2 = Matrix_D('m2', int_row, int_col)
    check_dense_mat(m2, int_row, int_col, 'm2')

    dim_row = Dimension([3, 2, 1, 4])
    dim_col = Dimension([4, 2, 0, 2])
    # labeled, 4-irrep matrix
    m3 = Matrix_D('m3', dim_row, dim_col)
    check_block_sparse_mat(m3, dim_row.n(), dim_row, dim_col, 'm3')

    # labeled, 4-irrep, symmetry-assigned matrix
    m4 = Matrix_D('m4', dim_row, dim_col, 2)
    check_block_sparse_mat(m4, dim_row.n(), dim_row, dim_col, 'm4', 2)


def build_random_mat(rdim, cdim, symmetry=0):
    m = Matrix_D('test', rdim, cdim, symmetry)
    for h in range(m.nirrep):
        block_shape = (m.rows(h), m.cols(h ^ m.symmetry))
        m[h][:, :] = np.random.randn(*block_shape)
    return m


def generate_result(A, B, transA, transB):
    """Generate the result of a doublet operation

    This function computes op(A) x op(B) by:
           -> Loop over blocks of C (i_c):
            1. Determine which block of A (i_a) is needed
            2. Determine which block of B (i_b) is needed
            3. Compute: C[i_c] = op(A[i_a]) x op(B[i_b])

    This a bit different than the code in `Matrix::gemm` which does:
           -> Loop over Blocks of A (i_a)
            1. Determine which block of B (i_b) is needed
            2. Determine which block of C (i_c) is the target
            3. C[i_c] = op(A[i_a]) x op(B[i_b])

    I chose to work out how to do it both ways so that this test is a bit
    stronger than just saying see I translated the function into python and both give the same result.
    """
    GA = A.symmetry
    GB = B.symmetry
    GC = GA ^ GB
    C_shapes = []
    players = []
    C_blocks = []
    if transA:
        C_rowdim = A.colspi
        link_dim_A = A.rowspi
    else:
        C_rowdim = A.rowspi
        link_dim_A = A.colspi
    if transB:
        C_coldim = B.rowspi
        link_dim_B = B.colspi
    else:
        C_coldim = B.colspi
        link_dim_B = B.rowspi

    def rowsym(G, h):
        "rowsym of block h for matrix that transforms as G"
        return h

    def colsym(G, h):
        "colsym of block h for a matrix that transforms as G"
        return G ^ h

    for C_blk_idx in range(A.nirrep):
        C_shapes.append((C_rowdim[rowsym(GC, C_blk_idx)], C_coldim[colsym(GC, C_blk_idx)]))
        # require a_blk_idx st rowsym(C, c_blk_idx) == rowsym(op(A), a_blk_idx)
        if transA:
            # if op(A) = A^T, rowsym(A^t, a_blk_idx) = colsym(A, a_blk_idx)
            # rowsym(C, c_blk_idx) = c_blk_idx
            # c_blk_idx = colsym(A, a_blk_idx)
            # c_blk_idx = GA ^ a_blk_idx
            # a_blk_idx = GA ^ c_blk_idx
            # a_blk_idx = colsym(GA, c_blk_idx)
            A_blk_idx = colsym(GA, C_blk_idx)
        else:
            # if op(A), rowsym(op(A), a_blk_idx) = a_blk_idx
            # rowsym(C, c_blk_idx) = a_blk_idx
            # c_blk_idx = a_blk_idx
            A_blk_idx = C_blk_idx

        # require b_blk_idx st colsym(C, c_blk_idx) == colsym(op(B), b_blk_idx)
        if transB:
            # if op(B) = B^T, colsym(op(B), b_blk_idx) = rowsym(B, b_blk_idx)
            # colsym(C, c_blk_idx) = rowsym(B, b_blk_idx)
            # colsym(C, c_blk_idx) = b_blk_idx
            B_blk_idx = colsym(GC, C_blk_idx)
        else:
            # op(B) = B, colsym(op(B), b_blk_idx) = colsym(B, b_blk_idx)
            # colsym(C, c_blk_idx) = colsym(B, b_blk_idx)
            # colsym(C, c_blk_idx) = GB ^ b_blk_idx
            # b_blk_idx = colsym(C, c_blk_idx) ^ GB
            B_blk_idx = colsym(GC, C_blk_idx) ^ GB
        players.append((A[A_blk_idx], B[B_blk_idx]))

    for C_blk_idx in range(A.nirrep):
        # to compute C[C_blk_idx] we take:
        # op(A[A_blk_idx]) x op(B[B_blk_idx])
        A_blk, B_blk = players[C_blk_idx]
        if transA:
            op_A_blk = A_blk.T
        else:
            op_A_blk = A_blk
        if transB:
            op_B_blk = B_blk.T
        else:
            op_B_blk = B_blk
        # we can make sure the shapes match up
        op_A_r, op_A_c = op_A_blk.shape
        op_B_r, op_B_c = op_B_blk.shape
        C_blk_r, C_blk_c = C_shapes[C_blk_idx]
        C_blocks.append(np.dot(op_A_blk, op_B_blk))

    return C_blocks


def name_doublet_test(ni, Ga, Gb, at, bt, sq_or_rec):
    gsz = ni
    if at:
        a_name = "A^T"
    else:
        a_name = "A  "
    if bt:
        b_name = "B^T"
    else:
        b_name = "B  "
    return "  N(G){} || G(A): {} || G(B): {} || doublet({} x {}) || {}".format(gsz, Ga, Gb, a_name, b_name,
                                                                               sq_or_rec.upper())


dim_choices1 = [2, 3, 4, 5, 6, 7, 8, 9]
dim_choices2 = [x + 1 for x in dim_choices1]
doublet_args = []
group_size = 4
for group_size in [1, 2, 4, 8]:
    d1 = Dimension([dim_choices1[x] for x in range(group_size)])
    d2 = Dimension([dim_choices2[x] for x in range(group_size)])

    a11_set = [(d1, d1, H) for H in range(group_size)]
    b11_set = [(d1, d1, H) for H in range(group_size)]

    for aargs, bargs, at, bt in itertools.product(a11_set, b11_set, [True, False], [True, False]):
        adl, adr, Ga = aargs
        bdl, bdr, Gb = bargs
        doublet_args.append((group_size, adl, adr, Ga, bdl, bdr, Gb, at, bt, 'square'))
    a12_set = [(d1, d2, H) for H in range(group_size)]
    b12_set = [(d1, d2, H) for H in range(group_size)]
    b21_set = [(d2, d1, H) for H in range(group_size)]

    for aargs, bargs, t in itertools.product(a12_set, b21_set, [False, True]):
        doublet_args.append((group_size, *aargs, *bargs, t, t, 'rectangular'))

    for aargs, bargs, at in itertools.product(a12_set, b12_set, [False, True]):
        bt = not at
        doublet_args.append((group_size, *aargs, *bargs, at, bt, 'rectangular'))


# If I try to prebuild the mats I run out of memory very fast, so I build the params, and create the mat w/in the test
@pytest.mark.parametrize("adl,adr,Ga,bdl,bdr,Gb,at,bt", [
    pytest.param(adl, adr, Ga, bdl, bdr, Gb, at, bt, id=name_doublet_test(ni, Ga, Gb, at, bt, sqrec))
    for ni, adl, adr, Ga, bdl, bdr, Gb, at, bt, sqrec in doublet_args
])
def test_doublets(adl, adr, Ga, bdl, bdr, Gb, at, bt):
    A = build_random_mat(adl, adr, Ga)
    B = build_random_mat(bdl, bdr, Gb)
    res = doublet(A, B, at, bt)
    expected = generate_result(A, B, at, bt)
    assert res.symmetry == A.symmetry ^ B.symmetry, "Symm mismatch {} x {} != {}".format(
        A.symmetry, B.symmetry, res.symmetry)
    block_checks = []
    for blk_idx in range(res.nirrep):
        assert compare_arrays(expected[blk_idx], res[blk_idx], 8, "Block[{}]".format(blk_idx))
