import itertools

import numpy as np
import pytest
from psi4.core import Dimension, Matrix, Slice, doublet
import psi4
from utils import compare_arrays

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]

def check_dense_mat(m, exp_r, exp_c, exp_name=None):
    assert m.rows() == exp_r
    assert m.cols() == exp_c
    if exp_name is not None:
        assert m.name == exp_name
    assert m.symmetry() == 0
    assert m.nirrep() == 1


def check_block_sparse_mat(m, exp_nirrep, exp_rdim, exp_cdim, exp_name=None, exp_sym=0):
    assert m.symmetry() == exp_sym
    assert m.nirrep() == exp_nirrep
    assert m.rowdim() == exp_rdim
    assert m.coldim() == exp_cdim
    for row_irr in range(m.nirrep()):
        col_irr = row_irr ^ m.symmetry()
        r, c = m.nph[row_irr].shape
        assert r == exp_rdim[row_irr]
        assert c == exp_cdim[col_irr]

    if exp_name is not None:
        assert m.name == exp_name


def test_constructors():
    int_row = 10
    int_col = 20
    # int row/col
    m1 = Matrix(int_row, int_col)
    check_dense_mat(m1, int_row, int_col)

    # int row/col w/ name
    m2 = Matrix("m2", int_row, int_col)
    check_dense_mat(m2, int_row, int_col, "m2")

    dim_row = Dimension([3, 2, 1, 4])
    dim_col = Dimension([4, 2, 0, 2])
    # dim row/col (default sym)
    m3 = Matrix("m3", dim_row, dim_col)
    check_block_sparse_mat(m3, 4, dim_row, dim_col, "m3")

    # dim row/col symm specified
    m4 = Matrix("m4", dim_row, dim_col, 2)
    check_block_sparse_mat(m4, 4, dim_row, dim_col, "m4", 2)


def build_random_mat(rdim, cdim, symmetry=0):
    m = Matrix("test", rdim, cdim, symmetry)
    for h in range(m.nirrep()):
        block_shape = (m.rows(h), m.cols(h ^ m.symmetry()))
        m.nph[h][:, :] = np.random.randn(*block_shape)
    return m


def generate_result(a, b, transa, transb):
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
    GA = a.symmetry()
    GB = b.symmetry()
    GC = GA ^ GB
    c_shapes = []
    players = []
    a_blocks = a.to_array()
    b_blocks = b.to_array()
    if isinstance(a_blocks, np.ndarray):
        a_blocks = [a_blocks]
    if isinstance(b_blocks, np.ndarray):
        b_blocks = [b_blocks]
    c_blocks = []
    if transa:
        c_rowdim = a.coldim()
        link_dim_a = a.rowdim()
    else:
        c_rowdim = a.rowdim()
        link_dim_a = a.coldim()
    if transb:
        c_coldim = b.rowdim()
        link_dim_b = b.coldim()
    else:
        c_coldim = b.coldim()
        link_dim_b = b.rowdim()

    def rowsym(G, h):
        "rowsym of block h for matrix that transforms a G"
        return h

    def colsym(G, h):
        "colsym of block h for a matrix that transforms as G"
        return G ^ h

    for c_blk_idx in range(a.nirrep()):
        c_shapes.append((c_rowdim[rowsym(GC, c_blk_idx)], c_coldim[colsym(GC, c_blk_idx)]))
        # require a_blk_idx st rowsym(C, c_blk_idx) == rowsym(op(A), a_blk_idx)
        if transa:
            # if op(A) = A^T, rowsym(A^t, a_blk_idx) = colsym(A, a_blk_idx)
            # rowsym(C, c_blk_idx) = c_blk_idx
            # c_blk_idx = colsym(A, a_blk_idx)
            # c_blk_idx = GA ^ a_blk_idx
            # a_blk_idx = GA ^ c_blk_idx
            # a_blk_idx = colsym(GA, c_blk_idx)
            a_blk_idx = colsym(GA, c_blk_idx)
        else:
            # if op(A), rowsym(op(A), a_blk_idx) = a_blk_idx
            # rowsym(C, c_blk_idx) = a_blk_idx
            # c_blk_idx = a_blk_idx
            a_blk_idx = c_blk_idx

        # require b_blk_idx st colsym(C, c_blk_idx) == colsym(op(B), b_blk_idx)
        if transb:
            # if op(B) = B^T, colsym(op(B), b_blk_idx) = rowsym(B, b_blk_idx)
            # colsym(C, c_blk_idx) = rowsym(B, b_blk_idx)
            # colsym(C, c_blk_idx) = b_blk_idx
            b_blk_idx = colsym(GC, c_blk_idx)
        else:
            # op(B) = B, colsym(op(B), b_blk_idx) = colsym(B, b_blk_idx)
            # colsym(C, c_blk_idx) = colsym(B, b_blk_idx)
            # colsym(C, c_blk_idx) = GB ^ b_blk_idx
            # b_blk_idx = colsym(C, c_blk_idx) ^ GB
            b_blk_idx = colsym(GC, c_blk_idx) ^ GB
        players.append((a_blocks[a_blk_idx], b_blocks[b_blk_idx]))

    for c_blk_idx in range(a.nirrep()):
        # to compute C[c_blk_idx] we take:
        # op(A[a_blk_idx]) x op(B[b_blk_idx])
        a_blk, b_blk = players[c_blk_idx]
        if transa:
            op_a_blk = a_blk.T
        else:
            op_a_blk = a_blk
        if transb:
            op_b_blk = b_blk.T
        else:
            op_b_blk = b_blk
        # we can make sure the shapes match up
        op_a_r, op_a_c = op_a_blk.shape
        op_b_r, op_b_c = op_b_blk.shape
        # assert op_a_c == op_b_r, "block matmul: col(op(A)) != row(op(B)) [A: {}x{} op:{}] [B: {}x{} op:{}]".format(a.rowspi(a)
        c_blk_r, c_blk_c = c_shapes[c_blk_idx]
        # assert op_a_r == c_blk_r, "block matmul: row(op(A)) != row(C)"
        # assert op_b_c == c_blk_c, "block matmul: col(op(B)) != col(C)"
        c_blocks.append(np.dot(op_a_blk, op_b_blk))

    return c_blocks


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
    return f"  N(G){gsz} || G(A): {Ga} || G(B): {Gb} || doublet({a_name} x {b_name}) || {sq_or_rec.upper()}"


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
    a = build_random_mat(adl, adr, Ga)
    b = build_random_mat(bdl, bdr, Gb)
    res = doublet(a, b, at, bt)
    expected = generate_result(a, b, at, bt)
    assert res.symmetry(
    ) == a.symmetry() ^ b.symmetry(), f"Symm mismatch {a.symmetry()} x {b.symmetry()} != {res.symmetry()}"
    res_blocks = res.to_array()
    if isinstance(res_blocks, np.ndarray):
        res_blocks = [res_blocks]
    block_checks = []
    for blk_idx in range(res.nirrep()):
        assert compare_arrays(expected[blk_idx], res_blocks[blk_idx], 8, f"Block[{blk_idx}]")

def test_transform():
    # transformed = transformer.T @ original @ transformer
    rdim = Dimension([3, 2])
    cdim = Dimension([4, 5])
    original = build_random_mat(rdim, rdim, 1)

    # Test case that transformed resize needed.
    transformer = build_random_mat(rdim, cdim)
    transformed = original.clone()
    transformed.transform(transformer)
    assert compare_arrays((transformer.nph[0]).T @ original.nph[0] @ transformer.nph[1], transformed.nph[0]) 
    assert compare_arrays((transformer.nph[1]).T @ original.nph[1] @ transformer.nph[0], transformed.nph[1])

    # Test case that transformed resize not needed.
    transformer = build_random_mat(rdim, rdim)
    transformed = original.clone()
    transformed.transform(transformer)
    assert compare_arrays((transformer.nph[0]).T @ original.nph[0] @ transformer.nph[1], transformed.nph[0]) 
    assert compare_arrays((transformer.nph[1]).T @ original.nph[1] @ transformer.nph[0], transformed.nph[1])

def test_back_transform():
    # transformed = transformer @ original @ transformer.T
    rdim = Dimension([3, 2])
    cdim = Dimension([4, 5])
    original = build_random_mat(cdim, cdim, 1)

    # Test case that transformed resize needed.
    transformer = build_random_mat(rdim, cdim)
    transformed = original.clone()
    transformed.back_transform(transformer)
    assert compare_arrays(transformer.nph[0] @ original.nph[0] @ (transformer.nph[1]).T, transformed.nph[0]) 
    assert compare_arrays(transformer.nph[1] @ original.nph[1] @ (transformer.nph[0]).T, transformed.nph[1])

    # Test case that transformed resize not needed.
    transformer = build_random_mat(cdim, cdim)
    transformed = original.clone()
    transformed.back_transform(transformer)
    assert compare_arrays(transformer.nph[0] @ original.nph[0] @ (transformer.nph[1]).T, transformed.nph[0]) 
    assert compare_arrays(transformer.nph[1] @ original.nph[1] @ (transformer.nph[0]).T, transformed.nph[1])

def test_get_matrix():
    # Use get_matrix to extract a block of a non-totally symmetric matrix.
    dim = Dimension([3, 2])
    mat = Matrix("Name", dim, dim, 1)
    new_dim = Dimension([2, 1])
    control_mat = Matrix("Name", new_dim, new_dim, 1)
    i = 0
    for j in range(3):
        for k in range(2):
            # Populate both non-zero irreps of the original matrix.
            mat.set(0, j, k, i)
            mat.set(1, k, j, 10 + i)
            # Populate the parts we extract of both non-zero irreps, for the control matrix.
            if i in {0, 2}:
                control_mat.set(0, j, k, i)
                control_mat.set(1, k, j, 10 + i)
            i += 1
    new_slice = Slice(Dimension([0, 0]), new_dim)
    new_mat = Matrix("Name", new_dim, new_dim, 1)
    block = mat.get_block(new_slice, new_slice)
    assert psi4.compare_matrices(control_mat, block)

def test_set_matrix():
    # Use set_matrix to zero a block of a non-totally symmetric matrix.
    dim = Dimension([3, 2])
    mat = Matrix("Name", dim, dim, 1)
    control_mat = Matrix("Name", dim, dim, 1)
    i = 0
    for j in range(3):
        for k in range(2):
            # Populate both non-zero irreps, for the original matrix.
            mat.set(0, j, k, i)
            mat.set(1, k, j, 10 + i)
            # Populate the parts we don't zero of both non-zero irreps, for the control matrix.
            if i not in {0, 2}:
                control_mat.set(0, j, k, i)
                control_mat.set(1, k, j, 10 + i)
            i += 1
    new_dim = Dimension([2, 1])
    new_slice = Slice(Dimension([0, 0]), new_dim)
    new_mat = Matrix("Name", new_dim, new_dim, 1)
    mat.set_block(new_slice, new_slice, new_mat)
    assert psi4.compare_matrices(mat, control_mat)

def test_transpose_this():
    dim1 = Dimension([3, 2])
    dim2 = Dimension([5, 7])
    mat = Matrix("Name", dim1, dim2, 1)
    with pytest.raises(RuntimeError):
        mat.transpose_this()

def test_transpose():
    dim1 = Dimension([2, 3])
    dim2 = Dimension([4, 1])
    mat = Matrix("Matrix", dim1, dim2, 1)
    transposed_mat = Matrix("Transposed Matrix", dim2, dim1, 1)
    i = 0
    for j in range(2):
        for k in range(1):
            mat.set(0, j, k, i)
            transposed_mat.set(1, k, j, i)
            i += 1
    for j in range(3):
        for k in range(4):
            mat.set(1, j, k, i)
            transposed_mat.set(0, k, j, i)
            i += 1
    assert psi4.compare_matrices(mat.transpose(), transposed_mat)
