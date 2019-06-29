"""
Tests for doublets using the matrix class based on xtensor.
"""

import itertools

import numpy as np
import pytest

from psi4.core import Dimension
from psi4.linalg import Matrix_, Operation, doublet

from .utils import compare_arrays

pytestmark = pytest.mark.quick


def build_random_mat(rdim, cdim, symmetry=0, dtype=np.float):
    m = Matrix_(label='test', rowspi=rdim, colspi=cdim, symmetry=symmetry, dtype=dtype)
    for h in range(m.nirrep):
        block_shape = (m.rows(h), m.cols(h ^ m.symmetry))
        m[h][:, :] = np.random.randn(*block_shape)
        if dtype == np.complex128:
            m[h][:, :] += np.random.randn(*block_shape) *1j
    return m


def generate_result(A, B, opA, opB):
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
    stronger than just saying see I translated the function into python and
    both give the same result.
    """
    GA = A.symmetry
    GB = B.symmetry
    GC = GA ^ GB
    C_shapes = []
    players = []
    C_blocks = []
    transA = (opA == Operation.transpose or opA == Operation.transpose_conj)
    C_rowdim = A.colspi if transA else A.rowspi
    transB = (opB == Operation.transpose or opB == Operation.transpose_conj)
    C_coldim = B.rowspi if transB else B.colspi

    def rowsym(G, h):
        "rowsym of block h for matrix that transforms as G"
        return h

    def colsym(G, h):
        "colsym of block h for a matrix that transforms as G"
        return G ^ h

    def do_op(blk, op):
        op_blk = blk.T if (op == Operation.transpose or op == Operation.transpose_conj) else blk
        return (np.conjugate(op_blk) if op == Operation.transpose_conj else op_blk)

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
        op_A_blk = do_op(A_blk, opA)
        op_B_blk = do_op(B_blk, opB)

        # we can make sure the shapes match up
        C_blocks.append(np.dot(op_A_blk, op_B_blk))

    return C_blocks


naming = {Operation.none: "  ", Operation.transpose: "^T", Operation.transpose_conj: "^H"}


def name_doublet_test(ni, Ga, Gb, opA, opB, sq_or_rec, dtype=np.float):
    return f"  N(G): {ni} || G(A): {Ga} || G(B): {Gb} || doublet(A{naming[opA]} x B{naming[opB]}) || {sq_or_rec.upper()} || {np.dtype(dtype).name}"


dim_choices1 = [2, 3, 4, 5, 6, 7, 8, 9]
dim_choices2 = [x + 1 for x in dim_choices1]
doublet_args = []
group_size = 4
for group_size in [1, 2, 4, 8]:
    d1 = Dimension([dim_choices1[x] for x in range(group_size)])
    d2 = Dimension([dim_choices2[x] for x in range(group_size)])

    a11_set = [(d1, d1, H) for H in range(group_size)]
    b11_set = [(d1, d1, H) for H in range(group_size)]

    # Prune doublet_args to remove occurrences of Operation.transpose_conj with dtype double
    for aargs, bargs, opA, opB, dtype in itertools.product(
            a11_set, b11_set, [Operation.transpose_conj, Operation.transpose, Operation.none],
        [Operation.transpose_conj, Operation.transpose, Operation.none], [np.float, np.complex128]):
        if (opA == Operation.transpose_conj or opB == Operation.transpose_conj) and (dtype == np.float):
            continue
        adl, adr, Ga = aargs
        bdl, bdr, Gb = bargs
        doublet_args.append((group_size, adl, adr, Ga, bdl, bdr, Gb, opA, opB, 'square', dtype))
    a12_set = [(d1, d2, H) for H in range(group_size)]
    b12_set = [(d1, d2, H) for H in range(group_size)]
    b21_set = [(d2, d1, H) for H in range(group_size)]

    for aargs, bargs, op, dtype in itertools.product(a12_set, b21_set, [Operation.none, Operation.transpose],
                                                     [np.float, np.complex128]):
        if (op == Operation.transpose_conj) and (dtype == np.float):
            continue
        doublet_args.append((group_size, *aargs, *bargs, op, op, 'rectangular', dtype))

    for aargs, bargs, opA, dtype in itertools.product(a12_set, b12_set, [Operation.none, Operation.transpose],
                                                      [np.float, np.complex128]):
        if (opA == Operation.transpose_conj) and (dtype == np.float):
            continue
        # TODO if opA is Operation.none how do I decide whether I want Operation.transpose or Operation.transpose_conj?
        opB = Operation.none if (opA == Operation.transpose
                                 or opA == Operation.transpose_conj) else Operation.transpose
        doublet_args.append((group_size, *aargs, *bargs, opA, opB, 'rectangular', dtype))


# If I try to prebuild the mats I run out of memory very fast, so I build the params, and create the mat w/in the test
@pytest.mark.parametrize("adl,adr,Ga,bdl,bdr,Gb,opA,opB,dtype", [
    pytest.param(adl, adr, Ga, bdl, bdr, Gb, opA, opB, dtype, id=name_doublet_test(ni, Ga, Gb, opA, opB, sqrec, dtype))
    for ni, adl, adr, Ga, bdl, bdr, Gb, opA, opB, sqrec, dtype in doublet_args
])
def test_doublets(adl, adr, Ga, bdl, bdr, Gb, opA, opB, dtype):
    A = build_random_mat(adl, adr, Ga, dtype)
    B = build_random_mat(bdl, bdr, Gb, dtype)
    if dtype == np.complex128:
        res = doublet(A, B, opA, opB)
    else:
        # Use the function accepting two bools as last arguments (for coverage)
        transA = (opA == Operation.transpose or opA == Operation.transpose_conj)
        transB = (opB == Operation.transpose or opB == Operation.transpose_conj)
        res = doublet(A, B, transA, transB)
    expected = generate_result(A, B, opA, opB)
    assert res.symmetry == A.symmetry ^ B.symmetry, f"Symm mismatch {A.symmetry} x {B.symmetry} != {res.symmetry}"
    for blk_idx in range(res.nirrep):
        assert compare_arrays(expected[blk_idx], res[blk_idx], 8, f"Block[{blk_idx}]")
