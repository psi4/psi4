from itertools import repeat

import numpy as np
import pytest

from psi4.core import Dimension
from psi4.linalg import (eig, eigvals, make_tensor_2d_from_block, name_unary_test)

from .utils import compare_values

pytestmark = pytest.mark.quick


@pytest.fixture
def linalg_eig_block():
    return np.array([[0.89342434, 0.96630682, 0.83113658, 0.9014204, 0.17622395],
                     [0.01114647, 0.93096724, 0.35509599, 0.35329223, 0.65759337],
                     [0.27868701, 0.376794, 0.63310696, 0.90892131, 0.35454718],
                     [0.02962539, 0.20561053, 0.2004051, 0.83641883, 0.08335324],
                     [0.76958296, 0.23132089, 0.33539779, 0.70616527, 0.40256713]])


eig_args = []
for group_size in [1, 2, 4, 8]:
    d = Dimension(list(repeat(10, group_size)))

    for dtype in [np.float64]:
        eig_args.append((group_size, d, dtype))


@pytest.mark.parametrize(
    "dim,dtype", [pytest.param(dim, dtype, id=name_unary_test(ni, "eig", dtype)) for ni, dim, dtype in eig_args])
def test_eig(linalg_eig_block, dim, dtype):
    A = make_tensor_2d_from_block(dim, dim, linalg_eig_block)

    ls, vs = eig(A)
    np_ls = []
    np_vs = []
    for blk in A:
        _ls, _vs = np.linalg.eig(blk)
        np_ls.append(_ls)
        np_vs.append(_vs)

    for blk_idx in range(A.nirrep):
        assert compare_values(np.sort(np_ls[blk_idx]), np.sort(ls[blk_idx]), 10, f"Eigenvalues of  block[{blk_idx}]")

        S = np.einsum("ki,kj->ij", np.conj(np_vs[blk_idx]), vs[blk_idx])
        assert compare_values(np_vs[blk_idx], S.diagonal() * vs[blk_idx], 10, f"Eigenvectors of block[{blk_idx}]")


@pytest.mark.parametrize(
    "dim,dtype", [pytest.param(dim, dtype, id=name_unary_test(ni, "eigvals", dtype)) for ni, dim, dtype in eig_args])
def test_eigvals(linalg_eig_block, dim, dtype):
    A = make_tensor_2d_from_block(dim, dim, linalg_eig_block)

    ls = eigvals(A)
    np_ls = [np.linalg.eigvals(blk) for blk in A]

    for blk_idx in range(A.nirrep):
        assert compare_values(np.sort(np_ls[blk_idx]), np.sort(ls[blk_idx]), 10, f"Eigenvalues of  block[{blk_idx}]")
