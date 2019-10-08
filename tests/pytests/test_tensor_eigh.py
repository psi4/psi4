import numpy as np
import pytest

from psi4.core import Dimension
from psi4.linalg import eigh, eigvalsh, make_random_tensor_2d, name_unary_test

from .utils import compare_values

pytestmark = pytest.mark.quick

dim_choices = [1, 2, 3, 4, 5, 6, 7, 8, 9]
eig_args = []
for group_size in [1, 2, 4, 8]:
    d = Dimension([dim_choices[x] for x in range(group_size)])

    for dtype in [np.float64, np.complex128]:
        eig_args.append((group_size, d, dtype))


@pytest.mark.parametrize(
    "dim,dtype", [pytest.param(dim, dtype, id=name_unary_test(ni, "eigh", dtype)) for ni, dim, dtype in eig_args])
def test_eigh(dim, dtype):
    A = make_random_tensor_2d(dim, dim, dtype=dtype, builder="hpd")

    ls, vs = eigh(A)
    np_ls = []
    np_vs = []
    for blk in A:
        _ls, _vs = np.linalg.eigh(blk)
        np_ls.append(_ls)
        np_vs.append(_vs)

    for blk_idx in range(A.nirrep):
        assert compare_values(np.sort(np_ls[blk_idx]), np.sort(ls[blk_idx]), 10, f"Eigenvalues of  block[{blk_idx}]")

        # Take (U^T * V) product of eigenvectors, where U is the matrix
        # computed by NumPy and V the one computed by xtensros this has +/- 1
        # on the diagonal, based on the relative phases.
        # We multiply  our eigenvectors by the diagonal, to get the same phase.
        S = np.einsum("ki,kj->ij", np.conj(np_vs[blk_idx]), vs[blk_idx])
        assert compare_values(np_vs[blk_idx], S.diagonal() * vs[blk_idx], 10, f"Eigenvectors of block[{blk_idx}]")


@pytest.mark.parametrize(
    "dim,dtype", [pytest.param(dim, dtype, id=name_unary_test(ni, "eigvalsh", dtype)) for ni, dim, dtype in eig_args])
def test_eigvalsh(dim, dtype):
    A = make_random_tensor_2d(dim, dim, dtype=dtype, builder="hpd")

    ls = eigvalsh(A)
    np_ls = [np.linalg.eigvalsh(blk) for blk in A]

    for blk_idx in range(A.nirrep):
        assert compare_values(np.sort(np_ls[blk_idx]), np.sort(ls[blk_idx]), 10, f"Eigenvalues of  block[{blk_idx}]")
