import itertools

import numpy as np
import pytest

from psi4.core import Dimension
from psi4.linalg import Matrix_D, cholesky, make_random_tensor_2d

from .utils import compare_arrays

pytestmark = pytest.mark.quick


def generate_result(A):
    C_blocks = []
    for blk_idx in range(A.nirrep):
        C_blocks.append(np.linalg.cholesky(A[blk_idx]))
    return C_blocks


def name_cholesky_test(ni, dtype):
    return f"  N(G): {ni} || G(A): 0 || cholesky(A) || {np.dtype(dtype).name}"


dim_choices = [2, 3, 4, 5, 6, 7, 8, 9]
cholesky_args = []
for group_size in [1, 2, 4, 8]:
    d = Dimension([dim_choices[x] for x in range(group_size)])

    for aargs, dtype in itertools.product([(d, d) for _ in range(group_size)], [np.float, np.complex128]):
        cholesky_args.append((group_size, aargs[0], aargs[1], dtype))


@pytest.mark.parametrize(
    "adl,adr,dtype",
    [pytest.param(adl, adr, dtype, id=name_cholesky_test(ni, dtype)) for ni, adl, adr, dtype in cholesky_args])
def test_cholesky(adl, adr, dtype):
    A = make_random_tensor_2d(adl, adr, symmetry=0, dtype=dtype, builder="hpd")
    res = cholesky(A)

    expected = generate_result(A)

    for blk_idx in range(res.nirrep):
        assert compare_arrays(expected[blk_idx], res[blk_idx], 10, f"Block[{blk_idx}]")
