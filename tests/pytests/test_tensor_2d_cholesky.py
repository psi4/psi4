import itertools

import numpy as np
import pytest

from psi4.core import Dimension
from psi4.linalg import cholesky, make_random_tensor_2d

from .utils import compare_values

pytestmark = pytest.mark.quick


def name_cholesky_test(ni, dtype):
    return f"  N(G): {ni} || G(A): 0 || cholesky(A) || {np.dtype(dtype).name}"


dim_choices = [2, 3, 4, 5, 6, 7, 8, 9]
cholesky_args = []
for group_size in [1, 2, 4, 8]:
    d = Dimension([dim_choices[x] for x in range(group_size)])

    for dtype in [np.float32, np.float64, np.complex128]:
        cholesky_args.append((group_size, d, dtype))


@pytest.mark.parametrize(
    "dim,dtype", [pytest.param(dim, dtype, id=name_cholesky_test(ni, dtype)) for ni, dim, dtype in cholesky_args])
def test_cholesky(dim, dtype):
    A = make_random_tensor_2d(dim, dim, dtype=dtype, builder="hpd")
    res = cholesky(A)

    expected = [np.linalg.cholesky(blk) for blk in A]

    if dtype == np.float32 or dtype == np.complex64:
        atol = 1.e-6
    else:
        atol = 1.e-10

    for blk_idx, blks in enumerate(zip(expected, res)):
        assert compare_values(blks[0], blks[1], f"Cholesky factor L for block[{blk_idx}]", atol=atol)
