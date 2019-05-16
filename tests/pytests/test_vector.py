"""
Tests for Vector class.
"""

import pytest

from psi4.core import Dimension, Vector

pytestmark = pytest.mark.quick


def check_dense_vec(v, exp_d, exp_name=None):
    assert v.dim() == exp_d
    if exp_name is not None:
        assert v.name == exp_name
    assert v.nirrep() == 1


def test_constructors():
    int_d = 10

    # unnamed 1-irrep vector
    v1 = Vector(int_d)
    check_dense_vec(v1, int_d)

    # named 1-irrep vector
    v2 = Vector("v2", int_d)
    check_dense_vec(v2, int_d, "v2")


def check_block_vec(v, exp_nirrep, exp_d, exp_name=None):
    assert v.nirrep() == exp_nirrep
    assert v.dimpi() == exp_d
    for irr in range(v.nirrep()):
        assert v.nph[irr].shape[0] == exp_d[irr]

    if exp_name is not None:
        assert v.name == exp_name


dim_choices1 = [2, 3, 4, 5, 6, 7, 8, 9]
test_args = []
for group_size in [1, 2, 4, 8]:
    name = "v{:d}".format(group_size)
    dim = Dimension([dim_choices1[x] for x in range(group_size)])
    test_args.append((name, dim))


@pytest.mark.parametrize("name,dim", [pytest.param(name, dim) for name, dim in test_args])
def test_constructors_w_symmetry(name, dim):
    v = Vector(name, dim)
    check_block_vec(v, dim.n(), dim, name)
