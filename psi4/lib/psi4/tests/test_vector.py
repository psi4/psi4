"""
Tests for Vector class.
"""

import pytest

from psi4.core import Dimension, Vector, IntVector, Slice

pytestmark = [pytest.mark.psi, pytest.mark.api, pytest.mark.quick]


def check_dense_vec(v, exp_d, exp_name=None):
    assert v.dim() == exp_d
    if exp_name is not None:
        assert v.name == exp_name
    assert v.nirrep() == 1


@pytest.mark.parametrize("tested_class", [pytest.param(i) for i in [Vector, IntVector]])
def test_constructors(tested_class):
    int_d = 10

    # unnamed 1-irrep vector
    v1 = tested_class(int_d)
    check_dense_vec(v1, int_d)

    # named 1-irrep vector
    v2 = tested_class("v2", int_d)
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

@pytest.mark.parametrize("tested_class", [pytest.param(i) for i in [Vector, IntVector]])
def test_clone(tested_class):
    dim = Dimension([1, 2, 3])
    vec = tested_class(dim)
    copy = vec.clone()
    assert copy.dimpi() == dim

@pytest.mark.parametrize("tested_class", [pytest.param(i) for i in [Vector, IntVector]])
def test_add(tested_class):
    dim = Dimension([1, 2, 3])
    vec = tested_class(dim)
    vec.add(0, 5)
    vec.add(0, 5)
    vec.add(2, 2, 7)
    assert vec.get(0) == 10
    assert vec.get(2, 2) == 7

@pytest.mark.parametrize("tested_class", [pytest.param(i) for i in [Vector, IntVector]])
def test_set(tested_class):
    dim = Dimension([1, 2, 3])
    vec = tested_class(dim)
    vec.set(0, 5)
    vec.set(0, 5) # Deliberately doing this twice.
    vec.set(2, 2, 7)
    assert vec.get(0) == 5
    assert vec.get(2, 2) == 7

def test_int_vs_float():
    dim = Dimension(1)
    foo = IntVector(dim)
    with pytest.raises(TypeError):
        foo.set(0, 0.1)
    with pytest.raises(TypeError):
        foo.add(0, 0.1)

def test_iota():
    dim = Dimension([5, 3, 2])
    iota_vec = IntVector.iota(dim)
    for h in range(iota_vec.nirrep()):
        for i in range(iota_vec.dim(h)):
            assert iota_vec.get(h, i) == i

@pytest.mark.parametrize("tested_class", [pytest.param(i) for i in [Vector, IntVector]])
def test_init(tested_class):
    dim = Dimension([1, 2, 3])
    vec = tested_class(dim)
    vec.set(0, 5)
    vec.set(0, 5) # Deliberately doing this twice.
    vec.set(2, 2, 7)
    vec.init(Dimension(5))
    assert vec.nirrep() == 5
    assert vec.dim(2) == 0

@pytest.mark.parametrize("tested_class", [pytest.param(i) for i in [Vector, IntVector]])
def test_copy(tested_class):
    dim = Dimension([1, 2, 3])
    vec = tested_class(dim)
    vec.set(0, 5)
    vec.set(0, 5) # Deliberately doing this twice.
    vec.set(2, 2, 7)
    vec2 = tested_class(Dimension(5))
    vec.copy(vec2)
    assert vec.nirrep() == 5
    assert vec.dim(2) == 0

@pytest.mark.parametrize("tested_class", [pytest.param(i) for i in [Vector, IntVector]])
def test_zero(tested_class):
    dim = Dimension([1, 2, 3])
    vec = tested_class(dim)
    vec.set(0, 5)
    vec.set(0, 5) # Deliberately doing this twice.
    vec.set(2, 2, 7)
    vec.zero()
    assert vec.get(0) == 0
    assert vec.get(2, 2) == 0

@pytest.mark.parametrize("tested_class", [pytest.param(i) for i in [Vector, IntVector]])
def test_get_block(tested_class):
    dim = Dimension([1, 2, 3])
    vec = tested_class(dim)
    for h in range(vec.nirrep()):
        for i in range(vec.dim(h)):
            vec.set(h, i, 10 * h + i)
    my_slice = Slice(Dimension([0, 1, 0]), Dimension([0, 2, 1]))
    block = vec.get_block(my_slice)
    assert block.dimpi() == Dimension([0, 1, 1])
    assert block.get(1, 0) == 11
    assert block.get(2, 0) == 20

@pytest.mark.parametrize("tested_class", [pytest.param(i) for i in [Vector, IntVector]])
def test_set_block(tested_class):
    dim = Dimension([1, 2, 3])
    vec = tested_class(dim)
    for h in range(vec.nirrep()):
        for i in range(vec.dim(h)):
            vec.set(h, i, 10 * h + i)
    my_slice = Slice(Dimension([0, 1, 0]), Dimension([0, 2, 1]))
    block = tested_class(Dimension([0, 1, 1]))
    block.set(1, 0, 5)
    block.set(2, 0, 10)
    vec.set_block(my_slice, block)
    assert vec.get(1, 1) == 5
    assert vec.get(2, 0) == 10

@pytest.mark.parametrize("tested_class", [pytest.param(i) for i in [Vector, IntVector]])
def test_bounds(tested_class):
    dim = Dimension([1, 2, 3])
    vec = tested_class(dim)
    with pytest.raises(RuntimeError):
        vec.add(1000, 0, 0)
    with pytest.raises(RuntimeError):
        vec.get(1000, 0)
    with pytest.raises(RuntimeError):
        vec.set(1000, 0, 0)

    with pytest.raises(RuntimeError):
        vec.add(0, 1000, 0)
    with pytest.raises(RuntimeError):
        vec.get(0, 1000)
    with pytest.raises(RuntimeError):
        vec.set(0, 1000, 0)

    vec.set(1, 1000)
    assert 1000 == vec.get(1)
    vec.add(1, 1000)
    assert 2000 == vec.get(1)
