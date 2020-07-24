"""
Tests for NewVector class.
"""

import numpy as np
import pytest
from psi4.core import Dimension
from psi4.core import Tensor_D_1 as NewVector
from psi4.core import Tensor_D_2 as NewMatrix
from utils import compare_arrays


def check_dense_vec(v, exp_d, exp_name=None):
    assert v.dim == exp_d
    if exp_name is not None:
        assert v.name == exp_name
    assert v.nirrep == 1

def check_block_vec(v, exp_nirrep, exp_d, exp_name=None):
    assert v.nirrep == exp_nirrep
    assert v.dimpi == exp_d
    for irr in range(v.nirrep):
        assert v.nph[irr].shape[0] == exp_d[irr]

    if exp_name is not None:
        assert v.name == exp_name

def test_constructors():
    int_d = 10
    dim = Dimension([3, 2, 1, 0])

    # unnamed 1-irrep vector
    v1 = NewVector(int_d)
    check_dense_vec(v1, int_d)

    # named 1-irrep vector
    v2 = NewVector('v2', int_d)
    check_dense_vec(v2, int_d, 'v2')

    ## unnamed 4-irrep vector, from blocks and shapes
    #v3 = NewVector('v3', dim)
    #check_block_vec(v3, dim.n(), dim, 'v3')

    ## named 4-irrep vector
    #v4 = Vector("v4", dim)
    #check_block_vec(v4, dim.n(), dim)
