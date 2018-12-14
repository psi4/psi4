"""
Tests for NewVector class.
"""

import numpy as np
import pytest
from psi4.core import Dimension, NewMatrix_D, NewVector_D
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
        assert v.nph[irr].shape[0] == exp_d[irr]

    if exp_label is not None:
        assert v.label == exp_label

def test_constructors():
    int_d = 10

    # unlabeled 1-irrep vector
    v1 = NewVector_D(int_d)
    check_dense_vec(v1, int_d)

    # labeled 1-irrep vector
    v2 = NewVector_D('v2', int_d)
    check_dense_vec(v2, int_d, 'v2')

    dim = Dimension([3, 2, 1, 0])
    # unlabeled 4-irrep vector, from blocks and shapes
    v3 = NewVector_D('v3', dim)
    check_block_vec(v3, dim.n(), dim, 'v3')

    ## labeled 4-irrep vector
    #v4 = Vector("v4", dim)
    #check_block_vec(v4, dim.n(), dim)
