#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""Linear algebra utilities"""

import numpy as np

from .factory import Matrix_


def make_random_tensor_2d(rdim, cdim, symmetry=0, dtype=np.float, builder="random"):
    """Generate a blocked matrix with random entries.

    Parameters
    ----------
    rdim : Dimension
        Row dimensions
    cdim : Dimension
        Column dimensions
    symmetry: int
        Overall tensor symmetry
    dtype :
        Type of the blocks
    builder : str
         How to build the blocks
    """
    m = Matrix_(label='test', rowspi=rdim, colspi=cdim, symmetry=symmetry, dtype=dtype)

    for h in range(m.nirrep):
        block_shape = (m.rows(h), m.cols(h ^ m.symmetry))
        m[h][:, :] = _make_block[builder](block_shape, dtype=dtype)
    return m


def _random_block(shape, dtype):
    """Generate a random Hermitian, positive-definite matrix.

    Parameters
    ----------
    shape :
        The matrix dimensions.
    dtype:
        Type of the block.

    Returns
    -------
    X : array of shape [n_dim, n_dim]
        The random Hermitian, positive-definite matrix.
    """
    if dtype == np.complex64 or dtype == np.complex128:
        return np.random.randn(*shape) + np.random.randn(*shape) * 1j
    else:
        return np.random.randn(*shape)


def _hpd_block(shape, dtype):
    """Generate a random Hermitian, positive-definite matrix.

    Parameters
    ----------
    shape :
        The matrix dimensions.
    dtype:
        Type of the block.

    Returns
    -------
    X : array of correct shape
        The random Hermitian, positive-definite matrix.

    Notes
    -----
    Adapted from: https://github.com/scikit-learn/scikit-learn/blob/1495f6924/sklearn/datasets/samples_generator.py#L1213
    """

    assert shape[0] == shape[1]
    A = _random_block(shape=shape, dtype=dtype)
    U, s, V = np.linalg.svd(np.dot(A.conj().T, A))
    X = np.dot(np.dot(U, 1.0 + np.diag(np.random.rand(shape[0]))), V)

    return X


_make_block = {"random": _random_block, "hpd": _hpd_block, "spd": _hpd_block}
