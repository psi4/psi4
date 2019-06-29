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

from ..core.linalg import Vector_F, Vector_D, Vector_CD, Matrix_F, Matrix_D, Matrix_CD, Tensor3_F, Tensor3_D, Tensor3_CD


class ObjectFactory:
    """An implementation of the Factory method pattern."""

    __slots__ = ["_name", "_builders"]

    def __init__(self, name, builders):
        self._name = name
        self._builders = builders

    def __call__(self, dtype, **kwargs):
        builder = self._builders.get(np.dtype(dtype))
        if not builder:
            raise ValueError(f"{self._name} of type {dtype} has no builder registered")
        return builder(**kwargs)


# NOTE for posterity. Remove the trailing "_" once the xtensor-based
# implementation supplants the old one

Vector_ = ObjectFactory(name="Vector",
                        builders={
                            np.dtype(np.float32): Vector_F,
                            np.dtype(np.float64): Vector_D,
                            np.dtype(np.complex128): Vector_CD
                        })

Matrix_ = ObjectFactory(name="Matrix",
                        builders={
                            np.dtype(np.float32): Matrix_F,
                            np.dtype(np.float64): Matrix_D,
                            np.dtype(np.complex128): Matrix_CD
                        })

Tensor_ = ObjectFactory(name="Tensor",
                        builders={
                            np.dtype(np.float32): Tensor3_F,
                            np.dtype(np.float64): Tensor3_D,
                            np.dtype(np.complex128): Tensor3_CD
                        })
