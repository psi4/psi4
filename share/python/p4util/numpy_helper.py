#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import numpy as np
import psi4
from exceptions import *


def _find_dim(arr, ndim):
    """
    Helper function to help deal with zero or sized arrays
    """

    # Zero arrays
    if (arr is None) or (arr is False):
        return [0] * ndim

    # Make sure this is a numpy array like thing
    try:
        arr.shape
    except:
        raise ValidationError("Expected numpy array, found object of type '%s'", type(arr))

    if len(arr.shape) == ndim:
        return [arr.shape[x] for x in range(ndim)]
    else:
        raise ValidationError("Input array does not have a valid shape.")

@classmethod
def array_to_matrix(self, arr, name="New Matrix"):
    """
    Converts a numpy array or list of numpy arrays into a Psi4 Matrix (irreped if list)
    """

    # What type is it? MRO can help.
    arr_type = self.__mro__[0]
    non_arr = np.array(None)

    # Irreped case
    if isinstance(arr, (list, tuple)):
        irreps = len(arr)
        if arr_type == psi4.Matrix:
            dim1 = psi4.Dimension(irreps)
            dim2 = psi4.Dimension(irreps)
        
            for i in range(irreps):
                d1, d2 = _find_dim(arr[i], 2)
                dim1[i] = d1
                dim2[i] = d2

            ret = self(name, dim1, dim2)

        elif arr_type == psi4.Vector:
            dim1 = psi4.Dimension(irreps)

            for i in range(irreps):
                d1 = _find_dim(arr[i], 1)
                dim1[i] = d1[0]

            ret = self(name, dim1)
        else:
            raise ValidationError("Array_to_Matrix: type '%s' is not recognized." % str(arr_type))

        views = []
        for interface in ret.array_interfaces():
            if 0 in interface.__array_interface__["shape"]:
                views.append(None)
            else:
                views.append(np.asarray(interface))

        for i in range(irreps):
            if views[i] is None: continue
            views[i][:] = arr[i]

        return ret            

    # No irreps
    else:
        if arr_type == psi4.Matrix:
            ret = self(name, arr.shape[0], arr.shape[1])
        elif arr_type == psi4.Vector:
            ret = self(name, arr.shape[0])
        else:
            raise ValidationError("Array_to_Matrix: type '%s' is not recognized." % str(arr_type))

        ret_view = np.asarray(ret)
        ret_view[:] = arr
        return ret

def to_array(matrix, copy=True, dense=False):
    """
    Converts a Psi4 Matrix or Vector to a numpy array. Either copies the data or simply
    consturcts a view.
    """
    if matrix.nirrep() > 1:

        # We will copy when we make a large matrix
        if dense:
            copy = False

        ret = []
        for h in matrix.array_interfaces():
            if 0 in h.__array_interface__["shape"]:
                ret.append(np.empty(shape = h.__array_interface__["shape"]))
            else:
                ret.append(np.array(h, copy=copy))

        # Return the list of arrays
        if dense is False:
            return ret
        
        # Build the dense matrix
        if isinstance(matrix, psi4.Vector):
            ret_type = '1D'
        elif isinstance(matrix, psi4.Matrix):
            ret_type = '2D'
        else:
            raise ValidationError("Array_to_Matrix: type '%s' is not recognized." % type(matrix))

        dim1 = []
        dim2 = []
        for h in ret:
            # Ignore zero dim irreps
            if 0 in h.shape:
                dim1.append(0)
                dim2.append(0)
            else:
                dim1.append(h.shape[0])
                if ret_type == '2D':
                    dim2.append(h.shape[1])

        ndim1 = np.sum(dim1)
        ndim2 = np.sum(dim2)
        if ret_type == '1D':
            dense_ret = np.empty(shape=(ndim1))
            start = 0
            for d1, arr in zip(dim1, ret):
                if d1 == 0: continue
                dense_ret[start: start + d1] = arr
                start += d1
        else:
            dense_ret = np.empty(shape=(ndim1, ndim2))
            start1 = 0
            start2 = 0
            for d1, d2, arr in zip(dim1, dim2, ret):
                if d1 == 0: continue
                dense_ret[start1: start1 + d1, start2: start2 + d2] = arr
                start1 += d1
                start1 += d2

        return dense_ret

    else:
        if 0 in matrix.__array_interface__["shape"]:
            return np.empty(shape = matrix.__array_interface__["shape"])
        else:
            return np.array(matrix, copy=copy)


# Matirx attributes
psi4.Matrix.from_array = array_to_matrix
psi4.Matrix.to_array = to_array

# Vector attributes
psi4.Vector.from_array = array_to_matrix
psi4.Vector.to_array = to_array
