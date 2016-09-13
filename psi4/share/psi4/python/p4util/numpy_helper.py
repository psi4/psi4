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
import sys
if sys.version_info < (3,0):
    from exceptions import *

# The next three functions make me angry
def translate_interface(interface):
    """
    This is extra stupid with unicode
    """
    nouni_interface = {}
    for k, v in interface.items():
        if k == 'typestr':
            nouni_interface[k.encode('ascii', 'ignore')] = v.encode('ascii', 'ignore')
        else:
            nouni_interface[k.encode('ascii', 'ignore')] = v

    return nouni_interface

class numpy_holder(object):
    """
    Blank object, stupid. Apparently you cannot create a view directly from a dictionary
    """
    def __init__(self, interface):
        self.__array_interface__ = translate_interface(interface)

def _get_raw_views(self, copy=False):
    """
    Gets simple raw view of the passed in object.
    """
    ret = []
    for x in [numpy_holder(self.array_interface(h)) for h in xrange(self.nirrep())]:

        # Yet another hack
        if not isinstance(x.__array_interface__["shape"][0], (int, long, float)):
            x.__array_interface__["shape"] = tuple(x.__array_interface__["shape"][0])

        if 0 in x.__array_interface__["shape"]:
            ret.append(np.empty(shape=x.__array_interface__["shape"]))
        else:
            ret.append(np.array(x, copy=copy))
    return ret

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
def _dimension_from_list(self, dims, name="New Dimension"):
    """
    Builds a psi4.Dimension object from a python list or tuple. If a dimension
    object is passed a copy will be returned.
    """

    if isinstance(dims, (tuple, list)):
        irreps = len(dims)
    elif isinstance(dims, psi4.Dimension):
        irreps = dims.n()
    else:
        raise ValidationError("Dimension from list: Type '%s' not understood" % type(dims))

    ret = psi4.Dimension(irreps, name)
    for i in range(irreps):
        ret[i] = dims[i]
    return ret

def _dimension_to_tuple(dim):
    """
    Converts a psi4.Dimension object to a tuple.
    """

    if isinstance(dim, (tuple, list)):
        return tuple(dim)
    
    irreps = dim.n()
    ret = []
    for i in range(irreps):
        ret.append(dim[i])
    return tuple(ret)


@classmethod
def array_to_matrix(self, arr, name="New Matrix", dim1=None, dim2=None):
    """
    Converts a numpy array or list of numpy arrays into a Psi4 Matrix (irreped if list).
    
    Parameters
    ----------
    arr : array or list of arrays
        Numpy array or list of arrays to use as the data for a new psi4.Matrix
    name : str
        Name to give the new psi4.Matrix
    dim1 : list, tuple, or psi4.Dimension (optional)
        If a single dense numpy array is given, a dimension can be supplied to
        apply irreps to this array. Note that this discards all extra information
        given in the matrix besides the diagonal blocks determined by the passed
        dimension.
    dim2 :
        Same as dim1 only if using a Psi4.Dimension object.

    Returns
    -------
    ret : psi4.Vector or psi4.Matrix
       Returns the given Psi4 object 

    Notes
    -----
    This is a generalized function to convert a NumPy array to a Psi4 object 

    Examples
    --------
    
    >>> data = np.random.rand(20)
    >>> vector = array_to_matrix(data)

    >>> irrep_data = [np.random.rand(2, 2), np.empty(shape=(0,3)), np.random.rand(4, 4)]
    >>> matrix = array_to_matrix(irrep_data)
    >>> print matrix.rowspi().to_tuple()
    >>> (2, 0, 4)
    """

    # What type is it? MRO can help.
    arr_type = self.__mro__[0]

    # Irreped case
    if isinstance(arr, (list, tuple)):
        if (dim1 is not None) or (dim2 is not None):
            raise ValidationError("Array_to_Matrix: If passed input is list of arrays dimension cannot be specified.")

        irreps = len(arr)
        if arr_type == psi4.Matrix:
            sdim1 = psi4.Dimension(irreps)
            sdim2 = psi4.Dimension(irreps)
        
            for i in range(irreps):
                d1, d2 = _find_dim(arr[i], 2)
                sdim1[i] = d1
                sdim2[i] = d2

            ret = self(name, sdim1, sdim2)

        elif arr_type == psi4.Vector:
            sdim1 = psi4.Dimension(irreps)

            for i in range(irreps):
                d1 = _find_dim(arr[i], 1)
                sdim1[i] = d1[0]

            ret = self(name, sdim1)
        else:
            raise ValidationError("Array_to_Matrix: type '%s' is not recognized." % str(arr_type))

        for view, vals in zip(ret.nph, arr):
            if 0 in view.shape: continue
            view[:] = vals

        return ret            

    # No irreps implied by list
    else:
        if arr_type == psi4.Matrix:

            # Build an irreped array back out
            if dim1 is not None:
                if dim2 is None:
                    raise ValidationError ("Array_to_Matrix: If dim1 is supplied must supply dim2 also")
    
                dim1 = psi4.Dimension.from_list(dim1) 
                dim2 = psi4.Dimension.from_list(dim2) 

                if dim1.n() != dim2.n():
                    raise ValidationError("Array_to_Matrix: Length of passed dim1 must equal length of dim2.")
            
                ret = self(name, dim1, dim2)
            
                start1 = 0
                start2 = 0
                for num, interface in enumerate(ret.nph):
                    d1 = dim1[num]
                    d2 = dim2[num]
                    if (d1 == 0) or (d2 == 0):
                        continue

                    view = np.asarray(interface)
                    view[:] = arr[start1:start1 + d1, start2:start2 + d2]
                    start1 += d1 
                    start2 += d2

                return ret
        
            # Simple case without irreps
            else:
                ret = self(name, arr.shape[0], arr.shape[1])
                ret_view = np.asarray(numpy_holder(ret.array_interface(0)))
                ret_view[:] = arr
                return ret

        elif arr_type == psi4.Vector:
            # Build an irreped array back out
            if dim1 is not None:
                if dim2 is not None:
                    raise ValidationError ("Array_to_Matrix: If dim2 should not be supplied for 1D vectors.")
    
                dim1 = psi4.Dimension.from_list(dim1) 
                ret = self(name, dim1)

                start1 = 0
                for num, interface in enumerate(ret.nph):
                    d1 = dim1[num]
                    if (d1 == 0):
                        continue

                    view = np.asarray(interface)
                    view[:] = arr[start1:start1 + d1]
                    start1 += d1 

                return ret
        
            # Simple case without irreps
            else:
                ret = self(name, arr.shape[0])
                ret.np[:] = arr
                return ret

        else:
            raise ValidationError("Array_to_Matrix: type '%s' is not recognized." % str(arr_type))


def to_array(matrix, copy=True, dense=False):
    """
    Converts a Psi4 Matrix or Vector to a numpy array. Either copies the data or simply
    consturcts a view.
    """
    if matrix.nirrep() > 1:

        # We will copy when we make a large matrix
        if dense:
            copy = False

        ret = _get_raw_views(matrix, copy=copy)

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
            dense_ret = np.zeros(shape=(ndim1))
            start = 0
            for d1, arr in zip(dim1, ret):
                if d1 == 0: continue
                dense_ret[start: start + d1] = arr
                start += d1
        else:
            dense_ret = np.zeros(shape=(ndim1, ndim2))
            start1 = 0
            start2 = 0
            for d1, d2, arr in zip(dim1, dim2, ret):
                if d1 == 0: continue

                dense_ret[start1: start1 + d1, start2: start2 + d2] = arr
                start1 += d1
                start2 += d2

        return dense_ret

    else:
        return _get_raw_views(matrix, copy=copy)[0]

def _build_view(matrix):
    """
    Builds a view of the vector or matrix
    """
    views = to_array(matrix, copy=False, dense=False)
    if matrix.nirrep() > 1:
        return tuple(views)
    else:
        return views

@property
def _np_shape(self):
    if '_np_view_data' not in self.cdict.keys():
        self.cdict['_np_view_data'] = _build_view(self)

    view_data = self.cdict['_np_view_data']
    if self.nirrep() > 1:
        return tuple(view_data for x in range(self.nirrep()))
    else:
        return view_data.shape

@property
def _np_view(self):
    """
    View without only one irrep
    """

    if '_np_view_data' not in self.cdict.keys():
        self.cdict['_np_view_data'] = _build_view(self)

    return self.cdict['_np_view_data']

@property
def _nph_view(self):
    """
    View with irreps.
    """
    if '_np_view_data' not in self.cdict.keys():
        self.cdict['_np_view_data'] = _build_view(self)

    if self.nirrep() > 1:
        return self.cdict['_np_view_data']
    else:
        return self.cdict['_np_view_data'],

@property
def _array_conversion(self):
    if self.nirrep() > 1:
        raise ValidationError("__array__interface__ can only be called on Psi4 data holders with only one irrep!")
    else:
        return self.np.__array_interface__

# Matirx attributes
psi4.Matrix.from_array = array_to_matrix
psi4.Matrix.to_array = to_array
psi4.Matrix.shape = _np_shape
psi4.Matrix.np = _np_view
psi4.Matrix.nph = _nph_view
psi4.Matrix.__array_interface__ = _array_conversion

# Vector attributes
psi4.Vector.from_array = array_to_matrix
psi4.Vector.to_array = to_array
psi4.Vector.shape = _np_shape
psi4.Vector.np = _np_view
psi4.Vector.nph = _nph_view
psi4.Vector.__array_interface__ = _array_conversion

# Dimension attributes
psi4.Dimension.from_list = _dimension_from_list
psi4.Dimension.to_tuple = _dimension_to_tuple

# CIVector attributes

@property
def _civec_view(self):
    "Returns a view of the CIVector's buffer"
    return np.asarray(self)


psi4.CIVector.np = _civec_view
