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

import sys

import numpy as np

from psi4 import core
from .exceptions import ValidationError

### Matrix and Vector properties


def _get_raw_views(self, copy=False):
    """
    Gets simple raw view of the passed in object.
    """
    if copy:
        return tuple([np.array(x) for x in self.array_interface()])
    else:
        return tuple(self.array_interface())


def _find_dim(arr, ndim):
    """
    Helper function to help deal with zero or sized arrays
    """

    # Zero arrays
    if (arr is None) or (arr is False):
        return [0] * ndim

    # Make sure this is a numpy array like thing
    if not hasattr(arr, 'shape'):
        raise ValidationError("Expected numpy array, found object of type '%s'" % type(arr))

    if len(arr.shape) == ndim:
        return [arr.shape[x] for x in range(ndim)]
    else:
        raise ValidationError("Input array does not have a valid shape.")


def array_to_matrix(self, arr, name="New Matrix", dim1=None, dim2=None):
    """
    Converts a numpy array or list of numpy arrays into a Psi4 Matrix (irreped if list).

    Parameters
    ----------
    arr : array or list of arrays
        Numpy array or list of arrays to use as the data for a new core.Matrix
    name : str
        Name to give the new core.Matrix
    dim1 : list, tuple, or core.Dimension (optional)
        If a single dense numpy array is given, a dimension can be supplied to
        apply irreps to this array. Note that this discards all extra information
        given in the matrix besides the diagonal blocks determined by the passed
        dimension.
    dim2 :
        Same as dim1 only if using a psi4.core.Dimension object.

    Returns
    -------
    matrix : :py:class:`~psi4.core.Matrix` or :py:class:`~psi4.core.Vector`
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
    (2, 0, 4)
    """

    # What type is it? MRO can help.
    arr_type = self.__mro__[0]

    # Irreped case
    if isinstance(arr, (list, tuple)):
        if (dim1 is not None) or (dim2 is not None):
            raise ValidationError("Array_to_Matrix: If passed input is list of arrays dimension cannot be specified.")

        irreps = len(arr)
        if arr_type == core.Matrix:
            sdim1 = core.Dimension(irreps)
            sdim2 = core.Dimension(irreps)

            for i in range(irreps):
                d1, d2 = _find_dim(arr[i], 2)
                sdim1[i] = d1
                sdim2[i] = d2

            ret = self(name, sdim1, sdim2)

        elif arr_type == core.Vector:
            sdim1 = core.Dimension(irreps)

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
        if arr_type == core.Matrix:

            # Build an irreped array back out
            if dim1 is not None:
                if dim2 is None:
                    raise ValidationError("Array_to_Matrix: If dim1 is supplied must supply dim2 also")

                dim1 = core.Dimension.from_list(dim1)
                dim2 = core.Dimension.from_list(dim2)

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
                ret.np[:] = arr
                return ret

        elif arr_type == core.Vector:
            # Build an irreped array back out
            if dim1 is not None:
                if dim2 is not None:
                    raise ValidationError("Array_to_Matrix: If dim2 should not be supplied for 1D vectors.")

                dim1 = core.Dimension.from_list(dim1)
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


def _to_array(matrix, copy=True, dense=False):
    """
    Converts a Psi4 Matrix or Vector to a numpy array. Either copies the data or simply
    constructs a view.

    Parameters
    ----------
    matrix : :py:class:`~psi4.core.Matrix` or :py:class:`~psi4.core.Vector`
        Pointers to which Psi4 core class should be used in the construction.
    copy : bool, optional
        Copy the data if `True`, return a view otherwise
    dense : bool, optional
        Converts irreped Psi4 objects to diagonally blocked dense arrays if `True`. Returns a list of arrays otherwise.

    Returns
    -------
    array : ndarray or list of ndarray
       Returns either a list of np.array's or the base array depending on options.

    Notes
    -----
    This is a generalized function to convert a Psi4 object to a NumPy array

    Examples
    --------

    >>> data = psi4.Matrix(3, 3)
    >>> data.to_array()
    [[ 0.  0.  0.]
     [ 0.  0.  0.]
     [ 0.  0.  0.]]
    """
    if matrix.nirrep() > 1:

        # We will copy when we make a large matrix
        if dense:
            copy = False

        matrix_views = _get_raw_views(matrix, copy=copy)

        # Return the list of arrays
        if dense is False:
            return matrix_views

        # Build the dense matrix
        if isinstance(matrix, core.Vector):
            ret_type = '1D'
        elif isinstance(matrix, core.Matrix):
            ret_type = '2D'
        else:
            raise ValidationError("Array_to_Matrix: type '%s' is not recognized." % type(matrix))

        dim1 = []
        dim2 = []
        for h in matrix_views:
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
            for d1, arr in zip(dim1, matrix_views):
                if d1 == 0: continue
                dense_ret[start:start + d1] = arr
                start += d1
        else:
            dense_ret = np.zeros(shape=(ndim1, ndim2))
            start1 = 0
            start2 = 0
            for d1, d2, arr in zip(dim1, dim2, matrix_views):
                if (d1 == 0) or (d2 == 0): continue

                dense_ret[start1:start1 + d1, start2:start2 + d2] = arr
                start1 += d1
                start2 += d2

        return dense_ret

    else:
        return _get_raw_views(matrix, copy=copy)[0]


@property
def _np_shape(self):
    """
    Shape of the Psi4 data object
    """
    view_data = _get_raw_views(self)
    if self.nirrep() > 1:
        return tuple(view_data[x].shape for x in range(self.nirrep()))
    else:
        return view_data[0].shape


@property
def _np_view(self):
    """
    View without only one irrep
    """
    if self.nirrep() > 1:
        raise ValidationError("Attempted to call .np on a Psi4 data object with multiple irreps."
                              "Please use .nph for objects with irreps.")
    return _get_raw_views(self)[0]


@property
def _nph_view(self):
    """
    View with irreps.
    """
    return _get_raw_views(self)


@property
def _array_conversion(self):
    """
    Provides the array interface to simply classes so that np.array(core.Matrix(5, 5)) works flawlessly.
    """
    if self.nirrep() > 1:
        raise ValidationError("__array__interface__ can only be called on Psi4 data object with only one irrep!")
    else:
        return self.np.__array_interface__


def _np_write(self, filename=None, prefix=""):
    """
    Writes the irreped matrix to a NumPy zipped file.

    Can return the packed data for saving many matrices into the same file.
    """

    ret = {}
    ret[prefix + "Irreps"] = self.nirrep()
    ret[prefix + "Name"] = self.name
    for h, v in enumerate(self.nph):
        # If returning arrays to user, we want to return copies (snapshot), not
        # views of the core.Matrix's memory.
        if filename is None and not v.flags['OWNDATA']:
            v = np.copy(v)
        ret[prefix + "IrrepData" + str(h)] = v

    if isinstance(self, core.Matrix):
        ret[prefix + "Dim1"] = self.rowdim().to_tuple()
        ret[prefix + "Dim2"] = self.coldim().to_tuple()
    if isinstance(self, core.Vector):
        ret[prefix + "Dim"] = [self.dim(x) for x in range(self.nirrep())]

    if filename is None:
        return ret

    np.savez(filename, **ret)


def _np_read(self, filename, prefix=""):
    """
    Reads the data from a NumPy compress file.
    """

    if isinstance(filename, np.lib.npyio.NpzFile):
        data = filename
    elif isinstance(filename, str):
        if not filename.endswith('.npz'):
            filename = filename + '.npz'

        data = np.load(filename)
    else:
        raise Exception("Filename not understood: %s" % filename)

    ret_data = []

    if ((prefix + "Irreps") not in data.keys()) or ((prefix + "Name") not in data.keys()):
        raise ValidationError("File %s does not appear to be a numpyz save" % filename)

    for h in range(data[prefix + "Irreps"]):
        ret_data.append(data[prefix + "IrrepData" + str(h)])

    arr_type = self.__mro__[0]
    if arr_type == core.Matrix:
        dim1 = core.Dimension.from_list(data[prefix + "Dim1"])
        dim2 = core.Dimension.from_list(data[prefix + "Dim2"])
        ret = self(str(data[prefix + "Name"]), dim1, dim2)
    elif arr_type == core.Vector:
        dim1 = core.Dimension.from_list(data[prefix + "Dim"])
        ret = self(str(data[prefix + "Name"]), dim1)

    for h in range(data[prefix + "Irreps"]):
        ret.nph[h][:] = ret_data[h]

    return ret


def _to_serial(data):
    """
    Converts an object with a .nph accessor to a serialized dictionary
    """

    json_data = {}
    json_data["shape"] = []
    json_data["data"] = []

    for view in data.nph:
        json_data["shape"].append(view.shape)
        json_data["data"].append(view.tostring())

    if len(json_data["shape"][0]) == 1:
        json_data["type"] = "vector"
    elif len(json_data["shape"][0]) == 2:
        json_data["type"] = "matrix"
    else:
        raise ValidationError("_to_json is only used for vector and matrix objects.")

    return json_data


def _from_serial(self, json_data):
    """
    Converts serialized data to the correct Psi4 data type
    """

    if json_data["type"] == "vector":
        dim1 = core.Dimension.from_list([x[0] for x in json_data["shape"]])
        ret = self("Vector from JSON", dim1)
    elif json_data["type"] == "matrix":
        dim1 = core.Dimension.from_list([x[0] for x in json_data["shape"]])
        dim2 = core.Dimension.from_list([x[1] for x in json_data["shape"]])
        ret = self("Matrix from JSON", dim1, dim2)
    else:
        raise ValidationError("_from_json did not recognize type option of %s." % str(json_data["type"]))

    for n in range(len(ret.nph)):
        ret.nph[n].flat[:] = np.frombuffer(json_data["data"][n], dtype=np.double)

    return ret


def _chain_dot(*args, **kwargs):
    """
    Chains dot products together from a series of Psi4 Matrix classes.

    By default there is no transposes, an optional vector of booleans can be passed in.
    """

    trans = kwargs.pop("trans", None)
    if trans is None:
        trans = [False for x in range(len(args))]
    else:
        if len(trans) != len(args):
            raise ValidationError(
                "Chain dot: The length of the transpose arguements is not equal to the length of args.")

    # Setup chain
    ret = args[0]
    if trans[0]:
        ret = ret.transpose()

    # Run through
    for n, mat in enumerate(args[1:]):
        ret = core.doublet(ret, mat, False, trans[n + 1])

    return ret


def _irrep_access(self, *args, **kwargs):
    """
    Warns user when iterating/accessing an irreped object.
    """
    raise ValidationError("Attempted to access by index/iteration a Psi4 data object that supports multiple"
                          "irreps. Please use .np or .nph explicitly.")


# Matrix attributes
core.Matrix.from_array = classmethod(array_to_matrix)
core.Matrix.from_list = classmethod(lambda self, x: array_to_matrix(self, np.array(x)))
core.Matrix.to_array = _to_array
core.Matrix.shape = _np_shape
core.Matrix.np = _np_view
core.Matrix.nph = _nph_view
core.Matrix.__array_interface__ = _array_conversion
core.Matrix.np_write = _np_write
core.Matrix.np_read = classmethod(_np_read)
core.Matrix.to_serial = _to_serial
core.Matrix.from_serial = classmethod(_from_serial)
core.Matrix.chain_dot = _chain_dot
core.Matrix.__iter__ = _irrep_access
core.Matrix.__getitem__ = _irrep_access

# Vector attributes
core.Vector.from_array = classmethod(array_to_matrix)
core.Vector.from_list = classmethod(lambda self, x: array_to_matrix(self, np.array(x)))
core.Vector.to_array = _to_array
core.Vector.shape = _np_shape
core.Vector.np = _np_view
core.Vector.nph = _nph_view
core.Vector.__array_interface__ = _array_conversion
core.Vector.np_write = _np_write
core.Vector.np_read = classmethod(_np_read)
core.Vector.to_serial = _to_serial
core.Vector.from_serial = classmethod(_from_serial)
core.Vector.__iter__ = _irrep_access
core.Vector.__getitem__ = _irrep_access

### CIVector properties


@property
def _civec_view(self):
    """
    Returns a view of the CIVector's buffer
    """
    return np.asarray(self)


core.CIVector.np = _civec_view

### Dimension properties


@classmethod
def _dimension_from_list(self, dims, name="New Dimension"):
    """
    Builds a core.Dimension object from a python list or tuple. If a dimension
    object is passed a copy will be returned.
    """

    if isinstance(dims, (tuple, list, np.ndarray)):
        irreps = len(dims)
    elif isinstance(dims, core.Dimension):
        irreps = dims.n()
    else:
        raise ValidationError("Dimension from list: Type '%s' not understood" % type(dims))

    ret = core.Dimension(irreps, name)
    for i in range(irreps):
        ret[i] = dims[i]
    return ret


def _dimension_to_tuple(dim):
    """
    Converts a core.Dimension object to a tuple.
    """

    if isinstance(dim, (tuple, list)):
        return tuple(dim)

    irreps = dim.n()
    ret = []
    for i in range(irreps):
        ret.append(dim[i])
    return tuple(ret)


def _dimension_iter(dim):
    """
    Provides an iterator class for the Dimension object.

    Allows:
        dim = psi4.core.Dimension(...)
        list(dim)
    """

    for i in range(dim.n()):
        yield dim[i]


# Dimension attributes
core.Dimension.from_list = _dimension_from_list
core.Dimension.to_tuple = _dimension_to_tuple
core.Dimension.__iter__ = _dimension_iter


# General functions for NumPy array manipulation
def block_diagonal_array(*args):
    """
    Convert square NumPy array to a single block diagonal array.
    Mimic of SciPy's block_diag.
    """

    # Validate the input matrices.
    dim = 0
    for matrix in args:
        try:
            shape = matrix.shape
            dim += shape[0]
        except (AttributeError, TypeError):
            raise ValidationError("Cannot construct block diagonal from non-arrays.")
        if len(shape) != 2:
            raise ValidationError("Cannot construct block diagonal from non-2D arrays.")
        if shape[0] != shape[1]:
            raise ValidationError("Cannot construct block diagonal from non-square arrays.")

    # If this is too slow, try a sparse matrix?
    block_diag = np.zeros((dim, dim))
    start = 0
    for matrix in args:
        next_block = slice(start, start + matrix.shape[0])
        block_diag[next_block, next_block] = matrix
        start += matrix.shape[0]

    return block_diag
