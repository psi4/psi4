import math
import numpy as np

from grendel.util.units import isunit
from grendel.gmath.tensor import Tensor
from grendel.util.strings import indented


__all__ = [
    "Matrix"
]

# TODO remove dependence on np.matrix
class Matrix(Tensor):
    """
     Encapsulates a vector.  Most functionality gets passed up to Tensor, which in turn passes things up to
    ``numpy.ndarray``.  Special functionality for vectors gets implemented here.

    """


    ##################
    # Initialization #
    ##################

    def __new__(cls, *args, **kwargs):
        """
        For now, just passes up to ``Tensor`` and view casts.  When passed a ``str`` as the first argument,
        the ``numpy.matrix`` constructor is called, and the Tensor constructor is called with the result and
        any remaining arguments.

        :Examples:


        >>> Matrix([[1,2,3],[4,5,6]])
        Matrix([[ 1.,  2.,  3.],
                [ 4.,  5.,  6.]])
        >>> Matrix((1,2,3),[4,5,6])
        Matrix([[ 1.,  2.,  3.],
                [ 4.,  5.,  6.]])
        >>> Matrix(shape=(4,3),default_val=17.5)
        Matrix([[ 17.5,  17.5,  17.5],
                [ 17.5,  17.5,  17.5],
                [ 17.5,  17.5,  17.5],
                [ 17.5,  17.5,  17.5]])

        """
        # Resolve and remove any special args/kwargs here...
        # (nothing to do yet...)
        # Then move on up and view cast
        ret_val = None
        if len(args) >= 1 and isinstance(args[0], basestring):
            ret_val = np.matrix(args[0])
            newargs = args[1:] if len(args) > 1 else ()
            ret_val = Tensor(ret_val, *newargs, **kwargs)
        else:
            ret_val = Tensor(*args, **kwargs)
        ret_val = ret_val.view(cls)
        if len(ret_val.shape) == 1:
            ret_val = ret_val.reshape(1, ret_val.shape[0])
        elif len(ret_val.shape) > 2:
            raise StandardError("Matrix objects cannot have more than two dimensions.  Use the Tensor class instead")
        return ret_val


    ##############
    # Properties #
    ##############

    @property
    def nrows(self):
        """ Number of rows that `self` has
        """
        return self.shape[0]

    @property
    def ncols(self):
        """ Number of columns that `self` has
        """
        return self.shape[1]

    @property
    def I(self):
        return self.view(np.matrix).I.view(self.__class__)

    @property
    def row_iter(self):
        for i in range(self.shape[0]):
            yield self[i]
    rows = row_iter

    @property
    def col_iter(self):
        for i in range(self.shape[1]):
            yield self[:,i]
    columns = cols = col_iter

    ###################
    # Special Methods #
    ###################

    #---------------------#
    # Container Emulation #
    #---------------------#

    def __getitem__(self, arg):
        """
        For one integer, returns the row as a ``Vector`` object (or a ``Matrix`` object if the ``no_vector_cast``
        attribute is set).  For two indices or a two integer tuple, the
        ``numpy.ndarray`` behavior is used (which just returns the entry, which is a special numpy subclass of
        ``int`` or ``float``).  Note that, as with ``Tensor``, every ``Vector`` object returned is just a view
        into the data of the ``Matrix``.  In the special case where a ``slice`` object is given for the dimension
        that turns out to be 1, the result is not view cast as a ``Vector`` but remains a ``Matrix``.
        **Note that changing an entry in the returned object changes the corresponding entry in ``self``.**

        :Examples:


        >>> m = Matrix([[1,2,3],[4,5,6]])
        >>> m[0,0]
        1.0
        >>> type(m[0,0])
        <type 'numpy.float64'>
        >>> m[0]
        Vector([ 1.,  2.,  3.])


        *Note how assignments work*::
        >>> v = m[0]
        >>> v[2] = 3.14159
        >>> m
        Matrix([[ 1.     ,  2.     ,  3.14159],
                [ 4.     ,  5.     ,  6.     ]])
        >>> v
        Vector([ 1.     ,  2.     ,  3.14159])
        >>> vert = m[:,1]
        >>> vert
        Vector([ 2.,  5.])
        >>> vert[1]=17.5
        >>> m
        Matrix([[  1.     ,   2.     ,   3.14159],
                [  4.     ,  17.5    ,   6.     ]])

        *Slices*::
        >>> m[0:2,0:2]
        Matrix([[  1. ,   2. ],
                [  4. ,  17.5]])
        >>> m[0:2,0:1]
        Matrix([[ 1.],
                [ 4.]])
        >>> m[0:1,0:1]
        Matrix([[ 1.]])
        >>> m[0:1,1]
        Vector([ 2.])
        >>> m[1,0:1]
        Vector([ 4.])
        >>> m[1:2]
        Matrix([[  4. ,  17.5,   6. ]])
        >>> m[:,1:2]
        Matrix([[  2. ],
                [ 17.5]])

        """
        ret_val = super(Matrix, self).__getitem__(arg)
        if isinstance(ret_val, EinsumTensor):
            return ret_val
        shp = ret_val.shape
        if len(shp) == 0:
            return ret_val
        if len(shp) == 1:
            return ret_val.view(Vector)
        if len(shp) > 1:
            if shp[0] == 1 and (not isinstance(arg, tuple) or not isinstance(arg[1], slice)):
                return ret_val.view(Vector)
            elif shp[1] == 1 and (not isinstance(arg, tuple) or not isinstance(arg[0], slice)):
                return ret_val.T.view(Vector)
            else:
                return ret_val
        else:
            return ret_val

    #----------------------#
    # Arithmetic Operators #
    #----------------------#

    def __mul__(self, other):
        if isinstance(other, Vector):
            # Correct matrix-vector multiplication
            mself = self.view(np.matrix)
            # ravel flattens the array without copying it
            mother = other.ravel().view(np.matrix).T
            return np.matrix.__mul__(mself, mother).ravel().view(Vector)
        elif isinstance(other, LightVector):
            # TODO Optimize this?
            # Correct matrix-vector multiplication
            mself = self.view(np.matrix)
            # ravel flattens the array without copying it
            mother = other.ravel().view(np.matrix).T
            return np.matrix.__mul__(mself, mother).ravel().view(LightVector)
        elif isinstance(other, Matrix):
            # Correct matrix-matrix multiplication
            mself = self.view(np.matrix)
            mother = other.view(np.matrix)
            return np.matrix.__mul__(mself, mother).view(Matrix)
        else:
            return super(Matrix, self).__mul__(other)

    def __pow__(self, pow):
        return np.matrix.__pow__(self.view(np.matrix), pow).view(Matrix)

    #----------------------#
    # Comparison Operators #
    #----------------------#

    def __eq__(self, other):
        return Tensor.__eq__(self, other)

    def __ne__(self, other):
        return Tensor.__ne__(self, other)

    #################
    # Class Methods #
    #################

    @classmethod
    def diagonal(cls, iterable):
        if isinstance(iterable, (Vector, np.ndarray)):
            return np.diagflat(iterable.view(Matrix))
        else:
            return Matrix(np.diagflat(iterable))

    @classmethod
    def identity(cls, n):
        return cls.diagonal([1]*n)

    ###########
    # Methods #
    ###########

    def eigensystem(self, sort=False):
        """
        Wrapper to `numpy.eigh` and `numpy.eig` for symmetric and non-symmetric
        matrices, respectively.  Note that the eigenvectors are column vectors
        in the returned matrix.
        """
        if self.is_hermitian():
            evals, evecs = np.linalg.eigh(self)
        else:
            evals, evecs = np.linalg.eig(self)
            evals = np.real_if_close(evals)
            evecs = np.real_if_close(evecs)
        if sort:
            evals, evecs = zip(*sorted((val, tuple(vec)) for val, vec in zip(evals.ravel(), evecs.col_iter)))
            evecs = Matrix(evecs).T
            evals = Vector(evals)
        else:
            evals = evals.view(Vector)
            evecs = evecs.view(Matrix)
        return evals, evecs

    def sqrt_matrix(self):
        evals, evecs = self.eigensystem()
        sqrtevals = Matrix.diagonal(np.sqrt(evals))
        return sqrtevals.transformed(evecs)

    def inverse_sqrt_matrix(self):
        evals, evecs = self.eigensystem()
        inv_sqrtevals = Matrix.diagonal(1.0/np.sqrt(evals))
        return inv_sqrtevals.transformed(evecs)

    def transformed(self, other):
        return self.linearly_transformed(other) #other * self * other.T

    def back_transformed(self, other):
        return self.linearly_transformed(other.T) #other.T * self * other

    def symmetrized(self):
        return 0.5 * (self + self.T)

    #-----------------#
    # Inquiry methods #
    #-----------------#

    def is_square(self):
        return self.nrows == self.ncols

    def is_hermitian(self):
        return self.is_square() and (self.view(np.matrix).H == self).all()

    def is_symmetric(self):
        return self.is_square() and (self.T == self).all()


#####################
# Dependent Imports #
#####################

from grendel.gmath.vector import Vector, LightVector
from grendel.gmath.einsum import EinsumTensor

