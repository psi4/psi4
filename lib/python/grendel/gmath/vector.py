from copy import copy
from warnings import warn
import numpy as np
import math

from grendel import sanity_checking_enabled, show_warnings
from grendel.gmath.tensor import Tensor, LightTensor
from grendel.util.aliasing import function_alias
from grendel.util.exceptions import ProgrammerNeedsMoreCoffeeError

__all__ = [
    "Vector",
    "cross",
    "magnitude", "norm"
]


#####################
# Wrappers to numpy #
#####################

def cross(v1, v2, *args, **kwargs):
    """ Wrapper to `numpy.cross`
    which viewcasts the result as a `Vector` if the first two arguments are `Vector`  objects

    :Examples:

    >>> cross( Vector(1, 2, 3), Vector(4, 5, 6) )
    Vector([-3.,  6., -3.])
    >>> # Can also be called as a method of Vector
    >>> Vector(1, 2, 3).cross( Vector(4, 5, 6) )
    Vector([-3.,  6., -3.])

    """
    if not (isinstance(v1, Tensor) or isinstance(v2, Tensor)):
        # View as a numpy.ndarray
        return np.cross(v1, v2, *args, **kwargs)
    elif not (isinstance(v1, Vector) or isinstance(v2, Vector)):
        # View as a Tensor
        return np.cross(v1, v2, *args, **kwargs).view(Tensor)
    else:
        return np.cross(v1, v2, *args, **kwargs).view(Vector)

def l_cross(v1, v2):
    """ "Light" version of `cross()`.  Always returns a `LightVector`, so
    make sure you know what you're doing before you use this.
    """
    return np.cross(v1, v2).view(LightVector)


def magnitude(v):
    """ Wrapper for `Vector.magnitude()`

    :Examples:

    >>> magnitude(Vector(3.0, 4.0))
    5.0

    """
    return v.magnitude()
norm = magnitude

class LightVector(LightTensor):

    ##################
    # Initialization #
    ##################

    def __new__(cls, iterable):
        return np.array(iterable).view(cls)

    ##################
    # Static Methods #
    ##################

    @staticmethod
    def l_sub(a, b):
        return LightVector(a.view(np.ndarray) - b.view(np.ndarray))

    @staticmethod
    def l_add(a, b):
        return LightVector(a.view(np.ndarray) + b.view(np.ndarray))

    ###########
    # Methods #
    ###########

    #-----------------------------#
    # Methods that return scalars #
    #-----------------------------#

    def magnitude(self):
        """ The magnitude of the vector.
        """
        return math.sqrt(self.dot(self))
        # Aliases
    norm = function_alias('norm', magnitude)

    def size(self):
        """ The number of entries in the vector
        """
        return self.shape[0]

    #---------------------------------------#
    # "Light" versions of methods in Vector #
    #---------------------------------------#

    def normalize(self):
        """ Performs an in-place normalization of `self`.

        :Examples:


        >>> Vector(2.0, 0.0, 0.0).normalize()
        Vector([ 1.,  0.,  0.])
        >>> Vector(0.5, 0.5, 0.5).normalize()
        Vector([ 0.57735027,  0.57735027,  0.57735027])
        >>> v = Vector(-0.5, 0.5, -0.5).normalize()
        >>> v
        Vector([-0.57735027,  0.57735027, -0.57735027])
        >>> v.magnitude()
        1.0

        """
        if self.is_zero():
            raise ZeroDivisionError("can't normalize a zero vector.  (A zero vector, "
                                    "in this case, is one with a magnitude less than {0}.  "
                                    "The vector you tried to normalize has a magnitude of "
                                    "{1}".format(self.zero_cutoff, self.magnitude()))
        norm = self.magnitude()
        for i in xrange(self.size()):
            self[i] /= norm
        return self

    def l_normalized(self):
        """ Same as `l_normalize()`, but does not modify self
        """
        ret_val = copy(self)
        ret_val.normalize()
        return ret_val.view(LightVector)
    normalized=function_alias('normalized', l_normalized)

    def l_cross(self, other):
        """ Returns cross product of self with other.  Wrapper to `numpy.cross`
        """
        return l_cross(self, other)

    #----------------------------------------------------------------------#
    # Boolean-returning methods for determing aspects of the vector `self` #
    #----------------------------------------------------------------------#

    def is_cartesian(self):
        """ Returns True if the vector is Cartesian
        (i.e. if the vector is two or three dimensional).
        """
        return self.size() == 2 or self.size() == 3

    def is_zero_vector(self, cutoff=None):
        """
        Alias for Tensor.is_zero()
        """
        return self.is_zero(cutoff)


class Vector(LightVector, Tensor):
    """
    Encapsulates a vector.  Most functionality gets passed up to Tensor, which in turn passes things up to
    ``numpy.ndarray``.  Special functionality for vectors gets implemented here.

    :Examples:


    * Multiplication *
    >>> from grendel.gmath import Matrix, Vector
    >>> Matrix([1,2],[3,4]) * Vector([5,6])
    Vector([ 17.,  39.])

    * Cartesian components *
    >>> vtwo = Vector(1., 2.)
    >>> vthree = Vector(1.11, 2.22, 3.33)
    >>> vtwo[0]
    1.0
    >>> vtwo.x
    1.0
    >>> vtwo.y == vtwo[1]
    True
    >>> vthree[2] == vthree.z
    True
    >>> # Can also set:
    >>> vtwo.x = -1
    >>> vtwo.y = -2
    >>> vtwo
    Vector([-1., -2.])
    >>> vthree.x, vthree.z, vthree.y = 1, 3, 2
    >>> vthree
    Vector([ 1.,  2.,  3.])

    """

    ####################
    # Class Attributes #
    ####################


    ##################
    # Initialization #
    ##################

    def __new__(cls, *args, **kwargs):
        """
        For now, just passes up to Tensor and view casts
        Note: vectors can only be 1-dimensional.  No matter what.  If a multi-dimensional Tensor is view-cast as
        a vector, it will be flattened in row-major order (try not to do this, though)
        """
        # Resolve and remove any special args/kwargs here...

        # Then move on up and view cast
        ret_val = None
        if len(kwargs) == 0:
            ret_val = Tensor(*args)
        else:
            ret_val = Tensor(*args, **kwargs)

        ret_val = ret_val.view(cls)

        return ret_val

    ###################
    # Special Methods #
    ###################

    def __array_finalize__(self, obj):
        """
        See http://docs.scipy.org/doc/numpy/user/basics.subclassing.html#array-finalize for explanation
        Initialize any defaults that may need to be initialized in a case where the __new__ is not called explicitly
        (but remember that this also gets called when __new__ is called explicitly)
        Remember that both view casting and new-from-template must be handled here.
        """
        # Flatten self.  No matter what.
        if len(self.shape) > 1:
            self.shape = -1
        # Call super
        super(Vector, self).__array_finalize__(obj)


    def __rmul__(self, other):
        # check for the strange behavior that causes numpy to think of Vector objects as row vectors
        if sanity_checking_enabled:
            if hasattr(other, 'shape'):
                if len(other.shape) == 2:
                    if not isinstance(other, Matrix) and show_warnings:
                        warn("matrix-like object that is not a Matrix multiplied by vector.  Multiplication will be "
                             "element-wise, not standard matrix-vector multiplication.  Please catch this error within"
                             " grendel if you know what you are doing.", MatrixMultiplicationWarning)
                    else: # pragma: no cover
                        # This should never happen.  Matrix.__mul__ should always be called and handel this.  I'm
                        # putting this here to guard against future changes that could cause this to happen.
                        raise ProgrammerNeedsMoreCoffeeError
        #--------------------------------------------------------------------------------#
        # call super
        return super(Vector, self).__rmul__(other)


    ##############
    # Properties #
    ##############

    @property
    def x(self):
        """ The x component of a cartesian vector

        :Raises:

        IndexError
            If the vector is not cartesian (i.e. if it is not two or three dimensional)

        """
        if self.is_cartesian():
            return self.__getitem__(0)
        else:
            raise IndexError("Getting the x component of a vector that is not in 2 or 3 dimensions doesn't make sense")

    @x.setter
    def x(self, val):
        if self.is_cartesian():
            self.__setitem__(0, val)
        else:
            raise IndexError("Setting the x component of a vector that is not in 2 or 3 dimensions doesn't make sense")

    @property
    def y(self):
        """ The y component of a cartesian vector

        :Raises:

        IndexError
            If the vector is not cartesian (i.e. if it is not two or three dimensional)

        """
        if self.is_cartesian():
            return self.__getitem__(1)
        else:
            raise IndexError("Accessing the y component of a vector that is not in 2 or 3 dimensions doesn't make sense")

    @y.setter
    def y(self, val):
        if self.is_cartesian():
            self.__setitem__(1, val)
        else:
            raise IndexError("Setting the y component of a vector that is not in 2 or 3 dimensions doesn't make sense")

    @property
    def z(self):
        """ The z component of a cartesian vector

        :Raises:

        IndexError
            If the vector is not three dimensional

        """
        if self.size() == 3:
            return self.__getitem__(2)
        else:
            raise IndexError("Accessing the z component of a vector that is not in 3 dimensions doesn't make sense")

    @z.setter
    def z(self, val):
        if self.size() == 3:
            self.__setitem__(2, val)
        else:
            raise IndexError("Setting the z component of a vector that is not 3 dimensions doesn't make sense")


    @property
    def column(self):
        """ Return the column-vector version of `self`, which is currently implemented as a `Matrix` object.
        Note that matrix-vector multiplication will work without getting `column` first, so don't use this unless
        you have a good reason to do so.

        """
        return self.reshape((self.size(), 1))

    ###########
    # Methods #
    ###########

    def reshape(self, newshape, order='C'):
        """ Overrides the `numpy.ndarray` reshape function to make things make sense.

        :Examples:


        >>> v = Vector([1,2,3,4])
        >>> v.reshape((2,2))
        Matrix([[ 1.,  2.],
                [ 3.,  4.]])
        >>> # `v` is unchanged:
        >>> v
        Vector([ 1.,  2.,  3.,  4.])
        >>> v.reshape((4,1))
        Matrix([[ 1.],
                [ 2.],
                [ 3.],
                [ 4.]])
        >>> import numpy as np
        >>> np.reshape(v, (2,2))
        Matrix([[ 1.,  2.],
                [ 3.,  4.]])


        """
        if len(newshape) == 1:
            return np.ndarray.reshape(self, newshape, order)
        elif len(newshape) == 2:
            return np.reshape(self.view(Matrix), newshape, order)
        else:
            return np.reshape(self.view(Tensor), newshape, order)


    #-----------------------------#
    # Methods that return Vectors #
    #-----------------------------#


    def normalized(self):
        """ Same as `normalize()`, but does not modify self
        """
        ret_val = copy(self)
        ret_val.normalize()
        return ret_val

    def cross(self, other):
        """ Returns cross product of self with other.  Wrapper to `numpy.cross`
        """
        return cross(self, other)

    #------------------------------#
    # Methods that return Matrices #
    #------------------------------#

    # This is a bad idea!!!
    #def transpose(self):
    #    """ Returns a transposed version of self (as a `Matrix`).
    #    This simply returns the vector "viewed" as a column vector in matrix form.  Synonym for the property `Vector.column`

    #    .. WARNING::
    #        Since `Vector` is a subclass of `numpy.ndarray`, and the transpose property `numpy.ndarray.T`
    #        of `numpy.ndarray` calls the first method named `transpose` in the class's `__mro__` list, this will get called
    #        by the `T` property.  The default behavior of the `T` operator for one-dimensional `numpy.ndarray` objects,
    #        however, is to do nothing.  Therefore, you should not use the `Vector` class for anything that may depend on
    #        the default behavior in `numpy.ndarray`.  The `T` behavior of `Vector` may change in the future to return
    #        a row vector.  Thus, do not call this from outside of `grendel.math` unless you have a good reason to do so.
    #        The `Vector.column` property should almost always give you what you need.

    #    Examples
    #    --------

    #    >>> Vector(1.0, 2.0, 3.0).transpose()
    #    Matrix([[ 1.],
    #            [ 2.],
    #            [ 3.]])

    #    """
    #    return self.column

Vector.x_axis = Vector([1, 0, 0])
Vector.y_axis = Vector([0, 1, 0])
Vector.z_axis = Vector([0, 0, 1])



#####################
# Dependent Imports #
#####################

from grendel.gmath.matrix import Matrix
from grendel.gmath import MatrixMultiplicationWarning

