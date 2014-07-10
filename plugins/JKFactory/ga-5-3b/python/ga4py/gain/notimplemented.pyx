
import sys

import numpy as np
cimport numpy as np


def docstring(obj, docstring):
    """docstring(obj, docstring)

    Add a docstring to a built-in obj if possible.
    If the obj already has a docstring raise a RuntimeError
    If this routine does not know how to add a docstring to the object
    raise a TypeError
    
    """
    # BUILTIN
    raise NotImplementedError


def add_newdoc(place, obj, doc):
    """Adds documentation to obj which is in module place.

    If doc is a string add it to obj as a docstring

    If doc is a tuple, then the first element is interpreted as
       an attribute of obj and the second as the docstring
          (method, docstring)

    If doc is a list, then each element of the list should be a
       sequence of length two --> [(method1, docstring1),
       (method2, docstring2), ...]

    This routine never raises an error.
    
    """
    raise NotImplementedError


def all(a, axis=None, out=None):
    """Test whether all array elements along a given axis evaluate to True.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which a logical AND is performed.
        The default (`axis` = `None`) is to perform a logical AND
        over a flattened input array.  `axis` may be negative, in which
        case it counts from the last to the first axis.
    out : ndarray, optional
        Alternate output array in which to place the result.
        It must have the same shape as the expected output and its
        type is preserved (e.g., if ``dtype(out)`` is float, the result
        will consist of 0.0's and 1.0's).  See `doc.ufuncs` (Section
        "Output arguments") for more details.

    Returns
    -------
    all : ndarray, bool
        A new boolean or array is returned unless `out` is specified,
        in which case a reference to `out` is returned.

    See Also
    --------
    ndarray.all : equivalent method

    any : Test whether any element along a given axis evaluates to True.

    Notes
    -----
    Not a Number (NaN), positive infinity and negative infinity
    evaluate to `True` because these are not equal to zero.

    Examples
    --------
    >>> np.all([[True,False],[True,True]])
    False

    >>> np.all([[True,False],[True,True]], axis=0)
    array([ True, False], dtype=bool)

    >>> np.all([-1, 4, 5])
    True

    >>> np.all([1.0, np.nan])
    True

    >>> o=np.array([False])
    >>> z=np.all([-1, 4, 5], out=o)
    >>> id(z), id(o), z                             # doctest: +SKIP
    (28293632, 28293632, array([ True], dtype=bool))
    
    """
    raise NotImplementedError


def allclose(a, b, rtol=1e-05, atol=1e-08):
    """Returns True if two arrays are element-wise equal within a tolerance.

    The tolerance values are positive, typically very small numbers.  The
    relative difference (`rtol` * abs(`b`)) and the absolute difference
    `atol` are added together to compare against the absolute difference
    between `a` and `b`.

    Parameters
    ----------
    a, b : array_like
        Input arrays to compare.
    rtol : float
        The relative tolerance parameter (see Notes).
    atol : float
        The absolute tolerance parameter (see Notes).

    Returns
    -------
    y : bool
        Returns True if the two arrays are equal within the given
        tolerance; False otherwise. If either array contains NaN, then
        False is returned.

    See Also
    --------
    all, any, alltrue, sometrue

    Notes
    -----
    If the following equation is element-wise True, then allclose returns
    True.

     absolute(`a` - `b`) <= (`atol` + `rtol` * absolute(`b`))

    The above equation is not symmetric in `a` and `b`, so that
    `allclose(a, b)` might be different from `allclose(b, a)` in
    some rare cases.

    Examples
    --------
    >>> np.allclose([1e10,1e-7], [1.00001e10,1e-8])
    False
    >>> np.allclose([1e10,1e-8], [1.00001e10,1e-9])
    True
    >>> np.allclose([1e10,1e-8], [1.0001e10,1e-9])
    False
    >>> np.allclose([1.0, np.nan], [1.0, np.nan])
    False
    
    """
    raise NotImplementedError


def alltrue(a, axis=None, out=None):
    """Check if all elements of input array are true.

    See Also
    --------
    numpy.all : Equivalent function; see for details.
    
    """
    raise NotImplementedError


def alterdot():
    """Change `dot`, `vdot`, and `innerproduct` to use accelerated BLAS functions.

    Typically, as a user of Numpy, you do not explicitly call this function. If
    Numpy is built with an accelerated BLAS, this function is automatically
    called when Numpy is imported.

    When Numpy is built with an accelerated BLAS like ATLAS, these functions
    are replaced to make use of the faster implementations.  The faster
    implementations only affect float32, float64, complex64, and complex128
    arrays. Furthermore, the BLAS API only includes matrix-matrix,
    matrix-vector, and vector-vector products. Products of arrays with larger
    dimensionalities use the built in functions and are not accelerated.

    See Also
    --------
    restoredot : `restoredot` undoes the effects of `alterdot`.
    
    """
    # BUILTIN
    raise NotImplementedError


def amax(a, axis=None, out=None):
    """Return the maximum of an array or maximum along an axis.

    Parameters
    ----------
    a : array_like
        Input data.
    axis : int, optional
        Axis along which to operate.  By default flattened input is used.
    out : ndarray, optional
        Alternate output array in which to place the result.  Must be of
        the same shape and buffer length as the expected output.  See
        `doc.ufuncs` (Section "Output arguments") for more details.

    Returns
    -------
    amax : ndarray or scalar
        Maximum of `a`. If `axis` is None, the result is a scalar value.
        If `axis` is given, the result is an array of dimension
        ``a.ndim - 1``.

    See Also
    --------
    nanmax : NaN values are ignored instead of being propagated.
    fmax : same behavior as the C99 fmax function.
    argmax : indices of the maximum values.

    Notes
    -----
    NaN values are propagated, that is if at least one item is NaN, the
    corresponding max value will be NaN as well.  To ignore NaN values
    (MATLAB behavior), please use nanmax.

    Examples
    --------
    >>> a = np.arange(4).reshape((2,2))
    >>> a
    array([[0, 1],
           [2, 3]])
    >>> np.amax(a)
    3
    >>> np.amax(a, axis=0)
    array([2, 3])
    >>> np.amax(a, axis=1)
    array([1, 3])

    >>> b = np.arange(5, dtype=np.float)
    >>> b[2] = np.NaN
    >>> np.amax(b)
    nan
    >>> np.nanmax(b)
    4.0
    
    """
    raise NotImplementedError


def amin(a, axis=None, out=None):
    """Return the minimum of an array or minimum along an axis.

    Parameters
    ----------
    a : array_like
        Input data.
    axis : int, optional
        Axis along which to operate.  By default a flattened input is used.
    out : ndarray, optional
        Alternative output array in which to place the result.  Must
        be of the same shape and buffer length as the expected output.
        See `doc.ufuncs` (Section "Output arguments") for more details.

    Returns
    -------
    amin : ndarray
        A new array or a scalar array with the result.

    See Also
    --------
    nanmin: nan values are ignored instead of being propagated
    fmin: same behavior as the C99 fmin function
    argmin: Return the indices of the minimum values.

    amax, nanmax, fmax

    Notes
    -----
    NaN values are propagated, that is if at least one item is nan, the
    corresponding min value will be nan as well. To ignore NaN values (matlab
    behavior), please use nanmin.

    Examples
    --------
    >>> a = np.arange(4).reshape((2,2))
    >>> a
    array([[0, 1],
           [2, 3]])
    >>> np.amin(a)           # Minimum of the flattened array
    0
    >>> np.amin(a, axis=0)         # Minima along the first axis
    array([0, 1])
    >>> np.amin(a, axis=1)         # Minima along the second axis
    array([0, 2])

    >>> b = np.arange(5, dtype=np.float)
    >>> b[2] = np.NaN
    >>> np.amin(b)
    nan
    >>> np.nanmin(b)
    0.0
    
    """
    raise NotImplementedError


def angle(z, deg=0):
    """Return the angle of the complex argument.

    Parameters
    ----------
    z : array_like
        A complex number or sequence of complex numbers.
    deg : bool, optional
        Return angle in degrees if True, radians if False (default).

    Returns
    -------
    angle : {ndarray, scalar}
        The counterclockwise angle from the positive real axis on
        the complex plane, with dtype as numpy.float64.

    See Also
    --------
    arctan2
    absolute



    Examples
    --------
    >>> np.angle([1.0, 1.0j, 1+1j])               # in radians
    array([ 0.        ,  1.57079633,  0.78539816])
    >>> np.angle(1+1j, deg=True)                  # in degrees
    45.0
    
    """
    raise NotImplementedError


def any(a, axis=None, out=None):
    """Test whether any array element along a given axis evaluates to True.

    Returns single boolean unless `axis` is not ``None``

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which a logical OR is performed.  The default
        (`axis` = `None`) is to perform a logical OR over a flattened
        input array. `axis` may be negative, in which case it counts
        from the last to the first axis.
    out : ndarray, optional
        Alternate output array in which to place the result.  It must have
        the same shape as the expected output and its type is preserved
        (e.g., if it is of type float, then it will remain so, returning
        1.0 for True and 0.0 for False, regardless of the type of `a`).
        See `doc.ufuncs` (Section "Output arguments") for details.

    Returns
    -------
    any : bool or ndarray
        A new boolean or `ndarray` is returned unless `out` is specified,
        in which case a reference to `out` is returned.

    See Also
    --------
    ndarray.any : equivalent method

    all : Test whether all elements along a given axis evaluate to True.

    Notes
    -----
    Not a Number (NaN), positive infinity and negative infinity evaluate
    to `True` because these are not equal to zero.

    Examples
    --------
    >>> np.any([[True, False], [True, True]])
    True

    >>> np.any([[True, False], [False, False]], axis=0)
    array([ True, False], dtype=bool)

    >>> np.any([-1, 0, 5])
    True

    >>> np.any(np.nan)
    True

    >>> o=np.array([False])
    >>> z=np.any([-1, 4, 5], out=o)
    >>> z, o
    (array([ True], dtype=bool), array([ True], dtype=bool))
    >>> # Check now that z is a reference to o
    >>> z is o
    True
    >>> id(z), id(o) # identity of z and o              # doctest: +SKIP
    (191614240, 191614240)
    
    """
    raise NotImplementedError


def append(arr, values, axis=None):
    """Append values to the end of an array.

    Parameters
    ----------
    arr : array_like
        Values are appended to a copy of this array.
    values : array_like
        These values are appended to a copy of `arr`.  It must be of the
        correct shape (the same shape as `arr`, excluding `axis`).  If `axis`
        is not specified, `values` can be any shape and will be flattened
        before use.
    axis : int, optional
        The axis along which `values` are appended.  If `axis` is not given,
        both `arr` and `values` are flattened before use.

    Returns
    -------
    out : ndarray
        A copy of `arr` with `values` appended to `axis`.  Note that `append`
        does not occur in-place: a new array is allocated and filled.  If
        `axis` is None, `out` is a flattened array.

    See Also
    --------
    insert : Insert elements into an array.
    delete : Delete elements from an array.

    Examples
    --------
    >>> np.append([1, 2, 3], [[4, 5, 6], [7, 8, 9]])
    array([1, 2, 3, 4, 5, 6, 7, 8, 9])

    When `axis` is specified, `values` must have the correct shape.

    >>> np.append([[1, 2, 3], [4, 5, 6]], [[7, 8, 9]], axis=0)
    array([[1, 2, 3],
           [4, 5, 6],
           [7, 8, 9]])
    >>> np.append([[1, 2, 3], [4, 5, 6]], [7, 8, 9], axis=0)
    Traceback (most recent call last):
    ...
    ValueError: arrays must have same number of dimensions
    
    """
    raise NotImplementedError


def apply_along_axis(func1d, axis, arr):
    """Apply a function to 1-D slices along the given axis.

    Execute `func1d(a, *args)` where `func1d` operates on 1-D arrays and `a`
    is a 1-D slice of `arr` along `axis`.

    Parameters
    ----------
    func1d : function
        This function should accept 1-D arrays. It is applied to 1-D
        slices of `arr` along the specified axis.
    axis : integer
        Axis along which `arr` is sliced.
    arr : ndarray
        Input array.
    args : any
        Additional arguments to `func1d`.

    Returns
    -------
    outarr : ndarray
        The output array. The shape of `outarr` is identical to the shape of
        `arr`, except along the `axis` dimension, where the length of `outarr`
        is equal to the size of the return value of `func1d`.  If `func1d`
        returns a scalar `outarr` will have one fewer dimensions than `arr`.

    See Also
    --------
    apply_over_axes : Apply a function repeatedly over multiple axes.

    Examples
    --------
    >>> def my_func(a):
    ...     '''Average first and last element of a 1-D array'''
    ...     return (a[0] + a[-1]) * 0.5
    >>> b = np.array([[1,2,3], [4,5,6], [7,8,9]])
    >>> np.apply_along_axis(my_func, 0, b)
    array([ 4.,  5.,  6.])
    >>> np.apply_along_axis(my_func, 1, b)
    array([ 2.,  5.,  8.])

    For a function that doesn't return a scalar, the number of dimensions in
    `outarr` is the same as `arr`.

    >>> def new_func(a):
    ...     '''Divide elements of a by 2.'''
    ...     return a * 0.5
    >>> b = np.array([[1,2,3], [4,5,6], [7,8,9]])
    >>> np.apply_along_axis(new_func, 0, b)
    array([[ 0.5,  1. ,  1.5],
           [ 2. ,  2.5,  3. ],
           [ 3.5,  4. ,  4.5]])
    
    """
    raise NotImplementedError


def apply_over_axes(func, a, axes):
    """Apply a function repeatedly over multiple axes.

    `func` is called as `res = func(a, axis)`, where `axis` is the first
    element of `axes`.  The result `res` of the function call must have
    either the same dimensions as `a` or one less dimension.  If `res`
    has one less dimension than `a`, a dimension is inserted before
    `axis`.  The call to `func` is then repeated for each axis in `axes`,
    with `res` as the first argument.

    Parameters
    ----------
    func : function
        This function must take two arguments, `func(a, axis)`.
    a : array_like
        Input array.
    axes : array_like
        Axes over which `func` is applied; the elements must be integers.

    Returns
    -------
    val : ndarray
        The output array.  The number of dimensions is the same as `a`,
        but the shape can be different.  This depends on whether `func`
        changes the shape of its output with respect to its input.

    See Also
    --------
    apply_along_axis :
        Apply a function to 1-D slices of an array along the given axis.

    Examples
    --------
    >>> a = np.arange(24).reshape(2,3,4)
    >>> a
    array([[[ 0,  1,  2,  3],
            [ 4,  5,  6,  7],
            [ 8,  9, 10, 11]],
           [[12, 13, 14, 15],
            [16, 17, 18, 19],
            [20, 21, 22, 23]]])

    Sum over axes 0 and 2. The result has same number of dimensions
    as the original array:

    >>> np.apply_over_axes(np.sum, a, [0,2])
    array([[[ 60],
            [ 92],
            [124]]])
    
    """
    raise NotImplementedError


def argmax(a, axis=None):
    """Indices of the maximum values along an axis.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : int, optional
        By default, the index is into the flattened array, otherwise
        along the specified axis.

    Returns
    -------
    index_array : ndarray of ints
        Array of indices into the array. It has the same shape as `a.shape`
        with the dimension along `axis` removed.

    See Also
    --------
    ndarray.argmax, argmin
    amax : The maximum value along a given axis.
    unravel_index : Convert a flat index into an index tuple.

    Notes
    -----
    In case of multiple occurrences of the maximum values, the indices
    corresponding to the first occurrence are returned.

    Examples
    --------
    >>> a = np.arange(6).reshape(2,3)
    >>> a
    array([[0, 1, 2],
           [3, 4, 5]])
    >>> np.argmax(a)
    5
    >>> np.argmax(a, axis=0)
    array([1, 1, 1])
    >>> np.argmax(a, axis=1)
    array([2, 2])

    >>> b = np.arange(6)
    >>> b[1] = 5
    >>> b
    array([0, 5, 2, 3, 4, 5])
    >>> np.argmax(b) # Only the first occurrence is returned.
    1
    
    """
    raise NotImplementedError


def argmin(a, axis=None):
    """Return the indices of the minimum values along an axis.

    See Also
    --------
    argmax : Similar function.  Please refer to `numpy.argmax` for detailed
        documentation.
    
    """
    raise NotImplementedError


def argsort(a, axis=-1, kind='quicksort', order=None):
    """Returns the indices that would sort an array.

    Perform an indirect sort along the given axis using the algorithm specified
    by the `kind` keyword. It returns an array of indices of the same shape as
    `a` that index data along the given axis in sorted order.

    Parameters
    ----------
    a : array_like
        Array to sort.
    axis : int or None, optional
        Axis along which to sort.  The default is -1 (the last axis). If None,
        the flattened array is used.
    kind : {'quicksort', 'mergesort', 'heapsort'}, optional
        Sorting algorithm.
    order : list, optional
        When `a` is an array with fields defined, this argument specifies
        which fields to compare first, second, etc.  Not all fields need be
        specified.

    Returns
    -------
    index_array : ndarray, int
        Array of indices that sort `a` along the specified axis.
        In other words, ``a[index_array]`` yields a sorted `a`.

    See Also
    --------
    sort : Describes sorting algorithms used.
    lexsort : Indirect stable sort with multiple keys.
    ndarray.sort : Inplace sort.

    Notes
    -----
    See `sort` for notes on the different sorting algorithms.

    As of NumPy 1.4.0 `argsort` works with real/complex arrays containing
    nan values. The enhanced sort order is documented in `sort`.

    Examples
    --------
    One dimensional array:

    >>> x = np.array([3, 1, 2])
    >>> np.argsort(x)
    array([1, 2, 0])

    Two-dimensional array:

    >>> x = np.array([[0, 3], [2, 2]])
    >>> x
    array([[0, 3],
           [2, 2]])

    >>> np.argsort(x, axis=0)
    array([[0, 1],
           [1, 0]])

    >>> np.argsort(x, axis=1)
    array([[0, 1],
           [0, 1]])

    Sorting with keys:

    >>> x = np.array([(1, 0), (0, 1)], dtype=[('x', '<i4'), ('y', '<i4')])
    >>> x
    array([(1, 0), (0, 1)],
          dtype=[('x', '<i4'), ('y', '<i4')])

    >>> np.argsort(x, order=('x','y'))
    array([1, 0])

    >>> np.argsort(x, order=('y','x'))
    array([0, 1])
    
    """
    raise NotImplementedError


def argwhere(a):
    """Find the indices of array elements that are non-zero, grouped by element.

    Parameters
    ----------
    a : array_like
        Input data.

    Returns
    -------
    index_array : ndarray
        Indices of elements that are non-zero. Indices are grouped by element.

    See Also
    --------
    where, nonzero

    Notes
    -----
    ``np.argwhere(a)`` is the same as ``np.transpose(np.nonzero(a))``.

    The output of ``argwhere`` is not suitable for indexing arrays.
    For this purpose use ``where(a)`` instead.

    Examples
    --------
    >>> x = np.arange(6).reshape(2,3)
    >>> x
    array([[0, 1, 2],
           [3, 4, 5]])
    >>> np.argwhere(x>1)
    array([[0, 2],
           [1, 0],
           [1, 1],
           [1, 2]])
    
    """
    raise NotImplementedError


def around(a, decimals=0, out=None):
    """Evenly round to the given number of decimals.

    Parameters
    ----------
    a : array_like
        Input data.
    decimals : int, optional
        Number of decimal places to round to (default: 0).  If
        decimals is negative, it specifies the number of positions to
        the left of the decimal point.
    out : ndarray, optional
        Alternative output array in which to place the result. It must have
        the same shape as the expected output, but the type of the output
        values will be cast if necessary. See `doc.ufuncs` (Section
        "Output arguments") for details.

    Returns
    -------
    rounded_array : ndarray
        An array of the same type as `a`, containing the rounded values.
        Unless `out` was specified, a new array is created.  A reference to
        the result is returned.

        The real and imaginary parts of complex numbers are rounded
        separately.  The result of rounding a float is a float.

    See Also
    --------
    ndarray.round : equivalent method

    ceil, fix, floor, rint, trunc


    Notes
    -----
    For values exactly halfway between rounded decimal values, Numpy
    rounds to the nearest even value. Thus 1.5 and 2.5 round to 2.0,
    -0.5 and 0.5 round to 0.0, etc. Results may also be surprising due
    to the inexact representation of decimal fractions in the IEEE
    floating point standard [1]_ and errors introduced when scaling
    by powers of ten.

    References
    ----------
    .. [1] "Lecture Notes on the Status of  IEEE 754", William Kahan,
           http://www.cs.berkeley.edu/~wkahan/ieee754status/IEEE754.PDF
    .. [2] "How Futile are Mindless Assessments of
           Roundoff in Floating-Point Computation?", William Kahan,
           http://www.cs.berkeley.edu/~wkahan/Mindless.pdf

    Examples
    --------
    >>> np.around([0.37, 1.64])
    array([ 0.,  2.])
    >>> np.around([0.37, 1.64], decimals=1)
    array([ 0.4,  1.6])
    >>> np.around([.5, 1.5, 2.5, 3.5, 4.5]) # rounds to nearest even value
    array([ 0.,  2.,  2.,  4.,  4.])
    >>> np.around([1,2,3,11], decimals=1) # ndarray of ints is returned
    array([ 1,  2,  3, 11])
    >>> np.around([1,2,3,11], decimals=-1)
    array([ 0,  0,  0, 10])
    
    """
    raise NotImplementedError


def array(object, dtype=None, copy=True, order=None, subok=False, ndmin=0):
    """array(object, dtype=None, copy=True, order=None, subok=False, ndmin=0)

    Create an array.

    Parameters
    ----------
    object : array_like
        An array, any object exposing the array interface, an
        object whose __array__ method returns an array, or any
        (nested) sequence.
    dtype : data-type, optional
        The desired data-type for the array.  If not given, then
        the type will be determined as the minimum type required
        to hold the objects in the sequence.  This argument can only
        be used to 'upcast' the array.  For downcasting, use the
        .astype(t) method.
    copy : bool, optional
        If true (default), then the object is copied.  Otherwise, a copy
        will only be made if __array__ returns a copy, if obj is a
        nested sequence, or if a copy is needed to satisfy any of the other
        requirements (`dtype`, `order`, etc.).
    order : {'C', 'F', 'A'}, optional
        Specify the order of the array.  If order is 'C' (default), then the
        array will be in C-contiguous order (last-index varies the
        fastest).  If order is 'F', then the returned array
        will be in Fortran-contiguous order (first-index varies the
        fastest).  If order is 'A', then the returned array may
        be in any order (either C-, Fortran-contiguous, or even
        discontiguous).
    subok : bool, optional
        If True, then sub-classes will be passed-through, otherwise
        the returned array will be forced to be a base-class array (default).
    ndmin : int, optional
        Specifies the minimum number of dimensions that the resulting
        array should have.  Ones will be pre-pended to the shape as
        needed to meet this requirement.

    Returns
    -------
    out : ndarray
        An array object satisfying the specified requirements.

    See Also
    --------
    empty, empty_like, zeros, zeros_like, ones, ones_like, fill

    Examples
    --------
    >>> np.array([1, 2, 3])
    array([1, 2, 3])

    Upcasting:

    >>> np.array([1, 2, 3.0])
    array([ 1.,  2.,  3.])

    More than one dimension:

    >>> np.array([[1, 2], [3, 4]])
    array([[1, 2],
           [3, 4]])

    Minimum dimensions 2:

    >>> np.array([1, 2, 3], ndmin=2)
    array([[1, 2, 3]])

    Type provided:

    >>> np.array([1, 2, 3], dtype=complex)
    array([ 1.+0.j,  2.+0.j,  3.+0.j])

    Data-type consisting of more than one element:

    >>> x = np.array([(1,2),(3,4)],dtype=[('a','<i4'),('b','<i4')])
    >>> x['a']
    array([1, 3])

    Creating an array from sub-classes:

    >>> np.array(np.mat('1 2; 3 4'))
    array([[1, 2],
           [3, 4]])

    >>> np.array(np.mat('1 2; 3 4'), subok=True)
    matrix([[1, 2],
            [3, 4]])
    
    """
    # BUILTIN
    raise NotImplementedError


def array2string(a, max_line_width=None, precision=None, suppress_small=None, separator=' ', prefix='', style=repr):
    """Return a string representation of an array.

    Parameters
    ----------
    a : ndarray
        Input array.
    max_line_width : int, optional
        The maximum number of columns the string should span. Newline
        characters splits the string appropriately after array elements.
    precision : int, optional
        Floating point precision. Default is the current printing
        precision (usually 8), which can be altered using `set_printoptions`.
    suppress_small : bool, optional
        Represent very small numbers as zero. A number is "very small" if it
        is smaller than the current printing precision.
    separator : str, optional
        Inserted between elements.
    prefix : str, optional
        An array is typically printed as::

          'prefix(' + array2string(a) + ')'

        The length of the prefix string is used to align the
        output correctly.
    style : function, optional
        A function that accepts an ndarray and returns a string.  Used only
        when the shape of `a` is equal to ().

    Returns
    -------
    array_str : str
        String representation of the array.

    See Also
    --------
    array_str, array_repr, set_printoptions

    Examples
    --------
    >>> x = np.array([1e-16,1,2,3])
    >>> print np.array2string(x, precision=2, separator=',',
    ...                       suppress_small=True)
    [ 0., 1., 2., 3.]
    
    """
    raise NotImplementedError


def array_equal(a1, a2):
    """True if two arrays have the same shape and elements, False otherwise.

    Parameters
    ----------
    a1, a2 : array_like
        Input arrays.

    Returns
    -------
    b : bool
        Returns True if the arrays are equal.

    See Also
    --------
    allclose: Returns True if two arrays are element-wise equal within a
              tolerance.
    array_equiv: Returns True if input arrays are shape consistent and all
                 elements equal.

    Examples
    --------
    >>> np.array_equal([1, 2], [1, 2])
    True
    >>> np.array_equal(np.array([1, 2]), np.array([1, 2]))
    True
    >>> np.array_equal([1, 2], [1, 2, 3])
    False
    >>> np.array_equal([1, 2], [1, 4])
    False
    
    """
    raise NotImplementedError


def array_equiv(a1, a2):
    """Returns True if input arrays are shape consistent and all elements equal.

    Shape consistent means they are either the same shape, or one input array
    can be broadcasted to create the same shape as the other one.

    Parameters
    ----------
    a1, a2 : array_like
        Input arrays.

    Returns
    -------
    out : bool
        True if equivalent, False otherwise.

    Examples
    --------
    >>> np.array_equiv([1, 2], [1, 2])
    True
    >>> np.array_equiv([1, 2], [1, 3])
    False

    Showing the shape equivalence:

    >>> np.array_equiv([1, 2], [[1, 2], [1, 2]])
    True
    >>> np.array_equiv([1, 2], [[1, 2, 1, 2], [1, 2, 1, 2]])
    False

    >>> np.array_equiv([1, 2], [[1, 2], [1, 3]])
    False
    
    """
    raise NotImplementedError


def array_repr(arr, max_line_width=None, precision=None, suppress_small=None):
    """Return the string representation of an array.

    Parameters
    ----------
    arr : ndarray
        Input array.
    max_line_width : int, optional
        The maximum number of columns the string should span. Newline
        characters split the string appropriately after array elements.
    precision : int, optional
        Floating point precision. Default is the current printing precision
        (usually 8), which can be altered using `set_printoptions`.
    suppress_small : bool, optional
        Represent very small numbers as zero, default is False. Very small
        is defined by `precision`, if the precision is 8 then
        numbers smaller than 5e-9 are represented as zero.

    Returns
    -------
    string : str
      The string representation of an array.

    See Also
    --------
    array_str, array2string, set_printoptions

    Examples
    --------
    >>> np.array_repr(np.array([1,2]))
    'array([1, 2])'
    >>> np.array_repr(np.ma.array([0.]))
    'MaskedArray([ 0.])'
    >>> np.array_repr(np.array([], np.int32))
    'array([], dtype=int32)'

    >>> x = np.array([1e-6, 4e-7, 2, 3])
    >>> np.array_repr(x, precision=6, suppress_small=True)
    'array([ 0.000001,  0.      ,  2.      ,  3.      ])'
    
    """
    raise NotImplementedError


def array_split(ary, indices_or_sections, axis=0):
    """Split an array into multiple sub-arrays of equal or near-equal size.

    Please refer to the ``split`` documentation.  The only difference
    between these functions is that ``array_split`` allows
    `indices_or_sections` to be an integer that does *not* equally
    divide the axis.

    See Also
    --------
    split : Split array into multiple sub-arrays of equal size.

    Examples
    --------
    >>> x = np.arange(8.0)
    >>> np.array_split(x, 3)
        [array([ 0.,  1.,  2.]), array([ 3.,  4.,  5.]), array([ 6.,  7.])]
    
    """
    raise NotImplementedError


def array_str(a, max_line_width=None, precision=None, suppress_small=None):
    """Return a string representation of the data in an array.

    The data in the array is returned as a single string.  This function is
    similar to `array_repr`, the difference being that `array_repr` also
    returns information on the kind of array and its data type.

    Parameters
    ----------
    a : ndarray
        Input array.
    max_line_width : int, optional
        Inserts newlines if text is longer than `max_line_width`.  The
        default is, indirectly, 75.
    precision : int, optional
        Floating point precision.  Default is the current printing precision
        (usually 8), which can be altered using `set_printoptions`.
    suppress_small : bool, optional
        Represent numbers "very close" to zero as zero; default is False.
        Very close is defined by precision: if the precision is 8, e.g.,
        numbers smaller (in absolute value) than 5e-9 are represented as
        zero.

    See Also
    --------
    array2string, array_repr, set_printoptions

    Examples
    --------
    >>> np.array_str(np.arange(3))
    '[0 1 2]'
    
    """
    raise NotImplementedError


def asanyarray(a, dtype=None, order=None):
    """Convert the input to an ndarray, but pass ndarray subclasses through.

    Parameters
    ----------
    a : array_like
        Input data, in any form that can be converted to an array.  This
        includes scalars, lists, lists of tuples, tuples, tuples of tuples,
        tuples of lists, and ndarrays.
    dtype : data-type, optional
        By default, the data-type is inferred from the input data.
    order : {'C', 'F'}, optional
        Whether to use row-major ('C') or column-major ('F') memory
        representation.  Defaults to 'C'.

    Returns
    -------
    out : ndarray or an ndarray subclass
        Array interpretation of `a`.  If `a` is an ndarray or a subclass
        of ndarray, it is returned as-is and no copy is performed.

    See Also
    --------
    asarray : Similar function which always returns ndarrays.
    ascontiguousarray : Convert input to a contiguous array.
    asfarray : Convert input to a floating point ndarray.
    asfortranarray : Convert input to an ndarray with column-major
                     memory order.
    asarray_chkfinite : Similar function which checks input for NaNs and
                        Infs.
    fromiter : Create an array from an iterator.
    fromfunction : Construct an array by executing a function on grid
                   positions.

    Examples
    --------
    Convert a list into an array:

    >>> a = [1, 2]
    >>> np.asanyarray(a)
    array([1, 2])

    Instances of `ndarray` subclasses are passed through as-is:

    >>> a = np.matrix([1, 2])
    >>> np.asanyarray(a) is a
    True
    
    """
    raise NotImplementedError


def asarray_chkfinite(a):
    """Convert the input to an array, checking for NaNs or Infs.

    Parameters
    ----------
    a : array_like
        Input data, in any form that can be converted to an array.  This
        includes lists, lists of tuples, tuples, tuples of tuples, tuples
        of lists and ndarrays.  Success requires no NaNs or Infs.
    dtype : data-type, optional
        By default, the data-type is inferred from the input data.
    order : {'C', 'F'}, optional
        Whether to use row-major ('C') or column-major ('FORTRAN') memory
        representation.  Defaults to 'C'.

    Returns
    -------
    out : ndarray
        Array interpretation of `a`.  No copy is performed if the input
        is already an ndarray.  If `a` is a subclass of ndarray, a base
        class ndarray is returned.

    Raises
    ------
    ValueError
        Raises ValueError if `a` contains NaN (Not a Number) or Inf (Infinity).

    See Also
    --------
    asarray : Create and array.
    asanyarray : Similar function which passes through subclasses.
    ascontiguousarray : Convert input to a contiguous array.
    asfarray : Convert input to a floating point ndarray.
    asfortranarray : Convert input to an ndarray with column-major
                     memory order.
    fromiter : Create an array from an iterator.
    fromfunction : Construct an array by executing a function on grid
                   positions.

    Examples
    --------
    Convert a list into an array.  If all elements are finite
    ``asarray_chkfinite`` is identical to ``asarray``.

    >>> a = [1, 2]
    >>> np.asarray_chkfinite(a)
    array([1, 2])

    Raises ValueError if array_like contains Nans or Infs.

    >>> a = [1, 2, np.inf]
    >>> try:
    ...     np.asarray_chkfinite(a)
    ... except ValueError:
    ...     print 'ValueError'
    ...
    ValueError
    
    """
    raise NotImplementedError


def ascontiguousarray(a, dtype=None):
    """Return a contiguous array in memory (C order).

    Parameters
    ----------
    a : array_like
        Input array.
    dtype : str or dtype object, optional
        Data-type of returned array.

    Returns
    -------
    out : ndarray
        Contiguous array of same shape and content as `a`, with type `dtype`
        if specified.

    See Also
    --------
    asfortranarray : Convert input to an ndarray with column-major
                     memory order.
    require : Return an ndarray that satisfies requirements.
    ndarray.flags : Information about the memory layout of the array.

    Examples
    --------
    >>> x = np.arange(6).reshape(2,3)
    >>> np.ascontiguousarray(x, dtype=np.float32)
    array([[ 0.,  1.,  2.],
           [ 3.,  4.,  5.]], dtype=float32)
    >>> x.flags['C_CONTIGUOUS']
    True
    
    """
    raise NotImplementedError


def asfarray(a, dtype=np.float64):
    """Return an array converted to a float type.

    Parameters
    ----------
    a : array_like
        The input array.
    dtype : str or dtype object, optional
        Float type code to coerce input array `a`.  If `dtype` is one of the
        'int' dtypes, it is replaced with float64.

    Returns
    -------
    out : ndarray
        The input `a` as a float ndarray.

    Examples
    --------
    >>> np.asfarray([2, 3])
    array([ 2.,  3.])
    >>> np.asfarray([2, 3], dtype='float')
    array([ 2.,  3.])
    >>> np.asfarray([2, 3], dtype='int8')
    array([ 2.,  3.])
    
    """
    raise NotImplementedError


def asfortranarray(a, dtype=None):
    """Return an array laid out in Fortran order in memory.

    Parameters
    ----------
    a : array_like
        Input array.
    dtype : str or dtype object, optional
        By default, the data-type is inferred from the input data.

    Returns
    -------
    out : ndarray
        The input `a` in Fortran, or column-major, order.

    See Also
    --------
    ascontiguousarray : Convert input to a contiguous (C order) array.
    asanyarray : Convert input to an ndarray with either row or
        column-major memory order.
    require : Return an ndarray that satisfies requirements.
    ndarray.flags : Information about the memory layout of the array.

    Examples
    --------
    >>> x = np.arange(6).reshape(2,3)
    >>> y = np.asfortranarray(x)
    >>> x.flags['F_CONTIGUOUS']
    False
    >>> y.flags['F_CONTIGUOUS']
    True
    
    """
    raise NotImplementedError


def asmatrix(data, dtype=None):
    """Interpret the input as a matrix.

    Unlike `matrix`, `asmatrix` does not make a copy if the input is already
    a matrix or an ndarray.  Equivalent to ``matrix(data, copy=False)``.

    Parameters
    ----------
    data : array_like
        Input data.

    Returns
    -------
    mat : matrix
        `data` interpreted as a matrix.

    Examples
    --------
    >>> x = np.array([[1, 2], [3, 4]])

    >>> m = np.asmatrix(x)

    >>> x[0,0] = 5

    >>> m
    matrix([[5, 2],
            [3, 4]])
    
    """
    raise NotImplementedError


def asscalar(a):
    """Convert an array of size 1 to its scalar equivalent.

    Parameters
    ----------
    a : ndarray
        Input array of size 1.

    Returns
    -------
    out : scalar
        Scalar representation of `a`. The input data type is preserved.

    Examples
    --------
    >>> np.asscalar(np.array([24]))
    24
    
    """
    raise NotImplementedError


def atleast_1d():
    """Convert inputs to arrays with at least one dimension.

    Scalar inputs are converted to 1-dimensional arrays, whilst
    higher-dimensional inputs are preserved.

    Parameters
    ----------
    array1, array2, ... : array_like
        One or more input arrays.

    Returns
    -------
    ret : ndarray
        An array, or sequence of arrays, each with ``a.ndim >= 1``.
        Copies are made only if necessary.

    See Also
    --------
    atleast_2d, atleast_3d

    Examples
    --------
    >>> np.atleast_1d(1.0)
    array([ 1.])

    >>> x = np.arange(9.0).reshape(3,3)
    >>> np.atleast_1d(x)
    array([[ 0.,  1.,  2.],
           [ 3.,  4.,  5.],
           [ 6.,  7.,  8.]])
    >>> np.atleast_1d(x) is x
    True

    >>> np.atleast_1d(1, [3, 4])
    [array([1]), array([3, 4])]
    
    """
    raise NotImplementedError


def atleast_2d():
    """View inputs as arrays with at least two dimensions.

    Parameters
    ----------
    array1, array2, ... : array_like
        One or more array-like sequences.  Non-array inputs are converted
        to arrays.  Arrays that already have two or more dimensions are
        preserved.

    Returns
    -------
    res, res2, ... : ndarray
        An array, or tuple of arrays, each with ``a.ndim >= 2``.
        Copies are avoided where possible, and views with two or more
        dimensions are returned.

    See Also
    --------
    atleast_1d, atleast_3d

    Examples
    --------
    >>> np.atleast_2d(3.0)
    array([[ 3.]])

    >>> x = np.arange(3.0)
    >>> np.atleast_2d(x)
    array([[ 0.,  1.,  2.]])
    >>> np.atleast_2d(x).base is x
    True

    >>> np.atleast_2d(1, [1, 2], [[1, 2]])
    [array([[1]]), array([[1, 2]]), array([[1, 2]])]
    
    """
    raise NotImplementedError


def atleast_3d():
    """View inputs as arrays with at least three dimensions.

    Parameters
    ----------
    array1, array2, ... : array_like
        One or more array-like sequences.  Non-array inputs are converted to
        arrays.  Arrays that already have three or more dimensions are
        preserved.

    Returns
    -------
    res1, res2, ... : ndarray
        An array, or tuple of arrays, each with ``a.ndim >= 3``.  Copies are
        avoided where possible, and views with three or more dimensions are
        returned.  For example, a 1-D array of shape ``(N,)`` becomes a view
        of shape ``(1, N, 1)``, and a 2-D array of shape ``(M, N)`` becomes a
        view of shape ``(M, N, 1)``.

    See Also
    --------
    atleast_1d, atleast_2d

    Examples
    --------
    >>> np.atleast_3d(3.0)
    array([[[ 3.]]])

    >>> x = np.arange(3.0)
    >>> np.atleast_3d(x).shape
    (1, 3, 1)

    >>> x = np.arange(12.0).reshape(4,3)
    >>> np.atleast_3d(x).shape
    (4, 3, 1)
    >>> np.atleast_3d(x).base is x
    True

    >>> for arr in np.atleast_3d([1, 2], [[1, 2]], [[[1, 2]]]):
    ...     print arr, arr.shape
    ...
    [[[1]
      [2]]] (1, 2, 1)
    [[[1]
      [2]]] (1, 2, 1)
    [[[1 2]]] (1, 1, 2)
    
    """
    raise NotImplementedError


def average(a, axis=None, weights=None, returned=False):
    """Compute the weighted average along the specified axis.

    Parameters
    ----------
    a : array_like
        Array containing data to be averaged. If `a` is not an array, a
        conversion is attempted.
    axis : int, optional
        Axis along which to average `a`. If `None`, averaging is done over
        the flattened array.
    weights : array_like, optional
        An array of weights associated with the values in `a`. Each value in
        `a` contributes to the average according to its associated weight.
        The weights array can either be 1-D (in which case its length must be
        the size of `a` along the given axis) or of the same shape as `a`.
        If `weights=None`, then all data in `a` are assumed to have a
        weight equal to one.
    returned : bool, optional
        Default is `False`. If `True`, the tuple (`average`, `sum_of_weights`)
        is returned, otherwise only the average is returned.
        If `weights=None`, `sum_of_weights` is equivalent to the number of
        elements over which the average is taken.


    Returns
    -------
    average, [sum_of_weights] : {array_type, double}
        Return the average along the specified axis. When returned is `True`,
        return a tuple with the average as the first element and the sum
        of the weights as the second element. The return type is `Float`
        if `a` is of integer type, otherwise it is of the same type as `a`.
        `sum_of_weights` is of the same type as `average`.

    Raises
    ------
    ZeroDivisionError
        When all weights along axis are zero. See `numpy.ma.average` for a
        version robust to this type of error.
    TypeError
        When the length of 1D `weights` is not the same as the shape of `a`
        along axis.

    See Also
    --------
    mean

    ma.average : average for masked arrays

    Examples
    --------
    >>> data = range(1,5)
    >>> data
    [1, 2, 3, 4]
    >>> np.average(data)
    2.5
    >>> np.average(range(1,11), weights=range(10,0,-1))
    4.0

    >>> data = np.arange(6).reshape((3,2))
    >>> data
    array([[0, 1],
           [2, 3],
           [4, 5]])
    >>> np.average(data, axis=1, weights=[1./4, 3./4])
    array([ 0.75,  2.75,  4.75])
    >>> np.average(data, weights=[1./4, 3./4])
    Traceback (most recent call last):
    ...
    TypeError: Axis must be specified when shapes of a and weights differ.
    
    """
    raise NotImplementedError


def bartlett(M):
    """Return the Bartlett window.

    The Bartlett window is very similar to a triangular window, except
    that the end points are at zero.  It is often used in signal
    processing for tapering a signal, without generating too much
    ripple in the frequency domain.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an
        empty array is returned.

    Returns
    -------
    out : array
        The triangular window, normalized to one (the value one
        appears only if the number of samples is odd), with the first
        and last samples equal to zero.

    See Also
    --------
    blackman, hamming, hanning, kaiser

    Notes
    -----
    The Bartlett window is defined as

    .. math:: w(n) = \frac{2}{M-1} \left(
              \frac{M-1}{2} - \left|n - \frac{M-1}{2}\right|
              \right)

    Most references to the Bartlett window come from the signal
    processing literature, where it is used as one of many windowing
    functions for smoothing values.  Note that convolution with this
    window produces linear interpolation.  It is also known as an
    apodization (which means"removing the foot", i.e. smoothing
    discontinuities at the beginning and end of the sampled signal) or
    tapering function. The fourier transform of the Bartlett is the product
    of two sinc functions.
    Note the excellent discussion in Kanasewich.

    References
    ----------
    .. [1] M.S. Bartlett, "Periodogram Analysis and Continuous Spectra",
           Biometrika 37, 1-16, 1950.
    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics",
           The University of Alberta Press, 1975, pp. 109-110.
    .. [3] A.V. Oppenheim and R.W. Schafer, "Discrete-Time Signal
           Processing", Prentice-Hall, 1999, pp. 468-471.
    .. [4] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function
    .. [5] W.H. Press,  B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
           "Numerical Recipes", Cambridge University Press, 1986, page 429.


    Examples
    --------
    >>> np.bartlett(12)
    array([ 0.        ,  0.18181818,  0.36363636,  0.54545455,  0.72727273,
            0.90909091,  0.90909091,  0.72727273,  0.54545455,  0.36363636,
            0.18181818,  0.        ])

    Plot the window and its frequency response (requires SciPy and matplotlib):

    >>> from numpy import clip, log10, array, bartlett, linspace
    >>> from numpy.fft import fft, fftshift
    >>> import matplotlib.pyplot as plt

    >>> window = bartlett(51)
    >>> plt.plot(window)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Bartlett window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Sample")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    >>> plt.figure()
    <matplotlib.figure.Figure object at 0x...>
    >>> A = fft(window, 2048) / 25.5
    >>> mag = abs(fftshift(A))
    >>> freq = linspace(-0.5,0.5,len(A))
    >>> response = 20*log10(mag)
    >>> response = clip(response,-100,100)
    >>> plt.plot(freq, response)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Frequency response of Bartlett window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Magnitude [dB]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.axis('tight')
    (-0.5, 0.5, -100.0, ...)
    >>> plt.show()
    
    """
    raise NotImplementedError


def base_repr(number, base=2, padding=0):
    """Return a string representation of a number in the given base system.

    Parameters
    ----------
    number : int
        The value to convert. Only positive values are handled.
    base : int, optional
        Convert `number` to the `base` number system. The valid range is 2-36,
        the default value is 2.
    padding : int, optional
        Number of zeros padded on the left. Default is 0 (no padding).

    Returns
    -------
    out : str
        String representation of `number` in `base` system.

    See Also
    --------
    binary_repr : Faster version of `base_repr` for base 2.

    Examples
    --------
    >>> np.base_repr(5)
    '101'
    >>> np.base_repr(6, 5)
    '11'
    >>> np.base_repr(7, base=5, padding=3)
    '00012'

    >>> np.base_repr(10, base=16)
    'A'
    >>> np.base_repr(32, base=16)
    '20'
    
    """
    raise NotImplementedError


def binary_repr(num, width=None):
    """Return the binary representation of the input number as a string.

    For negative numbers, if width is not given, a minus sign is added to the
    front. If width is given, the two's complement of the number is
    returned, with respect to that width.

    In a two's-complement system negative numbers are represented by the two's
    complement of the absolute value. This is the most common method of
    representing signed integers on computers [1]_. A N-bit two's-complement
    system can represent every integer in the range
    :math:`-2^{N-1}` to :math:`+2^{N-1}-1`.

    Parameters
    ----------
    num : int
        Only an integer decimal number can be used.
    width : int, optional
        The length of the returned string if `num` is positive, the length of
        the two's complement if `num` is negative.

    Returns
    -------
    bin : str
        Binary representation of `num` or two's complement of `num`.

    See Also
    --------
    base_repr: Return a string representation of a number in the given base
               system.

    Notes
    -----
    `binary_repr` is equivalent to using `base_repr` with base 2, but about 25x
    faster.

    References
    ----------
    .. [1] Wikipedia, "Two's complement",
        http://en.wikipedia.org/wiki/Two's_complement

    Examples
    --------
    >>> np.binary_repr(3)
    '11'
    >>> np.binary_repr(-3)
    '-11'
    >>> np.binary_repr(3, width=4)
    '0011'

    The two's complement is returned when the input number is negative and
    width is specified:

    >>> np.binary_repr(-3, width=4)
    '1101'
    
    """
    raise NotImplementedError


def bincount(x, weights=None, minlength=None):
    """bincount(x, weights=None, minlength=None)

    Count number of occurrences of each value in array of non-negative ints.

    The number of bins (of size 1) is one larger than the largest value in
    `x`. If `minlength` is specified, there will be at least this number
    of bins in the output array (though it will be longer if necessary,
    depending on the contents of `x`).
    Each bin gives the number of occurrences of its index value in `x`.
    If `weights` is specified the input array is weighted by it, i.e. if a
    value ``n`` is found at position ``i``, ``out[n] += weight[i]`` instead
    of ``out[n] += 1``.

    Parameters
    ----------
    x : array_like, 1 dimension, nonnegative ints
        Input array.
    weights : array_like, optional
        Weights, array of the same shape as `x`.
    minlength : int, optional
        .. versionadded:: 1.6.0

        A minimum number of bins for the output array.

    Returns
    -------
    out : ndarray of ints
        The result of binning the input array.
        The length of `out` is equal to ``np.amax(x)+1``.

    Raises
    ------
    ValueError
        If the input is not 1-dimensional, or contains elements with negative
        values, or if `minlength` is non-positive.
    TypeError
        If the type of the input is float or complex.

    See Also
    --------
    histogram, digitize, unique

    Examples
    --------
    >>> np.bincount(np.arange(5))
    array([1, 1, 1, 1, 1])
    >>> np.bincount(np.array([0, 1, 1, 3, 2, 1, 7]))
    array([1, 3, 1, 1, 0, 0, 0, 1])

    >>> x = np.array([0, 1, 1, 3, 2, 1, 7, 23])
    >>> np.bincount(x).size == np.amax(x)+1
    True

    The input array needs to be of integer dtype, otherwise a
    TypeError is raised:

    >>> np.bincount(np.arange(5, dtype=np.float))
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: array cannot be safely cast to required type

    A possible use of ``bincount`` is to perform sums over
    variable-size chunks of an array, using the ``weights`` keyword.

    >>> w = np.array([0.3, 0.5, 0.2, 0.7, 1., -0.6]) # weights
    >>> x = np.array([0, 1, 1, 2, 2, 2])
    >>> np.bincount(x,  weights=w)
    array([ 0.3,  0.7,  1.1])
    
    """
    # BUILTIN
    raise NotImplementedError


def blackman(M):
    """Return the Blackman window.

    The Blackman window is a taper formed by using the the first three
    terms of a summation of cosines. It was designed to have close to the
    minimal leakage possible.  It is close to optimal, only slightly worse
    than a Kaiser window.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an empty
        array is returned.

    Returns
    -------
    out : ndarray
        The window, normalized to one (the value one appears only if the
        number of samples is odd).

    See Also
    --------
    bartlett, hamming, hanning, kaiser

    Notes
    -----
    The Blackman window is defined as

    .. math::  w(n) = 0.42 - 0.5 \cos(2\pi n/M) + 0.08 \cos(4\pi n/M)

    Most references to the Blackman window come from the signal processing
    literature, where it is used as one of many windowing functions for
    smoothing values.  It is also known as an apodization (which means
    "removing the foot", i.e. smoothing discontinuities at the beginning
    and end of the sampled signal) or tapering function. It is known as a
    "near optimal" tapering function, almost as good (by some measures)
    as the kaiser window.

    References
    ----------
    Blackman, R.B. and Tukey, J.W., (1958) The measurement of power spectra,
    Dover Publications, New York.

    Oppenheim, A.V., and R.W. Schafer. Discrete-Time Signal Processing.
    Upper Saddle River, NJ: Prentice-Hall, 1999, pp. 468-471.

    Examples
    --------
    >>> from numpy import blackman
    >>> blackman(12)
    array([ -1.38777878e-17,   3.26064346e-02,   1.59903635e-01,
             4.14397981e-01,   7.36045180e-01,   9.67046769e-01,
             9.67046769e-01,   7.36045180e-01,   4.14397981e-01,
             1.59903635e-01,   3.26064346e-02,  -1.38777878e-17])


    Plot the window and the frequency response:

    >>> from numpy import clip, log10, array, blackman, linspace
    >>> from numpy.fft import fft, fftshift
    >>> import matplotlib.pyplot as plt

    >>> window = blackman(51)
    >>> plt.plot(window)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Blackman window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Sample")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    >>> plt.figure()
    <matplotlib.figure.Figure object at 0x...>
    >>> A = fft(window, 2048) / 25.5
    >>> mag = abs(fftshift(A))
    >>> freq = linspace(-0.5,0.5,len(A))
    >>> response = 20*log10(mag)
    >>> response = clip(response,-100,100)
    >>> plt.plot(freq, response)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Frequency response of Blackman window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Magnitude [dB]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.axis('tight')
    (-0.5, 0.5, -100.0, ...)
    >>> plt.show()
    
    """
    raise NotImplementedError


def bmat(obj, ldict=None, gdict=None):
    """Build a matrix object from a string, nested sequence, or array.

    Parameters
    ----------
    obj : str or array_like
        Input data.  Names of variables in the current scope may be
        referenced, even if `obj` is a string.

    Returns
    -------
    out : matrix
        Returns a matrix object, which is a specialized 2-D array.

    See Also
    --------
    matrix

    Examples
    --------
    >>> A = np.mat('1 1; 1 1')
    >>> B = np.mat('2 2; 2 2')
    >>> C = np.mat('3 4; 5 6')
    >>> D = np.mat('7 8; 9 0')

    All the following expressions construct the same block matrix:

    >>> np.bmat([[A, B], [C, D]])
    matrix([[1, 1, 2, 2],
            [1, 1, 2, 2],
            [3, 4, 7, 8],
            [5, 6, 9, 0]])
    >>> np.bmat(np.r_[np.c_[A, B], np.c_[C, D]])
    matrix([[1, 1, 2, 2],
            [1, 1, 2, 2],
            [3, 4, 7, 8],
            [5, 6, 9, 0]])
    >>> np.bmat('A,B; C,D')
    matrix([[1, 1, 2, 2],
            [1, 1, 2, 2],
            [3, 4, 7, 8],
            [5, 6, 9, 0]])
    
    """
    raise NotImplementedError


def broadcast_arrays():
    """Broadcast any number of arrays against each other.

    Parameters
    ----------
    `*args` : array_likes
        The arrays to broadcast.

    Returns
    -------
    broadcasted : list of arrays
        These arrays are views on the original arrays.  They are typically
        not contiguous.  Furthermore, more than one element of a
        broadcasted array may refer to a single memory location.  If you
        need to write to the arrays, make copies first.

    Examples
    --------
    >>> x = np.array([[1,2,3]])
    >>> y = np.array([[1],[2],[3]])
    >>> np.broadcast_arrays(x, y)
    [array([[1, 2, 3],
           [1, 2, 3],
           [1, 2, 3]]), array([[1, 1, 1],
           [2, 2, 2],
           [3, 3, 3]])]

    Here is a useful idiom for getting contiguous copies instead of
    non-contiguous views.

    >>> map(np.array, np.broadcast_arrays(x, y))
    [array([[1, 2, 3],
           [1, 2, 3],
           [1, 2, 3]]), array([[1, 1, 1],
           [2, 2, 2],
           [3, 3, 3]])]
    
    """
    raise NotImplementedError


def byte_bounds(a):
    """Returns pointers to the end-points of an array.

    Parameters
    ----------
    a : ndarray
        Input array. It must conform to the Python-side of the array interface.

    Returns
    -------
    (low, high) : tuple of 2 integers
        The first integer is the first byte of the array, the second integer is
        just past the last byte of the array.  If `a` is not contiguous it
        will not use every byte between the (`low`, `high`) values.

    Examples
    --------
    >>> I = np.eye(2, dtype='f'); I.dtype
    dtype('float32')
    >>> low, high = np.byte_bounds(I)
    >>> high - low == I.size*I.itemsize
    True
    >>> I = np.eye(2, dtype='G'); I.dtype
    dtype('complex192')
    >>> low, high = np.byte_bounds(I)
    >>> high - low == I.size*I.itemsize
    True
    
    """
    raise NotImplementedError


def can_cast():
    """can_cast(from, totype, casting = 'safe')

    Returns True if cast between data types can occur according to the
    casting rule.  If from is a scalar or array scalar, also returns
    True if the scalar value can be cast without overflow or truncation
    to an integer.

    Parameters
    ----------
    from : dtype, dtype specifier, scalar, or array
        Data type, scalar, or array to cast from.
    totype : dtype or dtype specifier
        Data type to cast to.
    casting : {'no', 'equiv', 'safe', 'same_kind', 'unsafe'}, optional
        Controls what kind of data casting may occur.

          * 'no' means the data types should not be cast at all.
          * 'equiv' means only byte-order changes are allowed.
          * 'safe' means only casts which can preserve values are allowed.
          * 'same_kind' means only safe casts or casts within a kind,
            like float64 to float32, are allowed.
          * 'unsafe' means any data conversions may be done.

    Returns
    -------
    out : bool
        True if cast can occur according to the casting rule.

    See also
    --------
    dtype, result_type

    Examples
    --------

    Basic examples

    >>> np.can_cast(np.int32, np.int64)
    True
    >>> np.can_cast(np.float64, np.complex)
    True
    >>> np.can_cast(np.complex, np.float)
    False

    >>> np.can_cast('i8', 'f8')
    True
    >>> np.can_cast('i8', 'f4')
    False
    >>> np.can_cast('i4', 'S4')
    True

    Casting scalars

    >>> np.can_cast(100, 'i1')
    True
    >>> np.can_cast(150, 'i1')
    False
    >>> np.can_cast(150, 'u1')
    True

    >>> np.can_cast(3.5e100, np.float32)
    False
    >>> np.can_cast(1000.0, np.float32)
    True

    Array scalar checks the value, array does not

    >>> np.can_cast(np.array(1000.0), np.float32)
    True
    >>> np.can_cast(np.array([1000.0]), np.float32)
    False

    Using the casting rules

    >>> np.can_cast('i8', 'i8', 'no')
    True
    >>> np.can_cast('<i8', '>i8', 'no')
    False

    >>> np.can_cast('<i8', '>i8', 'equiv')
    True
    >>> np.can_cast('<i4', '>i8', 'equiv')
    False

    >>> np.can_cast('<i4', '>i8', 'safe')
    True
    >>> np.can_cast('<i8', '>i4', 'safe')
    False

    >>> np.can_cast('<i8', '>i4', 'same_kind')
    True
    >>> np.can_cast('<i8', '>u4', 'same_kind')
    False

    >>> np.can_cast('<i8', '>u4', 'unsafe')
    True
    
    """
    # BUILTIN
    raise NotImplementedError


def choose(a, choices, out=None, mode='raise'):
    """Construct an array from an index array and a set of arrays to choose from.

    First of all, if confused or uncertain, definitely look at the Examples -
    in its full generality, this function is less simple than it might
    seem from the following code description (below ndi =
    `numpy.lib.index_tricks`):

    ``np.choose(a,c) == np.array([c[a[I]][I] for I in ndi.ndindex(a.shape)])``.

    But this omits some subtleties.  Here is a fully general summary:

    Given an "index" array (`a`) of integers and a sequence of `n` arrays
    (`choices`), `a` and each choice array are first broadcast, as necessary,
    to arrays of a common shape; calling these *Ba* and *Bchoices[i], i =
    0,...,n-1* we have that, necessarily, ``Ba.shape == Bchoices[i].shape``
    for each `i`.  Then, a new array with shape ``Ba.shape`` is created as
    follows:

    * if ``mode=raise`` (the default), then, first of all, each element of
      `a` (and thus `Ba`) must be in the range `[0, n-1]`; now, suppose that
      `i` (in that range) is the value at the `(j0, j1, ..., jm)` position
      in `Ba` - then the value at the same position in the new array is the
      value in `Bchoices[i]` at that same position;

    * if ``mode=wrap``, values in `a` (and thus `Ba`) may be any (signed)
      integer; modular arithmetic is used to map integers outside the range
      `[0, n-1]` back into that range; and then the new array is constructed
      as above;

    * if ``mode=clip``, values in `a` (and thus `Ba`) may be any (signed)
      integer; negative integers are mapped to 0; values greater than `n-1`
      are mapped to `n-1`; and then the new array is constructed as above.

    Parameters
    ----------
    a : int array
        This array must contain integers in `[0, n-1]`, where `n` is the number
        of choices, unless ``mode=wrap`` or ``mode=clip``, in which cases any
        integers are permissible.
    choices : sequence of arrays
        Choice arrays. `a` and all of the choices must be broadcastable to the
        same shape.  If `choices` is itself an array (not recommended), then
        its outermost dimension (i.e., the one corresponding to
        ``choices.shape[0]``) is taken as defining the "sequence".
    out : array, optional
        If provided, the result will be inserted into this array. It should
        be of the appropriate shape and dtype.
    mode : {'raise' (default), 'wrap', 'clip'}, optional
        Specifies how indices outside `[0, n-1]` will be treated:

          * 'raise' : an exception is raised
          * 'wrap' : value becomes value mod `n`
          * 'clip' : values < 0 are mapped to 0, values > n-1 are mapped to n-1

    Returns
    -------
    merged_array : array
        The merged result.

    Raises
    ------
    ValueError: shape mismatch
        If `a` and each choice array are not all broadcastable to the same
        shape.

    See Also
    --------
    ndarray.choose : equivalent method

    Notes
    -----
    To reduce the chance of misinterpretation, even though the following
    "abuse" is nominally supported, `choices` should neither be, nor be
    thought of as, a single array, i.e., the outermost sequence-like container
    should be either a list or a tuple.

    Examples
    --------

    >>> choices = [[0, 1, 2, 3], [10, 11, 12, 13],
    ...   [20, 21, 22, 23], [30, 31, 32, 33]]
    >>> np.choose([2, 3, 1, 0], choices
    ... # the first element of the result will be the first element of the
    ... # third (2+1) "array" in choices, namely, 20; the second element
    ... # will be the second element of the fourth (3+1) choice array, i.e.,
    ... # 31, etc.
    ... )
    array([20, 31, 12,  3])
    >>> np.choose([2, 4, 1, 0], choices, mode='clip') # 4 goes to 3 (4-1)
    array([20, 31, 12,  3])
    >>> # because there are 4 choice arrays
    >>> np.choose([2, 4, 1, 0], choices, mode='wrap') # 4 goes to (4 mod 4)
    array([20,  1, 12,  3])
    >>> # i.e., 0

    A couple examples illustrating how choose broadcasts:

    >>> a = [[1, 0, 1], [0, 1, 0], [1, 0, 1]]
    >>> choices = [-10, 10]
    >>> np.choose(a, choices)
    array([[ 10, -10,  10],
           [-10,  10, -10],
           [ 10, -10,  10]])

    >>> # With thanks to Anne Archibald
    >>> a = np.array([0, 1]).reshape((2,1,1))
    >>> c1 = np.array([1, 2, 3]).reshape((1,3,1))
    >>> c2 = np.array([-1, -2, -3, -4, -5]).reshape((1,1,5))
    >>> np.choose(a, (c1, c2)) # result is 2x3x5, res[0,:,:]=c1, res[1,:,:]=c2
    array([[[ 1,  1,  1,  1,  1],
            [ 2,  2,  2,  2,  2],
            [ 3,  3,  3,  3,  3]],
           [[-1, -2, -3, -4, -5],
            [-1, -2, -3, -4, -5],
            [-1, -2, -3, -4, -5]]])
    
    """
    raise NotImplementedError


def column_stack(tup):
    """Stack 1-D arrays as columns into a 2-D array.

    Take a sequence of 1-D arrays and stack them as columns
    to make a single 2-D array. 2-D arrays are stacked as-is,
    just like with `hstack`.  1-D arrays are turned into 2-D columns
    first.

    Parameters
    ----------
    tup : sequence of 1-D or 2-D arrays.
        Arrays to stack. All of them must have the same first dimension.

    Returns
    -------
    stacked : 2-D array
        The array formed by stacking the given arrays.

    See Also
    --------
    hstack, vstack, concatenate

    Notes
    -----
    This function is equivalent to ``np.vstack(tup).T``.

    Examples
    --------
    >>> a = np.array((1,2,3))
    >>> b = np.array((2,3,4))
    >>> np.column_stack((a,b))
    array([[1, 2],
           [2, 3],
           [3, 4]])
    
    """
    raise NotImplementedError


def common_type():
    """Return a scalar type which is common to the input arrays.

    The return type will always be an inexact (i.e. floating point) scalar
    type, even if all the arrays are integer arrays. If one of the inputs is
    an integer array, the minimum precision type that is returned is a
    64-bit floating point dtype.

    All input arrays can be safely cast to the returned dtype without loss
    of information.

    Parameters
    ----------
    array1, array2, ... : ndarrays
        Input arrays.

    Returns
    -------
    out : data type code
        Data type code.

    See Also
    --------
    dtype, mintypecode

    Examples
    --------
    >>> np.common_type(np.arange(2, dtype=np.float32))
    <type 'numpy.float32'>
    >>> np.common_type(np.arange(2, dtype=np.float32), np.arange(2))
    <type 'numpy.float64'>
    >>> np.common_type(np.arange(4), np.array([45, 6.j]), np.array([45.0]))
    <type 'numpy.complex128'>
    
    """
    raise NotImplementedError


def compare_chararrays():
    """
    
    """
    # BUILTIN
    raise NotImplementedError


def compress(condition, a, axis=None, out=None):
    """Return selected slices of an array along given axis.

    When working along a given axis, a slice along that axis is returned in
    `output` for each index where `condition` evaluates to True. When
    working on a 1-D array, `compress` is equivalent to `extract`.

    Parameters
    ----------
    condition : 1-D array of bools
        Array that selects which entries to return. If len(condition)
        is less than the size of `a` along the given axis, then output is
        truncated to the length of the condition array.
    a : array_like
        Array from which to extract a part.
    axis : int, optional
        Axis along which to take slices. If None (default), work on the
        flattened array.
    out : ndarray, optional
        Output array.  Its type is preserved and it must be of the right
        shape to hold the output.

    Returns
    -------
    compressed_array : ndarray
        A copy of `a` without the slices along axis for which `condition`
        is false.

    See Also
    --------
    take, choose, diag, diagonal, select
    ndarray.compress : Equivalent method.
    numpy.doc.ufuncs : Section "Output arguments"

    Examples
    --------
    >>> a = np.array([[1, 2], [3, 4], [5, 6]])
    >>> a
    array([[1, 2],
           [3, 4],
           [5, 6]])
    >>> np.compress([0, 1], a, axis=0)
    array([[3, 4]])
    >>> np.compress([False, True, True], a, axis=0)
    array([[3, 4],
           [5, 6]])
    >>> np.compress([False, True], a, axis=1)
    array([[2],
           [4],
           [6]])

    Working on the flattened array does not return slices along an axis but
    selects elements.

    >>> np.compress([False, True], a)
    array([2])
    
    """
    raise NotImplementedError


def concatenate():
    """concatenate((a1, a2, ...), axis=0)

    Join a sequence of arrays together.

    Parameters
    ----------
    a1, a2, ... : sequence of array_like
        The arrays must have the same shape, except in the dimension
        corresponding to `axis` (the first, by default).
    axis : int, optional
        The axis along which the arrays will be joined.  Default is 0.

    Returns
    -------
    res : ndarray
        The concatenated array.

    See Also
    --------
    ma.concatenate : Concatenate function that preserves input masks.
    array_split : Split an array into multiple sub-arrays of equal or
                  near-equal size.
    split : Split array into a list of multiple sub-arrays of equal size.
    hsplit : Split array into multiple sub-arrays horizontally (column wise)
    vsplit : Split array into multiple sub-arrays vertically (row wise)
    dsplit : Split array into multiple sub-arrays along the 3rd axis (depth).
    hstack : Stack arrays in sequence horizontally (column wise)
    vstack : Stack arrays in sequence vertically (row wise)
    dstack : Stack arrays in sequence depth wise (along third dimension)

    Notes
    -----
    When one or more of the arrays to be concatenated is a MaskedArray,
    this function will return a MaskedArray object instead of an ndarray,
    but the input masks are *not* preserved. In cases where a MaskedArray
    is expected as input, use the ma.concatenate function from the masked
    array module instead.

    Examples
    --------
    >>> a = np.array([[1, 2], [3, 4]])
    >>> b = np.array([[5, 6]])
    >>> np.concatenate((a, b), axis=0)
    array([[1, 2],
           [3, 4],
           [5, 6]])
    >>> np.concatenate((a, b.T), axis=1)
    array([[1, 2, 5],
           [3, 4, 6]])

    This function will not preserve masking of MaskedArray inputs.

    >>> a = np.ma.arange(3)
    >>> a[1] = np.ma.masked
    >>> b = np.arange(2, 5)
    >>> a
    masked_array(data = [0 -- 2],
                 mask = [False  True False],
           fill_value = 999999)
    >>> b
    array([2, 3, 4])
    >>> np.concatenate([a, b])
    masked_array(data = [0 1 2 2 3 4],
                 mask = False,
           fill_value = 999999)
    >>> np.ma.concatenate([a, b])
    masked_array(data = [0 -- 2 2 3 4],
                 mask = [False  True False False False False],
           fill_value = 999999)
    
    """
    # BUILTIN
    raise NotImplementedError


def convolve(a, v, mode='full'):
    """Returns the discrete, linear convolution of two one-dimensional sequences.

    The convolution operator is often seen in signal processing, where it
    models the effect of a linear time-invariant system on a signal [1]_.  In
    probability theory, the sum of two independent random variables is
    distributed according to the convolution of their individual
    distributions.

    Parameters
    ----------
    a : (N,) array_like
        First one-dimensional input array.
    v : (M,) array_like
        Second one-dimensional input array.
    mode : {'full', 'valid', 'same'}, optional
        'full':
          By default, mode is 'full'.  This returns the convolution
          at each point of overlap, with an output shape of (N+M-1,). At
          the end-points of the convolution, the signals do not overlap
          completely, and boundary effects may be seen.

        'same':
          Mode `same` returns output of length ``max(M, N)``.  Boundary
          effects are still visible.

        'valid':
          Mode `valid` returns output of length
          ``max(M, N) - min(M, N) + 1``.  The convolution product is only given
          for points where the signals overlap completely.  Values outside
          the signal boundary have no effect.

    Returns
    -------
    out : ndarray
        Discrete, linear convolution of `a` and `v`.

    See Also
    --------
    scipy.signal.fftconvolve : Convolve two arrays using the Fast Fourier
                               Transform.
    scipy.linalg.toeplitz : Used to construct the convolution operator.

    Notes
    -----
    The discrete convolution operation is defined as

    .. math:: (f * g)[n] = \sum_{m = -\infty}^{\infty} f[m] g[n - m]

    It can be shown that a convolution :math:`x(t) * y(t)` in time/space
    is equivalent to the multiplication :math:`X(f) Y(f)` in the Fourier
    domain, after appropriate padding (padding is necessary to prevent
    circular convolution).  Since multiplication is more efficient (faster)
    than convolution, the function `scipy.signal.fftconvolve` exploits the
    FFT to calculate the convolution of large data-sets.

    References
    ----------
    .. [1] Wikipedia, "Convolution", http://en.wikipedia.org/wiki/Convolution.

    Examples
    --------
    Note how the convolution operator flips the second array
    before "sliding" the two across one another:

    >>> np.convolve([1, 2, 3], [0, 1, 0.5])
    array([ 0. ,  1. ,  2.5,  4. ,  1.5])

    Only return the middle values of the convolution.
    Contains boundary effects, where zeros are taken
    into account:

    >>> np.convolve([1,2,3],[0,1,0.5], 'same')
    array([ 1. ,  2.5,  4. ])

    The two arrays are of the same length, so there
    is only one position where they completely overlap:

    >>> np.convolve([1,2,3],[0,1,0.5], 'valid')
    array([ 2.5])
    
    """
    raise NotImplementedError


def copy(a):
    """Return an array copy of the given object.

    Parameters
    ----------
    a : array_like
        Input data.

    Returns
    -------
    arr : ndarray
        Array interpretation of `a`.

    Notes
    -----
    This is equivalent to

    >>> np.array(a, copy=True)                              #doctest: +SKIP

    Examples
    --------
    Create an array x, with a reference y and a copy z:

    >>> x = np.array([1, 2, 3])
    >>> y = x
    >>> z = np.copy(x)

    Note that, when we modify x, y changes, but not z:

    >>> x[0] = 10
    >>> x[0] == y[0]
    True
    >>> x[0] == z[0]
    False
    
    """
    raise NotImplementedError


def corrcoef(x, y=None, rowvar=1, bias=0, ddof=None):
    """Return correlation coefficients.

    Please refer to the documentation for `cov` for more detail.  The
    relationship between the correlation coefficient matrix, `P`, and the
    covariance matrix, `C`, is

    .. math:: P_{ij} = \frac{ C_{ij} } { \sqrt{ C_{ii} * C_{jj} } }

    The values of `P` are between -1 and 1, inclusive.

    Parameters
    ----------
    x : array_like
        A 1-D or 2-D array containing multiple variables and observations.
        Each row of `m` represents a variable, and each column a single
        observation of all those variables. Also see `rowvar` below.
    y : array_like, optional
        An additional set of variables and observations. `y` has the same
        shape as `m`.
    rowvar : int, optional
        If `rowvar` is non-zero (default), then each row represents a
        variable, with observations in the columns. Otherwise, the relationship
        is transposed: each column represents a variable, while the rows
        contain observations.
    bias : int, optional
        Default normalization is by ``(N - 1)``, where ``N`` is the number of
        observations (unbiased estimate). If `bias` is 1, then
        normalization is by ``N``. These values can be overridden by using
        the keyword ``ddof`` in numpy versions >= 1.5.
    ddof : {None, int}, optional
        .. versionadded:: 1.5
        If not ``None`` normalization is by ``(N - ddof)``, where ``N`` is
        the number of observations; this overrides the value implied by
        ``bias``. The default value is ``None``.

    Returns
    -------
    out : ndarray
        The correlation coefficient matrix of the variables.

    See Also
    --------
    cov : Covariance matrix
    
    """
    raise NotImplementedError


def correlate(a, v, mode='valid', old_behavior=False):
    """Cross-correlation of two 1-dimensional sequences.

    This function computes the correlation as generally defined in signal
    processing texts::

        z[k] = sum_n a[n] * conj(v[n+k])

    with a and v sequences being zero-padded where necessary and conj being
    the conjugate.

    Parameters
    ----------
    a, v : array_like
        Input sequences.
    mode : {'valid', 'same', 'full'}, optional
        Refer to the `convolve` docstring.  Note that the default
        is `valid`, unlike `convolve`, which uses `full`.
    old_behavior : bool
        If True, uses the old behavior from Numeric, (correlate(a,v) == correlate(v,
        a), and the conjugate is not taken for complex arrays). If False, uses
        the conventional signal processing definition (see note).

    See Also
    --------
    convolve : Discrete, linear convolution of two one-dimensional sequences.

    Examples
    --------
    >>> np.correlate([1, 2, 3], [0, 1, 0.5])
    array([ 3.5])
    >>> np.correlate([1, 2, 3], [0, 1, 0.5], "same")
    array([ 2. ,  3.5,  3. ])
    >>> np.correlate([1, 2, 3], [0, 1, 0.5], "full")
    array([ 0.5,  2. ,  3.5,  3. ,  0. ])
    
    """
    raise NotImplementedError


def count_nonzero(a):
    """count_nonzero(a)

    Counts the number of non-zero values in the array ``a``.

    Parameters
    ----------
    a : array_like
        The array for which to count non-zeros.

    Returns
    -------
    count : int
        Number of non-zero values in the array.

    See Also
    --------
    nonzero : Return the coordinates of all the non-zero values.

    Examples
    --------
    >>> np.count_nonzero(np.eye(4))
    4

    >>> np.count_nonzero([[0,1,7,0,0],[3,0,0,2,19]])
    5
    
    """
    # BUILTIN
    raise NotImplementedError


def cov(m, y=None, rowvar=1, bias=0, ddof=None):
    """Estimate a covariance matrix, given data.

    Covariance indicates the level to which two variables vary together.
    If we examine N-dimensional samples, :math:`X = [x_1, x_2, ... x_N]^T`,
    then the covariance matrix element :math:`C_{ij}` is the covariance of
    :math:`x_i` and :math:`x_j`. The element :math:`C_{ii}` is the variance
    of :math:`x_i`.

    Parameters
    ----------
    m : array_like
        A 1-D or 2-D array containing multiple variables and observations.
        Each row of `m` represents a variable, and each column a single
        observation of all those variables. Also see `rowvar` below.
    y : array_like, optional
        An additional set of variables and observations. `y` has the same
        form as that of `m`.
    rowvar : int, optional
        If `rowvar` is non-zero (default), then each row represents a
        variable, with observations in the columns. Otherwise, the relationship
        is transposed: each column represents a variable, while the rows
        contain observations.
    bias : int, optional
        Default normalization is by ``(N - 1)``, where ``N`` is the number of
        observations given (unbiased estimate). If `bias` is 1, then
        normalization is by ``N``. These values can be overridden by using
        the keyword ``ddof`` in numpy versions >= 1.5.
    ddof : int, optional
        .. versionadded:: 1.5
        If not ``None`` normalization is by ``(N - ddof)``, where ``N`` is
        the number of observations; this overrides the value implied by
        ``bias``. The default value is ``None``.

    Returns
    -------
    out : ndarray
        The covariance matrix of the variables.

    See Also
    --------
    corrcoef : Normalized covariance matrix

    Examples
    --------
    Consider two variables, :math:`x_0` and :math:`x_1`, which
    correlate perfectly, but in opposite directions:

    >>> x = np.array([[0, 2], [1, 1], [2, 0]]).T
    >>> x
    array([[0, 1, 2],
           [2, 1, 0]])

    Note how :math:`x_0` increases while :math:`x_1` decreases. The covariance
    matrix shows this clearly:

    >>> np.cov(x)
    array([[ 1., -1.],
           [-1.,  1.]])

    Note that element :math:`C_{0,1}`, which shows the correlation between
    :math:`x_0` and :math:`x_1`, is negative.

    Further, note how `x` and `y` are combined:

    >>> x = [-2.1, -1,  4.3]
    >>> y = [3,  1.1,  0.12]
    >>> X = np.vstack((x,y))
    >>> print np.cov(X)
    [[ 11.71        -4.286     ]
     [ -4.286        2.14413333]]
    >>> print np.cov(x, y)
    [[ 11.71        -4.286     ]
     [ -4.286        2.14413333]]
    >>> print np.cov(x)
    11.71
    
    """
    raise NotImplementedError


def cross(a, b, axisa=-1, axisb=-1, axisc=-1, axis=None):
    """Return the cross product of two (arrays of) vectors.

    The cross product of `a` and `b` in :math:`R^3` is a vector perpendicular
    to both `a` and `b`.  If `a` and `b` are arrays of vectors, the vectors
    are defined by the last axis of `a` and `b` by default, and these axes
    can have dimensions 2 or 3.  Where the dimension of either `a` or `b` is
    2, the third component of the input vector is assumed to be zero and the
    cross product calculated accordingly.  In cases where both input vectors
    have dimension 2, the z-component of the cross product is returned.

    Parameters
    ----------
    a : array_like
        Components of the first vector(s).
    b : array_like
        Components of the second vector(s).
    axisa : int, optional
        Axis of `a` that defines the vector(s).  By default, the last axis.
    axisb : int, optional
        Axis of `b` that defines the vector(s).  By default, the last axis.
    axisc : int, optional
        Axis of `c` containing the cross product vector(s).  By default, the
        last axis.
    axis : int, optional
        If defined, the axis of `a`, `b` and `c` that defines the vector(s)
        and cross product(s).  Overrides `axisa`, `axisb` and `axisc`.

    Returns
    -------
    c : ndarray
        Vector cross product(s).

    Raises
    ------
    ValueError
        When the dimension of the vector(s) in `a` and/or `b` does not
        equal 2 or 3.

    See Also
    --------
    inner : Inner product
    outer : Outer product.
    ix_ : Construct index arrays.

    Examples
    --------
    Vector cross-product.

    >>> x = [1, 2, 3]
    >>> y = [4, 5, 6]
    >>> np.cross(x, y)
    array([-3,  6, -3])

    One vector with dimension 2.

    >>> x = [1, 2]
    >>> y = [4, 5, 6]
    >>> np.cross(x, y)
    array([12, -6, -3])

    Equivalently:

    >>> x = [1, 2, 0]
    >>> y = [4, 5, 6]
    >>> np.cross(x, y)
    array([12, -6, -3])

    Both vectors with dimension 2.

    >>> x = [1,2]
    >>> y = [4,5]
    >>> np.cross(x, y)
    -3

    Multiple vector cross-products. Note that the direction of the cross
    product vector is defined by the `right-hand rule`.

    >>> x = np.array([[1,2,3], [4,5,6]])
    >>> y = np.array([[4,5,6], [1,2,3]])
    >>> np.cross(x, y)
    array([[-3,  6, -3],
           [ 3, -6,  3]])

    The orientation of `c` can be changed using the `axisc` keyword.

    >>> np.cross(x, y, axisc=0)
    array([[-3,  3],
           [ 6, -6],
           [-3,  3]])

    Change the vector definition of `x` and `y` using `axisa` and `axisb`.

    >>> x = np.array([[1,2,3], [4,5,6], [7, 8, 9]])
    >>> y = np.array([[7, 8, 9], [4,5,6], [1,2,3]])
    >>> np.cross(x, y)
    array([[ -6,  12,  -6],
           [  0,   0,   0],
           [  6, -12,   6]])
    >>> np.cross(x, y, axisa=0, axisb=0)
    array([[-24,  48, -24],
           [-30,  60, -30],
           [-36,  72, -36]])
    
    """
    raise NotImplementedError


def cumprod(a, axis=None, dtype=None, out=None):
    """Return the cumulative product of elements along a given axis.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : int, optional
        Axis along which the cumulative product is computed.  By default
        the input is flattened.
    dtype : dtype, optional
        Type of the returned array, as well as of the accumulator in which
        the elements are multiplied.  If *dtype* is not specified, it
        defaults to the dtype of `a`, unless `a` has an integer dtype with
        a precision less than that of the default platform integer.  In
        that case, the default platform integer is used instead.
    out : ndarray, optional
        Alternative output array in which to place the result. It must
        have the same shape and buffer length as the expected output
        but the type of the resulting values will be cast if necessary.

    Returns
    -------
    cumprod : ndarray
        A new array holding the result is returned unless `out` is
        specified, in which case a reference to out is returned.

    See Also
    --------
    numpy.doc.ufuncs : Section "Output arguments"

    Notes
    -----
    Arithmetic is modular when using integer types, and no error is
    raised on overflow.

    Examples
    --------
    >>> a = np.array([1,2,3])
    >>> np.cumprod(a) # intermediate results 1, 1*2
    ...               # total product 1*2*3 = 6
    array([1, 2, 6])
    >>> a = np.array([[1, 2, 3], [4, 5, 6]])
    >>> np.cumprod(a, dtype=float) # specify type of output
    array([   1.,    2.,    6.,   24.,  120.,  720.])

    The cumulative product for each column (i.e., over the rows) of `a`:

    >>> np.cumprod(a, axis=0)
    array([[ 1,  2,  3],
           [ 4, 10, 18]])

    The cumulative product for each row (i.e. over the columns) of `a`:

    >>> np.cumprod(a,axis=1)
    array([[  1,   2,   6],
           [  4,  20, 120]])
    
    """
    raise NotImplementedError


def cumproduct(a, axis=None, dtype=None, out=None):
    """Return the cumulative product over the given axis.


    See Also
    --------
    cumprod : equivalent function; see for details.
    
    """
    raise NotImplementedError


def cumsum(a, axis=None, dtype=None, out=None):
    """Return the cumulative sum of the elements along a given axis.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : int, optional
        Axis along which the cumulative sum is computed. The default
        (None) is to compute the cumsum over the flattened array.
    dtype : dtype, optional
        Type of the returned array and of the accumulator in which the
        elements are summed.  If `dtype` is not specified, it defaults
        to the dtype of `a`, unless `a` has an integer dtype with a
        precision less than that of the default platform integer.  In
        that case, the default platform integer is used.
    out : ndarray, optional
        Alternative output array in which to place the result. It must
        have the same shape and buffer length as the expected output
        but the type will be cast if necessary. See `doc.ufuncs`
        (Section "Output arguments") for more details.

    Returns
    -------
    cumsum_along_axis : ndarray.
        A new array holding the result is returned unless `out` is
        specified, in which case a reference to `out` is returned. The
        result has the same size as `a`, and the same shape as `a` if
        `axis` is not None or `a` is a 1-d array.


    See Also
    --------
    sum : Sum array elements.

    trapz : Integration of array values using the composite trapezoidal rule.

    Notes
    -----
    Arithmetic is modular when using integer types, and no error is
    raised on overflow.

    Examples
    --------
    >>> a = np.array([[1,2,3], [4,5,6]])
    >>> a
    array([[1, 2, 3],
           [4, 5, 6]])
    >>> np.cumsum(a)
    array([ 1,  3,  6, 10, 15, 21])
    >>> np.cumsum(a, dtype=float)     # specifies type of output value(s)
    array([  1.,   3.,   6.,  10.,  15.,  21.])

    >>> np.cumsum(a,axis=0)      # sum over rows for each of the 3 columns
    array([[1, 2, 3],
           [5, 7, 9]])
    >>> np.cumsum(a,axis=1)      # sum over columns for each of the 2 rows
    array([[ 1,  3,  6],
           [ 4,  9, 15]])
    
    """
    raise NotImplementedError


def datetime_data(dtype):
    """Return (unit, numerator, denominator, events) from a datetime dtype
    
    """
    raise NotImplementedError


def delete(arr, obj, axis=None):
    """Return a new array with sub-arrays along an axis deleted.

    Parameters
    ----------
    arr : array_like
      Input array.
    obj : slice, int or array of ints
      Indicate which sub-arrays to remove.
    axis : int, optional
      The axis along which to delete the subarray defined by `obj`.
      If `axis` is None, `obj` is applied to the flattened array.

    Returns
    -------
    out : ndarray
        A copy of `arr` with the elements specified by `obj` removed. Note
        that `delete` does not occur in-place. If `axis` is None, `out` is
        a flattened array.

    See Also
    --------
    insert : Insert elements into an array.
    append : Append elements at the end of an array.

    Examples
    --------
    >>> arr = np.array([[1,2,3,4], [5,6,7,8], [9,10,11,12]])
    >>> arr
    array([[ 1,  2,  3,  4],
           [ 5,  6,  7,  8],
           [ 9, 10, 11, 12]])
    >>> np.delete(arr, 1, 0)
    array([[ 1,  2,  3,  4],
           [ 9, 10, 11, 12]])

    >>> np.delete(arr, np.s_[::2], 1)
    array([[ 2,  4],
           [ 6,  8],
           [10, 12]])
    >>> np.delete(arr, [1,3,5], None)
    array([ 1,  3,  5,  7,  8,  9, 10, 11, 12])
    
    """
    raise NotImplementedError


def deprecate():
    """Issues a DeprecationWarning, adds warning to `old_name`'s
    docstring, rebinds ``old_name.__name__`` and returns the new
    function object.

    This function may also be used as a decorator.

    Parameters
    ----------
    func : function
        The function to be deprecated.
    old_name : str, optional
        The name of the function to be deprecated. Default is None, in which
        case the name of `func` is used.
    new_name : str, optional
        The new name for the function. Default is None, in which case
        the deprecation message is that `old_name` is deprecated. If given,
        the deprecation message is that `old_name` is deprecated and `new_name`
        should be used instead.
    message : str, optional
        Additional explanation of the deprecation.  Displayed in the docstring
        after the warning.

    Returns
    -------
    old_func : function
        The deprecated function.

    Examples
    --------
    Note that ``olduint`` returns a value after printing Deprecation Warning:

    >>> olduint = np.deprecate(np.uint)
    >>> olduint(6)
    /usr/lib/python2.5/site-packages/numpy/lib/utils.py:114:
    DeprecationWarning: uint32 is deprecated
      warnings.warn(str1, DeprecationWarning)
    6
    
    """
    raise NotImplementedError


def deprecate_with_doc(msg):
    """message
    
    """
    raise NotImplementedError


def diag_indices(n, ndim=2):
    """Return the indices to access the main diagonal of an array.

    This returns a tuple of indices that can be used to access the main
    diagonal of an array `a` with ``a.ndim >= 2`` dimensions and shape
    (n, n, ..., n). For ``a.ndim = 2`` this is the usual diagonal, for
    ``a.ndim > 2`` this is the set of indices to access ``a[i, i, ..., i]``
    for ``i = [0..n-1]``.

    Parameters
    ----------
    n : int
      The size, along each dimension, of the arrays for which the returned
      indices can be used.

    ndim : int, optional
      The number of dimensions.

    See also
    --------
    diag_indices_from

    Notes
    -----
    .. versionadded:: 1.4.0

    Examples
    --------
    Create a set of indices to access the diagonal of a (4, 4) array:

    >>> di = np.diag_indices(4)
    >>> di
    (array([0, 1, 2, 3]), array([0, 1, 2, 3]))
    >>> a = np.arange(16).reshape(4, 4)
    >>> a
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11],
           [12, 13, 14, 15]])
    >>> a[di] = 100
    >>> a
    array([[100,   1,   2,   3],
           [  4, 100,   6,   7],
           [  8,   9, 100,  11],
           [ 12,  13,  14, 100]])

    Now, we create indices to manipulate a 3-D array:

    >>> d3 = np.diag_indices(2, 3)
    >>> d3
    (array([0, 1]), array([0, 1]), array([0, 1]))

    And use it to set the diagonal of an array of zeros to 1:

    >>> a = np.zeros((2, 2, 2), dtype=np.int)
    >>> a[d3] = 1
    >>> a
    array([[[1, 0],
            [0, 0]],
           [[0, 0],
            [0, 1]]])
    
    """
    raise NotImplementedError


def diag_indices_from(arr):
    """Return the indices to access the main diagonal of an n-dimensional array.

    See `diag_indices` for full details.

    Parameters
    ----------
    arr : array, at least 2-D

    See Also
    --------
    diag_indices

    Notes
    -----
    .. versionadded:: 1.4.0
    
    """
    raise NotImplementedError


def diagflat(v, k=0):
    """Create a two-dimensional array with the flattened input as a diagonal.

    Parameters
    ----------
    v : array_like
        Input data, which is flattened and set as the `k`-th
        diagonal of the output.
    k : int, optional
        Diagonal to set; 0, the default, corresponds to the "main" diagonal,
        a positive (negative) `k` giving the number of the diagonal above
        (below) the main.

    Returns
    -------
    out : ndarray
        The 2-D output array.

    See Also
    --------
    diag : MATLAB work-alike for 1-D and 2-D arrays.
    diagonal : Return specified diagonals.
    trace : Sum along diagonals.

    Examples
    --------
    >>> np.diagflat([[1,2], [3,4]])
    array([[1, 0, 0, 0],
           [0, 2, 0, 0],
           [0, 0, 3, 0],
           [0, 0, 0, 4]])

    >>> np.diagflat([1,2], 1)
    array([[0, 1, 0],
           [0, 0, 2],
           [0, 0, 0]])
    
    """
    raise NotImplementedError


def diff(a, n=1, axis=-1):
    """Calculate the n-th order discrete difference along given axis.

    The first order difference is given by ``out[n] = a[n+1] - a[n]`` along
    the given axis, higher order differences are calculated by using `diff`
    recursively.

    Parameters
    ----------
    a : array_like
        Input array
    n : int, optional
        The number of times values are differenced.
    axis : int, optional
        The axis along which the difference is taken, default is the last axis.

    Returns
    -------
    out : ndarray
        The `n` order differences. The shape of the output is the same as `a`
        except along `axis` where the dimension is smaller by `n`.

    See Also
    --------
    gradient, ediff1d

    Examples
    --------
    >>> x = np.array([1, 2, 4, 7, 0])
    >>> np.diff(x)
    array([ 1,  2,  3, -7])
    >>> np.diff(x, n=2)
    array([  1,   1, -10])

    >>> x = np.array([[1, 3, 6, 10], [0, 5, 6, 8]])
    >>> np.diff(x)
    array([[2, 3, 4],
           [5, 1, 2]])
    >>> np.diff(x, axis=0)
    array([[-1,  2,  0, -2]])
    
    """
    raise NotImplementedError


def digitize(x, bins):
    """digitize(x, bins)

    Return the indices of the bins to which each value in input array belongs.

    Each index ``i`` returned is such that ``bins[i-1] <= x < bins[i]`` if
    `bins` is monotonically increasing, or ``bins[i-1] > x >= bins[i]`` if
    `bins` is monotonically decreasing. If values in `x` are beyond the
    bounds of `bins`, 0 or ``len(bins)`` is returned as appropriate.

    Parameters
    ----------
    x : array_like
        Input array to be binned. It has to be 1-dimensional.
    bins : array_like
        Array of bins. It has to be 1-dimensional and monotonic.

    Returns
    -------
    out : ndarray of ints
        Output array of indices, of same shape as `x`.

    Raises
    ------
    ValueError
        If the input is not 1-dimensional, or if `bins` is not monotonic.
    TypeError
        If the type of the input is complex.

    See Also
    --------
    bincount, histogram, unique

    Notes
    -----
    If values in `x` are such that they fall outside the bin range,
    attempting to index `bins` with the indices that `digitize` returns
    will result in an IndexError.

    Examples
    --------
    >>> x = np.array([0.2, 6.4, 3.0, 1.6])
    >>> bins = np.array([0.0, 1.0, 2.5, 4.0, 10.0])
    >>> inds = np.digitize(x, bins)
    >>> inds
    array([1, 4, 3, 2])
    >>> for n in range(x.size):
    ...   print bins[inds[n]-1], "<=", x[n], "<", bins[inds[n]]
    ...
    0.0 <= 0.2 < 1.0
    4.0 <= 6.4 < 10.0
    2.5 <= 3.0 < 4.0
    1.0 <= 1.6 < 2.5
    
    """
    # BUILTIN
    raise NotImplementedError


def disp(mesg, device=None, linefeed=True):
    """Display a message on a device.

    Parameters
    ----------
    mesg : str
        Message to display.
    device : object
        Device to write message. If None, defaults to ``sys.stdout`` which is
        very similar to ``print``. `device` needs to have ``write()`` and
        ``flush()`` methods.
    linefeed : bool, optional
        Option whether to print a line feed or not. Defaults to True.

    Raises
    ------
    AttributeError
        If `device` does not have a ``write()`` or ``flush()`` method.

    Examples
    --------
    Besides ``sys.stdout``, a file-like object can also be used as it has
    both required methods:

    >>> from StringIO import StringIO
    >>> buf = StringIO()
    >>> np.disp('"Display" in a file', device=buf)
    >>> buf.getvalue()
    '"Display" in a file\n'
    
    """
    raise NotImplementedError


def dsplit(ary, indices_or_sections):
    """Split array into multiple sub-arrays along the 3rd axis (depth).

    Please refer to the `split` documentation.  `dsplit` is equivalent
    to `split` with ``axis=2``, the array is always split along the third
    axis provided the array dimension is greater than or equal to 3.

    See Also
    --------
    split : Split an array into multiple sub-arrays of equal size.

    Examples
    --------
    >>> x = np.arange(16.0).reshape(2, 2, 4)
    >>> x
    array([[[  0.,   1.,   2.,   3.],
            [  4.,   5.,   6.,   7.]],
           [[  8.,   9.,  10.,  11.],
            [ 12.,  13.,  14.,  15.]]])
    >>> np.dsplit(x, 2)
    [array([[[  0.,   1.],
            [  4.,   5.]],
           [[  8.,   9.],
            [ 12.,  13.]]]),
     array([[[  2.,   3.],
            [  6.,   7.]],
           [[ 10.,  11.],
            [ 14.,  15.]]])]
    >>> np.dsplit(x, np.array([3, 6]))
    [array([[[  0.,   1.,   2.],
            [  4.,   5.,   6.]],
           [[  8.,   9.,  10.],
            [ 12.,  13.,  14.]]]),
     array([[[  3.],
            [  7.]],
           [[ 11.],
            [ 15.]]]),
     array([], dtype=float64)]
    
    """
    raise NotImplementedError


def dstack(tup):
    """Stack arrays in sequence depth wise (along third axis).

    Takes a sequence of arrays and stack them along the third axis
    to make a single array. Rebuilds arrays divided by `dsplit`.
    This is a simple way to stack 2D arrays (images) into a single
    3D array for processing.

    Parameters
    ----------
    tup : sequence of arrays
        Arrays to stack. All of them must have the same shape along all
        but the third axis.

    Returns
    -------
    stacked : ndarray
        The array formed by stacking the given arrays.

    See Also
    --------
    vstack : Stack along first axis.
    hstack : Stack along second axis.
    concatenate : Join arrays.
    dsplit : Split array along third axis.

    Notes
    -----
    Equivalent to ``np.concatenate(tup, axis=2)``.

    Examples
    --------
    >>> a = np.array((1,2,3))
    >>> b = np.array((2,3,4))
    >>> np.dstack((a,b))
    array([[[1, 2],
            [2, 3],
            [3, 4]]])

    >>> a = np.array([[1],[2],[3]])
    >>> b = np.array([[2],[3],[4]])
    >>> np.dstack((a,b))
    array([[[1, 2]],
           [[2, 3]],
           [[3, 4]]])
    
    """
    raise NotImplementedError


def ediff1d(ary, to_end=None, to_begin=None):
    """The differences between consecutive elements of an array.

    Parameters
    ----------
    ary : array_like
        If necessary, will be flattened before the differences are taken.
    to_end : array_like, optional
        Number(s) to append at the end of the returned differences.
    to_begin : array_like, optional
        Number(s) to prepend at the beginning of the returned differences.

    Returns
    -------
    ed : ndarray
        The differences. Loosely, this is ``ary.flat[1:] - ary.flat[:-1]``.

    See Also
    --------
    diff, gradient

    Notes
    -----
    When applied to masked arrays, this function drops the mask information
    if the `to_begin` and/or `to_end` parameters are used.

    Examples
    --------
    >>> x = np.array([1, 2, 4, 7, 0])
    >>> np.ediff1d(x)
    array([ 1,  2,  3, -7])

    >>> np.ediff1d(x, to_begin=-99, to_end=np.array([88, 99]))
    array([-99,   1,   2,   3,  -7,  88,  99])

    The returned array is always 1D.

    >>> y = [[1, 2, 4], [1, 6, 24]]
    >>> np.ediff1d(y)
    array([ 1,  2, -3,  5, 18])
    
    """
    raise NotImplementedError


def einsum():
    """einsum(subscripts, *operands, out=None, dtype=None, order='K', casting='safe')

    Evaluates the Einstein summation convention on the operands.

    Using the Einstein summation convention, many common multi-dimensional
    array operations can be represented in a simple fashion.  This function
    provides a way compute such summations. The best way to understand this
    function is to try the examples below, which show how many common NumPy
    functions can be implemented as calls to `einsum`.

    Parameters
    ----------
    subscripts : str
        Specifies the subscripts for summation.
    operands : list of array_like
        These are the arrays for the operation.
    out : ndarray, optional
        If provided, the calculation is done into this array.
    dtype : data-type, optional
        If provided, forces the calculation to use the data type specified.
        Note that you may have to also give a more liberal `casting`
        parameter to allow the conversions.
    order : {'C', 'F', 'A', or 'K'}, optional
        Controls the memory layout of the output. 'C' means it should
        be C contiguous. 'F' means it should be Fortran contiguous,
        'A' means it should be 'F' if the inputs are all 'F', 'C' otherwise.
        'K' means it should be as close to the layout as the inputs as
        is possible, including arbitrarily permuted axes.
        Default is 'K'.
    casting : {'no', 'equiv', 'safe', 'same_kind', 'unsafe'}, optional
        Controls what kind of data casting may occur.  Setting this to
        'unsafe' is not recommended, as it can adversely affect accumulations.

          * 'no' means the data types should not be cast at all.
          * 'equiv' means only byte-order changes are allowed.
          * 'safe' means only casts which can preserve values are allowed.
          * 'same_kind' means only safe casts or casts within a kind,
            like float64 to float32, are allowed.
          * 'unsafe' means any data conversions may be done.

    Returns
    -------
    output : ndarray
        The calculation based on the Einstein summation convention.

    See Also
    --------
    dot, inner, outer, tensordot

    Notes
    -----
    .. versionadded:: 1.6.0

    The subscripts string is a comma-separated list of subscript labels,
    where each label refers to a dimension of the corresponding operand.
    Repeated subscripts labels in one operand take the diagonal.  For example,
    ``np.einsum('ii', a)`` is equivalent to ``np.trace(a)``.

    Whenever a label is repeated, it is summed, so ``np.einsum('i,i', a, b)``
    is equivalent to ``np.inner(a,b)``.  If a label appears only once,
    it is not summed, so ``np.einsum('i', a)`` produces a view of ``a``
    with no changes.

    The order of labels in the output is by default alphabetical.  This
    means that ``np.einsum('ij', a)`` doesn't affect a 2D array, while
    ``np.einsum('ji', a)`` takes its transpose.

    The output can be controlled by specifying output subscript labels
    as well.  This specifies the label order, and allows summing to
    be disallowed or forced when desired.  The call ``np.einsum('i->', a)``
    is like ``np.sum(a, axis=-1)``, and ``np.einsum('ii->i', a)``
    is like ``np.diag(a)``.  The difference is that `einsum` does not
    allow broadcasting by default.

    To enable and control broadcasting, use an ellipsis.  Default
    NumPy-style broadcasting is done by adding an ellipsis
    to the left of each term, like ``np.einsum('...ii->...i', a)``.
    To take the trace along the first and last axes,
    you can do ``np.einsum('i...i', a)``, or to do a matrix-matrix
    product with the left-most indices instead of rightmost, you can do
    ``np.einsum('ij...,jk...->ik...', a, b)``.

    When there is only one operand, no axes are summed, and no output
    parameter is provided, a view into the operand is returned instead
    of a new array.  Thus, taking the diagonal as ``np.einsum('ii->i', a)``
    produces a view.

    An alternative way to provide the subscripts and operands is as
    ``einsum(op0, sublist0, op1, sublist1, ..., [sublistout])``. The examples
    below have corresponding `einsum` calls with the two parameter methods.

    Examples
    --------
    >>> a = np.arange(25).reshape(5,5)
    >>> b = np.arange(5)
    >>> c = np.arange(6).reshape(2,3)

    >>> np.einsum('ii', a)
    60
    >>> np.einsum(a, [0,0])
    60
    >>> np.trace(a)
    60

    >>> np.einsum('ii->i', a)
    array([ 0,  6, 12, 18, 24])
    >>> np.einsum(a, [0,0], [0])
    array([ 0,  6, 12, 18, 24])
    >>> np.diag(a)
    array([ 0,  6, 12, 18, 24])

    >>> np.einsum('ij,j', a, b)
    array([ 30,  80, 130, 180, 230])
    >>> np.einsum(a, [0,1], b, [1])
    array([ 30,  80, 130, 180, 230])
    >>> np.dot(a, b)
    array([ 30,  80, 130, 180, 230])

    >>> np.einsum('ji', c)
    array([[0, 3],
           [1, 4],
           [2, 5]])
    >>> np.einsum(c, [1,0])
    array([[0, 3],
           [1, 4],
           [2, 5]])
    >>> c.T
    array([[0, 3],
           [1, 4],
           [2, 5]])

    >>> np.einsum('..., ...', 3, c)
    array([[ 0,  3,  6],
           [ 9, 12, 15]])
    >>> np.einsum(3, [Ellipsis], c, [Ellipsis])
    array([[ 0,  3,  6],
           [ 9, 12, 15]])
    >>> np.multiply(3, c)
    array([[ 0,  3,  6],
           [ 9, 12, 15]])

    >>> np.einsum('i,i', b, b)
    30
    >>> np.einsum(b, [0], b, [0])
    30
    >>> np.inner(b,b)
    30

    >>> np.einsum('i,j', np.arange(2)+1, b)
    array([[0, 1, 2, 3, 4],
           [0, 2, 4, 6, 8]])
    >>> np.einsum(np.arange(2)+1, [0], b, [1])
    array([[0, 1, 2, 3, 4],
           [0, 2, 4, 6, 8]])
    >>> np.outer(np.arange(2)+1, b)
    array([[0, 1, 2, 3, 4],
           [0, 2, 4, 6, 8]])

    >>> np.einsum('i...->...', a)
    array([50, 55, 60, 65, 70])
    >>> np.einsum(a, [0,Ellipsis], [Ellipsis])
    array([50, 55, 60, 65, 70])
    >>> np.sum(a, axis=0)
    array([50, 55, 60, 65, 70])

    >>> a = np.arange(60.).reshape(3,4,5)
    >>> b = np.arange(24.).reshape(4,3,2)
    >>> np.einsum('ijk,jil->kl', a, b)
    array([[ 4400.,  4730.],
           [ 4532.,  4874.],
           [ 4664.,  5018.],
           [ 4796.,  5162.],
           [ 4928.,  5306.]])
    >>> np.einsum(a, [0,1,2], b, [1,0,3], [2,3])
    array([[ 4400.,  4730.],
           [ 4532.,  4874.],
           [ 4664.,  5018.],
           [ 4796.,  5162.],
           [ 4928.,  5306.]])
    >>> np.tensordot(a,b, axes=([1,0],[0,1]))
    array([[ 4400.,  4730.],
           [ 4532.,  4874.],
           [ 4664.,  5018.],
           [ 4796.,  5162.],
           [ 4928.,  5306.]])
    
    """
    # BUILTIN
    raise NotImplementedError


def expand_dims(a, axis):
    """Expand the shape of an array.

    Insert a new axis, corresponding to a given position in the array shape.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : int
        Position (amongst axes) where new axis is to be inserted.

    Returns
    -------
    res : ndarray
        Output array. The number of dimensions is one greater than that of
        the input array.

    See Also
    --------
    doc.indexing, atleast_1d, atleast_2d, atleast_3d

    Examples
    --------
    >>> x = np.array([1,2])
    >>> x.shape
    (2,)

    The following is equivalent to ``x[np.newaxis,:]`` or ``x[np.newaxis]``:

    >>> y = np.expand_dims(x, axis=0)
    >>> y
    array([[1, 2]])
    >>> y.shape
    (1, 2)

    >>> y = np.expand_dims(x, axis=1)  # Equivalent to x[:,newaxis]
    >>> y
    array([[1],
           [2]])
    >>> y.shape
    (2, 1)

    Note that some examples may use ``None`` instead of ``np.newaxis``.  These
    are the same objects:

    >>> np.newaxis is None
    True
    
    """
    raise NotImplementedError


def extract(condition, arr):
    """Return the elements of an array that satisfy some condition.

    This is equivalent to ``np.compress(ravel(condition), ravel(arr))``.  If
    `condition` is boolean ``np.extract`` is equivalent to ``arr[condition]``.

    Parameters
    ----------
    condition : array_like
        An array whose nonzero or True entries indicate the elements of `arr`
        to extract.
    arr : array_like
        Input array of the same size as `condition`.

    See Also
    --------
    take, put, putmask, compress

    Examples
    --------
    >>> arr = np.arange(12).reshape((3, 4))
    >>> arr
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11]])
    >>> condition = np.mod(arr, 3)==0
    >>> condition
    array([[ True, False, False,  True],
           [False, False,  True, False],
           [False,  True, False, False]], dtype=bool)
    >>> np.extract(condition, arr)
    array([0, 3, 6, 9])


    If `condition` is boolean:

    >>> arr[condition]
    array([0, 3, 6, 9])
    
    """
    raise NotImplementedError


def _fastCopyAndTranspose(a):
    """_fastCopyAndTranspose(a)
    
    """
    # BUILTIN
    raise NotImplementedError


def fill_diagonal(a, val):
    """Fill the main diagonal of the given array of any dimensionality.

    For an array `a` with ``a.ndim > 2``, the diagonal is the list of
    locations with indices ``a[i, i, ..., i]`` all identical. This function
    modifies the input array in-place, it does not return a value.

    Parameters
    ----------
    a : array, at least 2-D.
      Array whose diagonal is to be filled, it gets modified in-place.

    val : scalar
      Value to be written on the diagonal, its type must be compatible with
      that of the array a.

    See also
    --------
    diag_indices, diag_indices_from

    Notes
    -----
    .. versionadded:: 1.4.0

    This functionality can be obtained via `diag_indices`, but internally
    this version uses a much faster implementation that never constructs the
    indices and uses simple slicing.

    Examples
    --------
    >>> a = np.zeros((3, 3), int)
    >>> np.fill_diagonal(a, 5)
    >>> a
    array([[5, 0, 0],
           [0, 5, 0],
           [0, 0, 5]])

    The same function can operate on a 4-D array:

    >>> a = np.zeros((3, 3, 3, 3), int)
    >>> np.fill_diagonal(a, 4)

    We only show a few blocks for clarity:

    >>> a[0, 0]
    array([[4, 0, 0],
           [0, 0, 0],
           [0, 0, 0]])
    >>> a[1, 1]
    array([[0, 0, 0],
           [0, 4, 0],
           [0, 0, 0]])
    >>> a[2, 2]
    array([[0, 0, 0],
           [0, 0, 0],
           [0, 0, 4]])
    
    """
    raise NotImplementedError


def find_common_type(array_types, scalar_types):
    """Determine common type following standard coercion rules.

    Parameters
    ----------
    array_types : sequence
        A list of dtypes or dtype convertible objects representing arrays.
    scalar_types : sequence
        A list of dtypes or dtype convertible objects representing scalars.

    Returns
    -------
    datatype : dtype
        The common data type, which is the maximum of `array_types` ignoring
        `scalar_types`, unless the maximum of `scalar_types` is of a
        different kind (`dtype.kind`). If the kind is not understood, then
        None is returned.

    See Also
    --------
    dtype, common_type, can_cast, mintypecode

    Examples
    --------
    >>> np.find_common_type([], [np.int64, np.float32, np.complex])
    dtype('complex128')
    >>> np.find_common_type([np.int64, np.float32], [])
    dtype('float64')

    The standard casting rules ensure that a scalar cannot up-cast an
    array unless the scalar is of a fundamentally different kind of data
    (i.e. under a different hierarchy in the data type hierarchy) then
    the array:

    >>> np.find_common_type([np.float32], [np.int64, np.float64])
    dtype('float32')

    Complex is of a different type, so it up-casts the float in the
    `array_types` argument:

    >>> np.find_common_type([np.float32], [np.complex])
    dtype('complex128')

    Type specifier strings are convertible to dtypes and can therefore
    be used instead of dtypes:

    >>> np.find_common_type(['f4', 'f4', 'i4'], ['c8'])
    dtype('complex128')
    
    """
    raise NotImplementedError


def fix(x, y=None):
    """Round to nearest integer towards zero.

    Round an array of floats element-wise to nearest integer towards zero.
    The rounded values are returned as floats.

    Parameters
    ----------
    x : array_like
        An array of floats to be rounded
    y : ndarray, optional
        Output array

    Returns
    -------
    out : ndarray of floats
        The array of rounded numbers

    See Also
    --------
    trunc, floor, ceil
    around : Round to given number of decimals

    Examples
    --------
    >>> np.fix(3.14)
    3.0
    >>> np.fix(3)
    3.0
    >>> np.fix([2.1, 2.9, -2.1, -2.9])
    array([ 2.,  2., -2., -2.])
    
    """
    raise NotImplementedError


def flatnonzero(a):
    """Return indices that are non-zero in the flattened version of a.

    This is equivalent to a.ravel().nonzero()[0].

    Parameters
    ----------
    a : ndarray
        Input array.

    Returns
    -------
    res : ndarray
        Output array, containing the indices of the elements of `a.ravel()`
        that are non-zero.

    See Also
    --------
    nonzero : Return the indices of the non-zero elements of the input array.
    ravel : Return a 1-D array containing the elements of the input array.

    Examples
    --------
    >>> x = np.arange(-2, 3)
    >>> x
    array([-2, -1,  0,  1,  2])
    >>> np.flatnonzero(x)
    array([0, 1, 3, 4])

    Use the indices of the non-zero elements as an index array to extract
    these elements:

    >>> x.ravel()[np.flatnonzero(x)]
    array([-2, -1,  1,  2])
    
    """
    raise NotImplementedError


def fliplr(m):
    """Flip array in the left/right direction.

    Flip the entries in each row in the left/right direction.
    Columns are preserved, but appear in a different order than before.

    Parameters
    ----------
    m : array_like
        Input array.

    Returns
    -------
    f : ndarray
        A view of `m` with the columns reversed.  Since a view
        is returned, this operation is :math:`\mathcal O(1)`.

    See Also
    --------
    flipud : Flip array in the up/down direction.
    rot90 : Rotate array counterclockwise.

    Notes
    -----
    Equivalent to A[:,::-1]. Does not require the array to be
    two-dimensional.

    Examples
    --------
    >>> A = np.diag([1.,2.,3.])
    >>> A
    array([[ 1.,  0.,  0.],
           [ 0.,  2.,  0.],
           [ 0.,  0.,  3.]])
    >>> np.fliplr(A)
    array([[ 0.,  0.,  1.],
           [ 0.,  2.,  0.],
           [ 3.,  0.,  0.]])

    >>> A = np.random.randn(2,3,5)
    >>> np.all(np.fliplr(A)==A[:,::-1,...])
    True
    
    """
    raise NotImplementedError


def flipud(m):
    """Flip array in the up/down direction.

    Flip the entries in each column in the up/down direction.
    Rows are preserved, but appear in a different order than before.

    Parameters
    ----------
    m : array_like
        Input array.

    Returns
    -------
    out : array_like
        A view of `m` with the rows reversed.  Since a view is
        returned, this operation is :math:`\mathcal O(1)`.

    See Also
    --------
    fliplr : Flip array in the left/right direction.
    rot90 : Rotate array counterclockwise.

    Notes
    -----
    Equivalent to ``A[::-1,...]``.
    Does not require the array to be two-dimensional.

    Examples
    --------
    >>> A = np.diag([1.0, 2, 3])
    >>> A
    array([[ 1.,  0.,  0.],
           [ 0.,  2.,  0.],
           [ 0.,  0.,  3.]])
    >>> np.flipud(A)
    array([[ 0.,  0.,  3.],
           [ 0.,  2.,  0.],
           [ 1.,  0.,  0.]])

    >>> A = np.random.randn(2,3,5)
    >>> np.all(np.flipud(A)==A[::-1,...])
    True

    >>> np.flipud([1,2])
    array([2, 1])
    
    """
    raise NotImplementedError


def frombuffer(buffer, dtype=float, count=-1, offset=0):
    """frombuffer(buffer, dtype=float, count=-1, offset=0)

    Interpret a buffer as a 1-dimensional array.

    Parameters
    ----------
    buffer : buffer_like
        An object that exposes the buffer interface.
    dtype : data-type, optional
        Data-type of the returned array; default: float.
    count : int, optional
        Number of items to read. ``-1`` means all data in the buffer.
    offset : int, optional
        Start reading the buffer from this offset; default: 0.

    Notes
    -----
    If the buffer has data that is not in machine byte-order, this should
    be specified as part of the data-type, e.g.::

      >>> dt = np.dtype(int)
      >>> dt = dt.newbyteorder('>')
      >>> np.frombuffer(buf, dtype=dt)

    The data of the resulting array will not be byteswapped, but will be
    interpreted correctly.

    Examples
    --------
    >>> s = 'hello world'
    >>> np.frombuffer(s, dtype='S1', count=5, offset=6)
    array(['w', 'o', 'r', 'l', 'd'],
          dtype='|S1')
    
    """
    # BUILTIN
    raise NotImplementedError


def fromfile(file, dtype=float, count=-1, sep=''):
    """fromfile(file, dtype=float, count=-1, sep='')

    Construct an array from data in a text or binary file.

    A highly efficient way of reading binary data with a known data-type,
    as well as parsing simply formatted text files.  Data written using the
    `tofile` method can be read using this function.

    Parameters
    ----------
    file : file or str
        Open file object or filename.
    dtype : data-type
        Data type of the returned array.
        For binary files, it is used to determine the size and byte-order
        of the items in the file.
    count : int
        Number of items to read. ``-1`` means all items (i.e., the complete
        file).
    sep : str
        Separator between items if file is a text file.
        Empty ("") separator means the file should be treated as binary.
        Spaces (" ") in the separator match zero or more whitespace characters.
        A separator consisting only of spaces must match at least one
        whitespace.

    See also
    --------
    load, save
    ndarray.tofile
    loadtxt : More flexible way of loading data from a text file.

    Notes
    -----
    Do not rely on the combination of `tofile` and `fromfile` for
    data storage, as the binary files generated are are not platform
    independent.  In particular, no byte-order or data-type information is
    saved.  Data can be stored in the platform independent ``.npy`` format
    using `save` and `load` instead.

    Examples
    --------
    Construct an ndarray:

    >>> dt = np.dtype([('time', [('min', int), ('sec', int)]),
    ...                ('temp', float)])
    >>> x = np.zeros((1,), dtype=dt)
    >>> x['time']['min'] = 10; x['temp'] = 98.25
    >>> x
    array([((10, 0), 98.25)],
          dtype=[('time', [('min', '<i4'), ('sec', '<i4')]), ('temp', '<f8')])

    Save the raw data to disk:

    >>> import os
    >>> fname = os.tmpnam()
    >>> x.tofile(fname)

    Read the raw data from disk:

    >>> np.fromfile(fname, dtype=dt)
    array([((10, 0), 98.25)],
          dtype=[('time', [('min', '<i4'), ('sec', '<i4')]), ('temp', '<f8')])

    The recommended way to store and load data:

    >>> np.save(fname, x)
    >>> np.load(fname + '.npy')
    array([((10, 0), 98.25)],
          dtype=[('time', [('min', '<i4'), ('sec', '<i4')]), ('temp', '<f8')])
    
    """
    # BUILTIN
    raise NotImplementedError


def fromiter(iterable, dtype, count=-1):
    """fromiter(iterable, dtype, count=-1)

    Create a new 1-dimensional array from an iterable object.

    Parameters
    ----------
    iterable : iterable object
        An iterable object providing data for the array.
    dtype : data-type
        The data-type of the returned array.
    count : int, optional
        The number of items to read from *iterable*.  The default is -1,
        which means all data is read.

    Returns
    -------
    out : ndarray
        The output array.

    Notes
    -----
    Specify `count` to improve performance.  It allows ``fromiter`` to
    pre-allocate the output array, instead of resizing it on demand.

    Examples
    --------
    >>> iterable = (x*x for x in range(5))
    >>> np.fromiter(iterable, np.float)
    array([  0.,   1.,   4.,   9.,  16.])
    
    """
    # BUILTIN
    raise NotImplementedError


def frompyfunc(func, nin, nout):
    """frompyfunc(func, nin, nout)

    Takes an arbitrary Python function and returns a Numpy ufunc.

    Can be used, for example, to add broadcasting to a built-in Python
    function (see Examples section).

    Parameters
    ----------
    func : Python function object
        An arbitrary Python function.
    nin : int
        The number of input arguments.
    nout : int
        The number of objects returned by `func`.

    Returns
    -------
    out : ufunc
        Returns a Numpy universal function (``ufunc``) object.

    Notes
    -----
    The returned ufunc always returns PyObject arrays.

    Examples
    --------
    Use frompyfunc to add broadcasting to the Python function ``oct``:

    >>> oct_array = np.frompyfunc(oct, 1, 1)
    >>> oct_array(np.array((10, 30, 100)))
    array([012, 036, 0144], dtype=object)
    >>> np.array((oct(10), oct(30), oct(100))) # for comparison
    array(['012', '036', '0144'],
          dtype='|S4')
    
    """
    # BUILTIN
    raise NotImplementedError


def fromregex(file, regexp, dtype):
    """Construct an array from a text file, using regular expression parsing.

    The returned array is always a structured array, and is constructed from
    all matches of the regular expression in the file. Groups in the regular
    expression are converted to fields of the structured array.

    Parameters
    ----------
    file : str or file
        File name or file object to read.
    regexp : str or regexp
        Regular expression used to parse the file.
        Groups in the regular expression correspond to fields in the dtype.
    dtype : dtype or list of dtypes
        Dtype for the structured array.

    Returns
    -------
    output : ndarray
        The output array, containing the part of the content of `file` that
        was matched by `regexp`. `output` is always a structured array.

    Raises
    ------
    TypeError
        When `dtype` is not a valid dtype for a structured array.

    See Also
    --------
    fromstring, loadtxt

    Notes
    -----
    Dtypes for structured arrays can be specified in several forms, but all
    forms specify at least the data type and field name. For details see
    `doc.structured_arrays`.

    Examples
    --------
    >>> f = open('test.dat', 'w')
    >>> f.write("1312 foo\n1534  bar\n444   qux")
    >>> f.close()

    >>> regexp = r"(\d+)\s+(...)"  # match [digits, whitespace, anything]
    >>> output = np.fromregex('test.dat', regexp,
    ...                       [('num', np.int64), ('key', 'S3')])
    >>> output
    array([(1312L, 'foo'), (1534L, 'bar'), (444L, 'qux')],
          dtype=[('num', '<i8'), ('key', '|S3')])
    >>> output['num']
    array([1312, 1534,  444], dtype=int64)
    
    """
    raise NotImplementedError


def fromstring(string, dtype=float, count=-1, sep=''):
    """fromstring(string, dtype=float, count=-1, sep='')

    A new 1-D array initialized from raw binary or text data in a string.

    Parameters
    ----------
    string : str
        A string containing the data.
    dtype : data-type, optional
        The data type of the array; default: float.  For binary input data,
        the data must be in exactly this format.
    count : int, optional
        Read this number of `dtype` elements from the data.  If this is
        negative (the default), the count will be determined from the
        length of the data.
    sep : str, optional
        If not provided or, equivalently, the empty string, the data will
        be interpreted as binary data; otherwise, as ASCII text with
        decimal numbers.  Also in this latter case, this argument is
        interpreted as the string separating numbers in the data; extra
        whitespace between elements is also ignored.

    Returns
    -------
    arr : ndarray
        The constructed array.

    Raises
    ------
    ValueError
        If the string is not the correct size to satisfy the requested
        `dtype` and `count`.

    See Also
    --------
    frombuffer, fromfile, fromiter

    Examples
    --------
    >>> np.fromstring('\x01\x02', dtype=np.uint8)
    array([1, 2], dtype=uint8)
    >>> np.fromstring('1 2', dtype=int, sep=' ')
    array([1, 2])
    >>> np.fromstring('1, 2', dtype=int, sep=',')
    array([1, 2])
    >>> np.fromstring('\x01\x02\x03\x04\x05', dtype=np.uint8, count=3)
    array([1, 2, 3], dtype=uint8)
    
    """
    # BUILTIN
    raise NotImplementedError


def fv(rate, nper, pmt, pv, when='end'):
    """Compute the future value.

    Given:
     * a present value, `pv`
     * an interest `rate` compounded once per period, of which
       there are
     * `nper` total
     * a (fixed) payment, `pmt`, paid either
     * at the beginning (`when` = {'begin', 1}) or the end
       (`when` = {'end', 0}) of each period

    Return:
       the value at the end of the `nper` periods

    Parameters
    ----------
    rate : scalar or array_like of shape(M, )
        Rate of interest as decimal (not per cent) per period
    nper : scalar or array_like of shape(M, )
        Number of compounding periods
    pmt : scalar or array_like of shape(M, )
        Payment
    pv : scalar or array_like of shape(M, )
        Present value
    when : {{'begin', 1}, {'end', 0}}, {string, int}, optional
        When payments are due ('begin' (1) or 'end' (0)).
        Defaults to {'end', 0}.

    Returns
    -------
    out : ndarray
        Future values.  If all input is scalar, returns a scalar float.  If
        any input is array_like, returns future values for each input element.
        If multiple inputs are array_like, they all must have the same shape.

    Notes
    -----
    The future value is computed by solving the equation::

     fv +
     pv*(1+rate)**nper +
     pmt*(1 + rate*when)/rate*((1 + rate)**nper - 1) == 0

    or, when ``rate == 0``::

     fv + pv + pmt * nper == 0

    References
    ----------
    .. [WRW] Wheeler, D. A., E. Rathke, and R. Weir (Eds.) (2009, May).
       Open Document Format for Office Applications (OpenDocument)v1.2,
       Part 2: Recalculated Formula (OpenFormula) Format - Annotated Version,
       Pre-Draft 12. Organization for the Advancement of Structured Information
       Standards (OASIS). Billerica, MA, USA. [ODT Document].
       Available:
       http://www.oasis-open.org/committees/documents.php?wg_abbrev=office-formula
       OpenDocument-formula-20090508.odt

    Examples
    --------
    What is the future value after 10 years of saving $100 now, with
    an additional monthly savings of $100.  Assume the interest rate is
    5% (annually) compounded monthly?

    >>> np.fv(0.05/12, 10*12, -100, -100)
    15692.928894335748

    By convention, the negative sign represents cash flow out (i.e. money not
    available today).  Thus, saving $100 a month at 5% annual interest leads
    to $15,692.93 available to spend in 10 years.

    If any input is array_like, returns an array of equal shape.  Let's
    compare different interest rates from the example above.

    >>> a = np.array((0.05, 0.06, 0.07))/12
    >>> np.fv(a, 10*12, -100, -100)
    array([ 15692.92889434,  16569.87435405,  17509.44688102])
    
    """
    raise NotImplementedError


def genfromtxt(fname, dtype=float, comments='#', delimiter=None, skiprows=0, skip_header=0, skip_footer=0, converters=None, missing='', missing_values=None, filling_values=None, usecols=None, names=None, excludelist=None, deletechars=None, replace_space='_', autostrip=False, case_sensitive=True, defaultfmt='f%i', unpack=None, usemask=False, loose=True, invalid_raise=True):
    """Load data from a text file, with missing values handled as specified.

    Each line past the first `skip_header` lines is split at the `delimiter`
    character, and characters following the `comments` character are discarded.

    Parameters
    ----------
    fname : file or str
        File, filename, or generator to read.  If the filename extension is
        `.gz` or `.bz2`, the file is first decompressed. Note that
        generators must return byte strings in Python 3k.
    dtype : dtype, optional
        Data type of the resulting array.
        If None, the dtypes will be determined by the contents of each
        column, individually.
    comments : str, optional
        The character used to indicate the start of a comment.
        All the characters occurring on a line after a comment are discarded
    delimiter : str, int, or sequence, optional
        The string used to separate values.  By default, any consecutive
        whitespaces act as delimiter.  An integer or sequence of integers
        can also be provided as width(s) of each field.
    skip_header : int, optional
        The numbers of lines to skip at the beginning of the file.
    skip_footer : int, optional
        The numbers of lines to skip at the end of the file
    converters : variable, optional
        The set of functions that convert the data of a column to a value.
        The converters can also be used to provide a default value
        for missing data: ``converters = {3: lambda s: float(s or 0)}``.
    missing_values : variable, optional
        The set of strings corresponding to missing data.
    filling_values : variable, optional
        The set of values to be used as default when the data are missing.
    usecols : sequence, optional
        Which columns to read, with 0 being the first.  For example,
        ``usecols = (1, 4, 5)`` will extract the 2nd, 5th and 6th columns.
    names : {None, True, str, sequence}, optional
        If `names` is True, the field names are read from the first valid line
        after the first `skip_header` lines.
        If `names` is a sequence or a single-string of comma-separated names,
        the names will be used to define the field names in a structured dtype.
        If `names` is None, the names of the dtype fields will be used, if any.
    excludelist : sequence, optional
        A list of names to exclude. This list is appended to the default list
        ['return','file','print']. Excluded names are appended an underscore:
        for example, `file` would become `file_`.
    deletechars : str, optional
        A string combining invalid characters that must be deleted from the
        names.
    defaultfmt : str, optional
        A format used to define default field names, such as "f%i" or "f_%02i".
    autostrip : bool, optional
        Whether to automatically strip white spaces from the variables.
    replace_space : char, optional
        Character(s) used in replacement of white spaces in the variables names.
        By default, use a '_'.
    case_sensitive : {True, False, 'upper', 'lower'}, optional
        If True, field names are case sensitive.
        If False or 'upper', field names are converted to upper case.
        If 'lower', field names are converted to lower case.
    unpack : bool, optional
        If True, the returned array is transposed, so that arguments may be
        unpacked using ``x, y, z = loadtxt(...)``
    usemask : bool, optional
        If True, return a masked array.
        If False, return a regular array.
    invalid_raise : bool, optional
        If True, an exception is raised if an inconsistency is detected in the
        number of columns.
        If False, a warning is emitted and the offending lines are skipped.

    Returns
    -------
    out : ndarray
        Data read from the text file. If `usemask` is True, this is a
        masked array.

    See Also
    --------
    numpy.loadtxt : equivalent function when no data is missing.

    Notes
    -----
    * When spaces are used as delimiters, or when no delimiter has been given
      as input, there should not be any missing data between two fields.
    * When the variables are named (either by a flexible dtype or with `names`,
      there must not be any header in the file (else a ValueError
      exception is raised).
    * Individual values are not stripped of spaces by default.
      When using a custom converter, make sure the function does remove spaces.

    Examples
    ---------
    >>> from StringIO import StringIO
    >>> import numpy as np

    Comma delimited file with mixed dtype

    >>> s = StringIO("1,1.3,abcde")
    >>> data = np.genfromtxt(s, dtype=[('myint','i8'),('myfloat','f8'),
    ... ('mystring','S5')], delimiter=",")
    >>> data
    array((1, 1.3, 'abcde'),
          dtype=[('myint', '<i8'), ('myfloat', '<f8'), ('mystring', '|S5')])

    Using dtype = None

    >>> s.seek(0) # needed for StringIO example only
    >>> data = np.genfromtxt(s, dtype=None,
    ... names = ['myint','myfloat','mystring'], delimiter=",")
    >>> data
    array((1, 1.3, 'abcde'),
          dtype=[('myint', '<i8'), ('myfloat', '<f8'), ('mystring', '|S5')])

    Specifying dtype and names

    >>> s.seek(0)
    >>> data = np.genfromtxt(s, dtype="i8,f8,S5",
    ... names=['myint','myfloat','mystring'], delimiter=",")
    >>> data
    array((1, 1.3, 'abcde'),
          dtype=[('myint', '<i8'), ('myfloat', '<f8'), ('mystring', '|S5')])

    An example with fixed-width columns

    >>> s = StringIO("11.3abcde")
    >>> data = np.genfromtxt(s, dtype=None, names=['intvar','fltvar','strvar'],
    ...     delimiter=[1,3,5])
    >>> data
    array((1, 1.3, 'abcde'),
          dtype=[('intvar', '<i8'), ('fltvar', '<f8'), ('strvar', '|S5')])
    
    """
    raise NotImplementedError


def get_array_wrap():
    """Find the wrapper for the array with the highest priority.

    In case of ties, leftmost wins. If no wrapper is found, return None
    
    """
    raise NotImplementedError


def get_include():
    """Return the directory that contains the NumPy \*.h header files.

    Extension modules that need to compile against NumPy should use this
    function to locate the appropriate include directory.

    Notes
    -----
    When using ``distutils``, for example in ``setup.py``.
    ::

        import numpy as np
        ...
        Extension('extension_name', ...
                include_dirs=[np.get_include()])
        ...
    
    """
    raise NotImplementedError


def get_numarray_include(type=None):
    """Return the directory that contains the numarray \*.h header files.

    Extension modules that need to compile against numarray should use this
    function to locate the appropriate include directory.

    Parameters
    ----------
    type : any, optional
        If `type` is not None, the location of the NumPy headers is returned
        as well.

    Returns
    -------
    dirs : str or list of str
        If `type` is None, `dirs` is a string containing the path to the
        numarray headers.
        If `type` is not None, `dirs` is a list of strings with first the
        path(s) to the numarray headers, followed by the path to the NumPy
        headers.

    Notes
    -----
    Useful when using ``distutils``, for example in ``setup.py``.
    ::

        import numpy as np
        ...
        Extension('extension_name', ...
                include_dirs=[np.get_numarray_include()])
        ...
    
    """
    raise NotImplementedError


def get_printoptions():
    """Return the current print options.

    Returns
    -------
    print_opts : dict
        Dictionary of current print options with keys

          - precision : int
          - threshold : int
          - edgeitems : int
          - linewidth : int
          - suppress : bool
          - nanstr : str
          - infstr : str

        For a full description of these options, see `set_printoptions`.

    See Also
    --------
    set_printoptions, set_string_function
    
    """
    raise NotImplementedError


def getbuffer():
    """getbuffer(obj [,offset[, size]])

    Create a buffer object from the given object referencing a slice of
    length size starting at offset.

    Default is the entire buffer. A read-write buffer is attempted followed
    by a read-only buffer.

    Parameters
    ----------
    obj : object

    offset : int, optional

    size : int, optional

    Returns
    -------
    buffer_obj : buffer

    Examples
    --------
    >>> buf = np.getbuffer(np.ones(5), 1, 3)
    >>> len(buf)
    3
    >>> buf[0]
    '\x00'
    >>> buf
    <read-write buffer for 0x8af1e70, size 3, offset 1 at 0x8ba4ec0>
    
    """
    # BUILTIN
    raise NotImplementedError


def getbufsize():
    """Return the size of the buffer used in ufuncs.
    
    """
    raise NotImplementedError


def geterr():
    """Get the current way of handling floating-point errors.

    Returns
    -------
    res : dict
        A dictionary with keys "divide", "over", "under", and "invalid",
        whose values are from the strings "ignore", "print", "log", "warn",
        "raise", and "call". The keys represent possible floating-point
        exceptions, and the values define how these exceptions are handled.

    See Also
    --------
    geterrcall, seterr, seterrcall

    Notes
    -----
    For complete documentation of the types of floating-point exceptions and
    treatment options, see `seterr`.

    Examples
    --------
    >>> np.geterr()
    {'over': 'warn', 'divide': 'warn', 'invalid': 'warn',
    'under': 'ignore'}
    >>> np.arange(3.) / np.arange(3.)
    array([ NaN,   1.,   1.])

    >>> oldsettings = np.seterr(all='warn', over='raise')
    >>> np.geterr()
    {'over': 'raise', 'divide': 'warn', 'invalid': 'warn', 'under': 'warn'}
    >>> np.arange(3.) / np.arange(3.)
    __main__:1: RuntimeWarning: invalid value encountered in divide
    array([ NaN,   1.,   1.])
    
    """
    raise NotImplementedError


def geterrcall():
    """Return the current callback function used on floating-point errors.

    When the error handling for a floating-point error (one of "divide",
    "over", "under", or "invalid") is set to 'call' or 'log', the function
    that is called or the log instance that is written to is returned by
    `geterrcall`. This function or log instance has been set with
    `seterrcall`.

    Returns
    -------
    errobj : callable, log instance or None
        The current error handler. If no handler was set through `seterrcall`,
        ``None`` is returned.

    See Also
    --------
    seterrcall, seterr, geterr

    Notes
    -----
    For complete documentation of the types of floating-point exceptions and
    treatment options, see `seterr`.

    Examples
    --------
    >>> np.geterrcall()  # we did not yet set a handler, returns None

    >>> oldsettings = np.seterr(all='call')
    >>> def err_handler(type, flag):
    ...     print "Floating point error (%s), with flag %s" % (type, flag)
    >>> oldhandler = np.seterrcall(err_handler)
    >>> np.array([1, 2, 3]) / 0.0
    Floating point error (divide by zero), with flag 1
    array([ Inf,  Inf,  Inf])

    >>> cur_handler = np.geterrcall()
    >>> cur_handler is err_handler
    True
    
    """
    raise NotImplementedError


def geterrobj():
    """geterrobj()

    Return the current object that defines floating-point error handling.

    The error object contains all information that defines the error handling
    behavior in Numpy. `geterrobj` is used internally by the other
    functions that get and set error handling behavior (`geterr`, `seterr`,
    `geterrcall`, `seterrcall`).

    Returns
    -------
    errobj : list
        The error object, a list containing three elements:
        [internal numpy buffer size, error mask, error callback function].

        The error mask is a single integer that holds the treatment information
        on all four floating point errors. The information for each error type
        is contained in three bits of the integer. If we print it in base 8, we
        can see what treatment is set for "invalid", "under", "over", and
        "divide" (in that order). The printed string can be interpreted with

        * 0 : 'ignore'
        * 1 : 'warn'
        * 2 : 'raise'
        * 3 : 'call'
        * 4 : 'print'
        * 5 : 'log'

    See Also
    --------
    seterrobj, seterr, geterr, seterrcall, geterrcall
    getbufsize, setbufsize

    Notes
    -----
    For complete documentation of the types of floating-point exceptions and
    treatment options, see `seterr`.

    Examples
    --------
    >>> np.geterrobj()  # first get the defaults
    [10000, 0, None]

    >>> def err_handler(type, flag):
    ...     print "Floating point error (%s), with flag %s" % (type, flag)
    ...
    >>> old_bufsize = np.setbufsize(20000)
    >>> old_err = np.seterr(divide='raise')
    >>> old_handler = np.seterrcall(err_handler)
    >>> np.geterrobj()
    [20000, 2, <function err_handler at 0x91dcaac>]

    >>> old_err = np.seterr(all='ignore')
    >>> np.base_repr(np.geterrobj()[1], 8)
    '0'
    >>> old_err = np.seterr(divide='warn', over='log', under='call',
                            invalid='print')
    >>> np.base_repr(np.geterrobj()[1], 8)
    '4351'
    
    """
    # BUILTIN
    raise NotImplementedError


def gradient(f):
    """Return the gradient of an N-dimensional array.

    The gradient is computed using central differences in the interior
    and first differences at the boundaries. The returned gradient hence has
    the same shape as the input array.

    Parameters
    ----------
    f : array_like
      An N-dimensional array containing samples of a scalar function.
    `*varargs` : scalars
      0, 1, or N scalars specifying the sample distances in each direction,
      that is: `dx`, `dy`, `dz`, ... The default distance is 1.


    Returns
    -------
    g : ndarray
      N arrays of the same shape as `f` giving the derivative of `f` with
      respect to each dimension.

    Examples
    --------
    >>> x = np.array([1, 2, 4, 7, 11, 16], dtype=np.float)
    >>> np.gradient(x)
    array([ 1. ,  1.5,  2.5,  3.5,  4.5,  5. ])
    >>> np.gradient(x, 2)
    array([ 0.5 ,  0.75,  1.25,  1.75,  2.25,  2.5 ])

    >>> np.gradient(np.array([[1, 2, 6], [3, 4, 5]], dtype=np.float))
    [array([[ 2.,  2., -1.],
           [ 2.,  2., -1.]]),
    array([[ 1. ,  2.5,  4. ],
           [ 1. ,  1. ,  1. ]])]
    
    """
    raise NotImplementedError


def hamming(M):
    """Return the Hamming window.

    The Hamming window is a taper formed by using a weighted cosine.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an
        empty array is returned.

    Returns
    -------
    out : ndarray
        The window, normalized to one (the value one
        appears only if the number of samples is odd).

    See Also
    --------
    bartlett, blackman, hanning, kaiser

    Notes
    -----
    The Hamming window is defined as

    .. math::  w(n) = 0.54 + 0.46cos\left(\frac{2\pi{n}}{M-1}\right)
               \qquad 0 \leq n \leq M-1

    The Hamming was named for R. W. Hamming, an associate of J. W. Tukey and
    is described in Blackman and Tukey. It was recommended for smoothing the
    truncated autocovariance function in the time domain.
    Most references to the Hamming window come from the signal processing
    literature, where it is used as one of many windowing functions for
    smoothing values.  It is also known as an apodization (which means
    "removing the foot", i.e. smoothing discontinuities at the beginning
    and end of the sampled signal) or tapering function.

    References
    ----------
    .. [1] Blackman, R.B. and Tukey, J.W., (1958) The measurement of power
           spectra, Dover Publications, New York.
    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics", The
           University of Alberta Press, 1975, pp. 109-110.
    .. [3] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function
    .. [4] W.H. Press,  B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
           "Numerical Recipes", Cambridge University Press, 1986, page 425.

    Examples
    --------
    >>> np.hamming(12)
    array([ 0.08      ,  0.15302337,  0.34890909,  0.60546483,  0.84123594,
            0.98136677,  0.98136677,  0.84123594,  0.60546483,  0.34890909,
            0.15302337,  0.08      ])

    Plot the window and the frequency response:

    >>> from numpy.fft import fft, fftshift
    >>> import matplotlib.pyplot as plt

    >>> window = np.hamming(51)
    >>> plt.plot(window)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Hamming window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Sample")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    >>> plt.figure()
    <matplotlib.figure.Figure object at 0x...>
    >>> A = fft(window, 2048) / 25.5
    >>> mag = np.abs(fftshift(A))
    >>> freq = np.linspace(-0.5, 0.5, len(A))
    >>> response = 20 * np.log10(mag)
    >>> response = np.clip(response, -100, 100)
    >>> plt.plot(freq, response)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Frequency response of Hamming window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Magnitude [dB]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.axis('tight')
    (-0.5, 0.5, -100.0, ...)
    >>> plt.show()
    
    """
    raise NotImplementedError


def hanning(M):
    """Return the Hanning window.

    The Hanning window is a taper formed by using a weighted cosine.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an
        empty array is returned.

    Returns
    -------
    out : ndarray, shape(M,)
        The window, normalized to one (the value one
        appears only if `M` is odd).

    See Also
    --------
    bartlett, blackman, hamming, kaiser

    Notes
    -----
    The Hanning window is defined as

    .. math::  w(n) = 0.5 - 0.5cos\left(\frac{2\pi{n}}{M-1}\right)
               \qquad 0 \leq n \leq M-1

    The Hanning was named for Julius van Hann, an Austrian meterologist. It is
    also known as the Cosine Bell. Some authors prefer that it be called a
    Hann window, to help avoid confusion with the very similar Hamming window.

    Most references to the Hanning window come from the signal processing
    literature, where it is used as one of many windowing functions for
    smoothing values.  It is also known as an apodization (which means
    "removing the foot", i.e. smoothing discontinuities at the beginning
    and end of the sampled signal) or tapering function.

    References
    ----------
    .. [1] Blackman, R.B. and Tukey, J.W., (1958) The measurement of power
           spectra, Dover Publications, New York.
    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics",
           The University of Alberta Press, 1975, pp. 106-108.
    .. [3] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function
    .. [4] W.H. Press,  B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling,
           "Numerical Recipes", Cambridge University Press, 1986, page 425.

    Examples
    --------
    >>> from numpy import hanning
    >>> hanning(12)
    array([ 0.        ,  0.07937323,  0.29229249,  0.57115742,  0.82743037,
            0.97974649,  0.97974649,  0.82743037,  0.57115742,  0.29229249,
            0.07937323,  0.        ])

    Plot the window and its frequency response:

    >>> from numpy.fft import fft, fftshift
    >>> import matplotlib.pyplot as plt

    >>> window = np.hanning(51)
    >>> plt.plot(window)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Hann window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Sample")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    >>> plt.figure()
    <matplotlib.figure.Figure object at 0x...>
    >>> A = fft(window, 2048) / 25.5
    >>> mag = abs(fftshift(A))
    >>> freq = np.linspace(-0.5,0.5,len(A))
    >>> response = 20*np.log10(mag)
    >>> response = np.clip(response,-100,100)
    >>> plt.plot(freq, response)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Frequency response of the Hann window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Magnitude [dB]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.axis('tight')
    (-0.5, 0.5, -100.0, ...)
    >>> plt.show()
    
    """
    raise NotImplementedError


def histogram(a, bins=10, range=None, normed=False, weights=None, density=None):
    """Compute the histogram of a set of data.

    Parameters
    ----------
    a : array_like
        Input data. The histogram is computed over the flattened array.
    bins : int or sequence of scalars, optional
        If `bins` is an int, it defines the number of equal-width
        bins in the given range (10, by default). If `bins` is a sequence,
        it defines the bin edges, including the rightmost edge, allowing
        for non-uniform bin widths.
    range : (float, float), optional
        The lower and upper range of the bins.  If not provided, range
        is simply ``(a.min(), a.max())``.  Values outside the range are
        ignored.
    normed : bool, optional
        This keyword is deprecated in Numpy 1.6 due to confusing/buggy
        behavior. It will be removed in Numpy 2.0. Use the density keyword
        instead.
        If False, the result will contain the number of samples
        in each bin.  If True, the result is the value of the
        probability *density* function at the bin, normalized such that
        the *integral* over the range is 1. Note that this latter behavior is
        known to be buggy with unequal bin widths; use `density` instead.
    weights : array_like, optional
        An array of weights, of the same shape as `a`.  Each value in `a`
        only contributes its associated weight towards the bin count
        (instead of 1).  If `normed` is True, the weights are normalized,
        so that the integral of the density over the range remains 1
    density : bool, optional
        If False, the result will contain the number of samples
        in each bin.  If True, the result is the value of the
        probability *density* function at the bin, normalized such that
        the *integral* over the range is 1. Note that the sum of the
        histogram values will not be equal to 1 unless bins of unity
        width are chosen; it is not a probability *mass* function.
        Overrides the `normed` keyword if given.

    Returns
    -------
    hist : array
        The values of the histogram. See `normed` and `weights` for a
        description of the possible semantics.
    bin_edges : array of dtype float
        Return the bin edges ``(length(hist)+1)``.


    See Also
    --------
    histogramdd, bincount, searchsorted, digitize

    Notes
    -----
    All but the last (righthand-most) bin is half-open.  In other words, if
    `bins` is::

      [1, 2, 3, 4]

    then the first bin is ``[1, 2)`` (including 1, but excluding 2) and the
    second ``[2, 3)``.  The last bin, however, is ``[3, 4]``, which *includes*
    4.

    Examples
    --------
    >>> np.histogram([1, 2, 1], bins=[0, 1, 2, 3])
    (array([0, 2, 1]), array([0, 1, 2, 3]))
    >>> np.histogram(np.arange(4), bins=np.arange(5), density=True)
    (array([ 0.25,  0.25,  0.25,  0.25]), array([0, 1, 2, 3, 4]))
    >>> np.histogram([[1, 2, 1], [1, 0, 1]], bins=[0,1,2,3])
    (array([1, 4, 1]), array([0, 1, 2, 3]))

    >>> a = np.arange(5)
    >>> hist, bin_edges = np.histogram(a, density=True)
    >>> hist
    array([ 0.5,  0. ,  0.5,  0. ,  0. ,  0.5,  0. ,  0.5,  0. ,  0.5])
    >>> hist.sum()
    2.4999999999999996
    >>> np.sum(hist*np.diff(bin_edges))
    1.0
    
    """
    raise NotImplementedError


def histogram2d(x, y, bins=10, range=None, normed=False, weights=None):
    """Compute the bi-dimensional histogram of two data samples.

    Parameters
    ----------
    x : array_like, shape(N,)
        A sequence of values to be histogrammed along the first dimension.
    y : array_like, shape(M,)
        A sequence of values to be histogrammed along the second dimension.
    bins : int or [int, int] or array_like or [array, array], optional
        The bin specification:

          * If int, the number of bins for the two dimensions (nx=ny=bins).
          * If [int, int], the number of bins in each dimension (nx, ny = bins).
          * If array_like, the bin edges for the two dimensions (x_edges=y_edges=bins).
          * If [array, array], the bin edges in each dimension (x_edges, y_edges = bins).

    range : array_like, shape(2,2), optional
        The leftmost and rightmost edges of the bins along each dimension
        (if not specified explicitly in the `bins` parameters):
        ``[[xmin, xmax], [ymin, ymax]]``. All values outside of this range
        will be considered outliers and not tallied in the histogram.
    normed : bool, optional
        If False, returns the number of samples in each bin. If True, returns
        the bin density, i.e. the bin count divided by the bin area.
    weights : array_like, shape(N,), optional
        An array of values ``w_i`` weighing each sample ``(x_i, y_i)``. Weights
        are normalized to 1 if `normed` is True. If `normed` is False, the
        values of the returned histogram are equal to the sum of the weights
        belonging to the samples falling into each bin.

    Returns
    -------
    H : ndarray, shape(nx, ny)
        The bi-dimensional histogram of samples `x` and `y`. Values in `x`
        are histogrammed along the first dimension and values in `y` are
        histogrammed along the second dimension.
    xedges : ndarray, shape(nx,)
        The bin edges along the first dimension.
    yedges : ndarray, shape(ny,)
        The bin edges along the second dimension.

    See Also
    --------
    histogram: 1D histogram
    histogramdd: Multidimensional histogram

    Notes
    -----
    When `normed` is True, then the returned histogram is the sample density,
    defined such that:

    .. math::
      \sum_{i=0}^{nx-1} \sum_{j=0}^{ny-1} H_{i,j} \Delta x_i \Delta y_j = 1

    where `H` is the histogram array and :math:`\Delta x_i \Delta y_i`
    the area of bin `{i,j}`.

    Please note that the histogram does not follow the Cartesian convention
    where `x` values are on the abcissa and `y` values on the ordinate axis.
    Rather, `x` is histogrammed along the first dimension of the array
    (vertical), and `y` along the second dimension of the array (horizontal).
    This ensures compatibility with `histogramdd`.

    Examples
    --------
    >>> x, y = np.random.randn(2, 100)
    >>> H, xedges, yedges = np.histogram2d(x, y, bins=(5, 8))
    >>> H.shape, xedges.shape, yedges.shape
    ((5, 8), (6,), (9,))

    We can now use the Matplotlib to visualize this 2-dimensional histogram:

    >>> extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    >>> import matplotlib.pyplot as plt
    >>> plt.imshow(H, extent=extent, interpolation='nearest')
    <matplotlib.image.AxesImage object at ...>
    >>> plt.colorbar()
    <matplotlib.colorbar.Colorbar instance at ...>
    >>> plt.show()
    
    """
    raise NotImplementedError


def histogramdd(sample, bins=10, range=None, normed=False, weights=None):
    """Compute the multidimensional histogram of some data.

    Parameters
    ----------
    sample : array_like
        The data to be histogrammed. It must be an (N,D) array or data
        that can be converted to such. The rows of the resulting array
        are the coordinates of points in a D dimensional polytope.
    bins : sequence or int, optional
        The bin specification:

        * A sequence of arrays describing the bin edges along each dimension.
        * The number of bins for each dimension (nx, ny, ... =bins)
        * The number of bins for all dimensions (nx=ny=...=bins).

    range : sequence, optional
        A sequence of lower and upper bin edges to be used if the edges are
        not given explicitely in `bins`. Defaults to the minimum and maximum
        values along each dimension.
    normed : bool, optional
        If False, returns the number of samples in each bin. If True, returns
        the bin density, ie, the bin count divided by the bin hypervolume.
    weights : array_like (N,), optional
        An array of values `w_i` weighing each sample `(x_i, y_i, z_i, ...)`.
        Weights are normalized to 1 if normed is True. If normed is False, the
        values of the returned histogram are equal to the sum of the weights
        belonging to the samples falling into each bin.

    Returns
    -------
    H : ndarray
        The multidimensional histogram of sample x. See normed and weights for
        the different possible semantics.
    edges : list
        A list of D arrays describing the bin edges for each dimension.

    See Also
    --------
    histogram: 1-D histogram
    histogram2d: 2-D histogram

    Examples
    --------
    >>> r = np.random.randn(100,3)
    >>> H, edges = np.histogramdd(r, bins = (5, 8, 4))
    >>> H.shape, edges[0].size, edges[1].size, edges[2].size
    ((5, 8, 4), 6, 9, 5)
    
    """
    raise NotImplementedError


def hsplit(ary, indices_or_sections):
    """Split an array into multiple sub-arrays horizontally (column-wise).

    Please refer to the `split` documentation.  `hsplit` is equivalent
    to `split` with ``axis=1``, the array is always split along the second
    axis regardless of the array dimension.

    See Also
    --------
    split : Split an array into multiple sub-arrays of equal size.

    Examples
    --------
    >>> x = np.arange(16.0).reshape(4, 4)
    >>> x
    array([[  0.,   1.,   2.,   3.],
           [  4.,   5.,   6.,   7.],
           [  8.,   9.,  10.,  11.],
           [ 12.,  13.,  14.,  15.]])
    >>> np.hsplit(x, 2)
    [array([[  0.,   1.],
           [  4.,   5.],
           [  8.,   9.],
           [ 12.,  13.]]),
     array([[  2.,   3.],
           [  6.,   7.],
           [ 10.,  11.],
           [ 14.,  15.]])]
    >>> np.hsplit(x, np.array([3, 6]))
    [array([[  0.,   1.,   2.],
           [  4.,   5.,   6.],
           [  8.,   9.,  10.],
           [ 12.,  13.,  14.]]),
     array([[  3.],
           [  7.],
           [ 11.],
           [ 15.]]),
     array([], dtype=float64)]

    With a higher dimensional array the split is still along the second axis.

    >>> x = np.arange(8.0).reshape(2, 2, 2)
    >>> x
    array([[[ 0.,  1.],
            [ 2.,  3.]],
           [[ 4.,  5.],
            [ 6.,  7.]]])
    >>> np.hsplit(x, 2)
    [array([[[ 0.,  1.]],
           [[ 4.,  5.]]]),
     array([[[ 2.,  3.]],
           [[ 6.,  7.]]])]
    
    """
    raise NotImplementedError


def hstack(tup):
    """Stack arrays in sequence horizontally (column wise).

    Take a sequence of arrays and stack them horizontally to make
    a single array. Rebuild arrays divided by `hsplit`.

    Parameters
    ----------
    tup : sequence of ndarrays
        All arrays must have the same shape along all but the second axis.

    Returns
    -------
    stacked : ndarray
        The array formed by stacking the given arrays.

    See Also
    --------
    vstack : Stack arrays in sequence vertically (row wise).
    dstack : Stack arrays in sequence depth wise (along third axis).
    concatenate : Join a sequence of arrays together.
    hsplit : Split array along second axis.

    Notes
    -----
    Equivalent to ``np.concatenate(tup, axis=1)``

    Examples
    --------
    >>> a = np.array((1,2,3))
    >>> b = np.array((2,3,4))
    >>> np.hstack((a,b))
    array([1, 2, 3, 2, 3, 4])
    >>> a = np.array([[1],[2],[3]])
    >>> b = np.array([[2],[3],[4]])
    >>> np.hstack((a,b))
    array([[1, 2],
           [2, 3],
           [3, 4]])
    
    """
    raise NotImplementedError


def i0(x):
    """Modified Bessel function of the first kind, order 0.

    Usually denoted :math:`I_0`.  This function does broadcast, but will *not*
    "up-cast" int dtype arguments unless accompanied by at least one float or
    complex dtype argument (see Raises below).

    Parameters
    ----------
    x : array_like, dtype float or complex
        Argument of the Bessel function.

    Returns
    -------
    out : ndarray, shape = x.shape, dtype = x.dtype
        The modified Bessel function evaluated at each of the elements of `x`.

    Raises
    ------
    TypeError: array cannot be safely cast to required type
        If argument consists exclusively of int dtypes.

    See Also
    --------
    scipy.special.iv, scipy.special.ive

    Notes
    -----
    We use the algorithm published by Clenshaw [1]_ and referenced by
    Abramowitz and Stegun [2]_, for which the function domain is partitioned
    into the two intervals [0,8] and (8,inf), and Chebyshev polynomial
    expansions are employed in each interval. Relative error on the domain
    [0,30] using IEEE arithmetic is documented [3]_ as having a peak of 5.8e-16
    with an rms of 1.4e-16 (n = 30000).

    References
    ----------
    .. [1] C. W. Clenshaw, "Chebyshev series for mathematical functions," in
           *National Physical Laboratory Mathematical Tables*, vol. 5, London:
           Her Majesty's Stationery Office, 1962.
    .. [2] M. Abramowitz and I. A. Stegun, *Handbook of Mathematical
           Functions*, 10th printing, New York: Dover, 1964, pp. 379.
           http://www.math.sfu.ca/~cbm/aands/page_379.htm
    .. [3] http://kobesearch.cpan.org/htdocs/Math-Cephes/Math/Cephes.html

    Examples
    --------
    >>> np.i0([0.])
    array(1.0)
    >>> np.i0([0., 1. + 2j])
    array([ 1.00000000+0.j        ,  0.18785373+0.64616944j])
    
    """
    raise NotImplementedError


def imag(val):
    """Return the imaginary part of the elements of the array.

    Parameters
    ----------
    val : array_like
        Input array.

    Returns
    -------
    out : ndarray
        Output array. If `val` is real, the type of `val` is used for the
        output.  If `val` has complex elements, the returned type is float.

    See Also
    --------
    real, angle, real_if_close

    Examples
    --------
    >>> a = np.array([1+2j, 3+4j, 5+6j])
    >>> a.imag
    array([ 2.,  4.,  6.])
    >>> a.imag = np.array([8, 10, 12])
    >>> a
    array([ 1. +8.j,  3.+10.j,  5.+12.j])
    
    """
    raise NotImplementedError


def in1d(ar1, ar2, assume_unique=False):
    """Test whether each element of a 1D array is also present in a second array.

    Returns a boolean array the same length as `ar1` that is True
    where an element of `ar1` is in `ar2` and False otherwise.

    Parameters
    ----------
    ar1 : array_like, shape (M,)
        Input array.
    ar2 : array_like
        The values against which to test each value of `ar1`.
    assume_unique : bool, optional
        If True, the input arrays are both assumed to be unique, which
        can speed up the calculation.  Default is False.

    Returns
    -------
    mask : ndarray of bools, shape(M,)
        The values `ar1[mask]` are in `ar2`.

    See Also
    --------
    numpy.lib.arraysetops : Module with a number of other functions for
                            performing set operations on arrays.

    Notes
    -----
    `in1d` can be considered as an element-wise function version of the
    python keyword `in`, for 1D sequences. ``in1d(a, b)`` is roughly
    equivalent to ``np.array([item in b for item in a])``.

    .. versionadded:: 1.4.0

    Examples
    --------
    >>> test = np.array([0, 1, 2, 5, 0])
    >>> states = [0, 2]
    >>> mask = np.in1d(test, states)
    >>> mask
    array([ True, False,  True, False,  True], dtype=bool)
    >>> test[mask]
    array([0, 2, 0])
    
    """
    raise NotImplementedError


def indices(dimensions, dtype=int):
    """Return an array representing the indices of a grid.

    Compute an array where the subarrays contain index values 0,1,...
    varying only along the corresponding axis.

    Parameters
    ----------
    dimensions : sequence of ints
        The shape of the grid.
    dtype : dtype, optional
        Data type of the result.

    Returns
    -------
    grid : ndarray
        The array of grid indices,
        ``grid.shape = (len(dimensions),) + tuple(dimensions)``.

    See Also
    --------
    mgrid, meshgrid

    Notes
    -----
    The output shape is obtained by prepending the number of dimensions
    in front of the tuple of dimensions, i.e. if `dimensions` is a tuple
    ``(r0, ..., rN-1)`` of length ``N``, the output shape is
    ``(N,r0,...,rN-1)``.

    The subarrays ``grid[k]`` contains the N-D array of indices along the
    ``k-th`` axis. Explicitly::

        grid[k,i0,i1,...,iN-1] = ik

    Examples
    --------
    >>> grid = np.indices((2, 3))
    >>> grid.shape
    (2, 2, 3)
    >>> grid[0]        # row indices
    array([[0, 0, 0],
           [1, 1, 1]])
    >>> grid[1]        # column indices
    array([[0, 1, 2],
           [0, 1, 2]])

    The indices can be used as an index into an array.

    >>> x = np.arange(20).reshape(5, 4)
    >>> row, col = np.indices((2, 3))
    >>> x[row, col]
    array([[0, 1, 2],
           [4, 5, 6]])

    Note that it would be more straightforward in the above example to
    extract the required elements directly with ``x[:2, :3]``.
    
    """
    raise NotImplementedError


def info(object=None, maxwidth=76, output=sys.stdout, toplevel='numpy'):
    """Get help information for a function, class, or module.

    Parameters
    ----------
    object : object or str, optional
        Input object or name to get information about. If `object` is a
        numpy object, its docstring is given. If it is a string, available
        modules are searched for matching objects.
        If None, information about `info` itself is returned.
    maxwidth : int, optional
        Printing width.
    output : file like object, optional
        File like object that the output is written to, default is ``stdout``.
        The object has to be opened in 'w' or 'a' mode.
    toplevel : str, optional
        Start search at this level.

    See Also
    --------
    source, lookfor

    Notes
    -----
    When used interactively with an object, ``np.info(obj)`` is equivalent to
    ``help(obj)`` on the Python prompt or ``obj?`` on the IPython prompt.

    Examples
    --------
    >>> np.info(np.polyval) # doctest: +SKIP
       polyval(p, x)
         Evaluate the polynomial p at x.
         ...

    When using a string for `object` it is possible to get multiple results.

    >>> np.info('fft') # doctest: +SKIP
         *** Found in numpy ***
    Core FFT routines
    ...
         *** Found in numpy.fft ***
     fft(a, n=None, axis=-1)
    ...
         *** Repeat reference found in numpy.fft.fftpack ***
         *** Total of 3 references found. ***
    
    """
    raise NotImplementedError


def inner(a, b):
    """inner(a, b)

    Inner product of two arrays.

    Ordinary inner product of vectors for 1-D arrays (without complex
    conjugation), in higher dimensions a sum product over the last axes.

    Parameters
    ----------
    a, b : array_like
        If `a` and `b` are nonscalar, their last dimensions of must match.

    Returns
    -------
    out : ndarray
        `out.shape = a.shape[:-1] + b.shape[:-1]`

    Raises
    ------
    ValueError
        If the last dimension of `a` and `b` has different size.

    See Also
    --------
    tensordot : Sum products over arbitrary axes.
    dot : Generalised matrix product, using second last dimension of `b`.
    einsum : Einstein summation convention.

    Notes
    -----
    For vectors (1-D arrays) it computes the ordinary inner-product::

        np.inner(a, b) = sum(a[:]*b[:])

    More generally, if `ndim(a) = r > 0` and `ndim(b) = s > 0`::

        np.inner(a, b) = np.tensordot(a, b, axes=(-1,-1))

    or explicitly::

        np.inner(a, b)[i0,...,ir-1,j0,...,js-1]
             = sum(a[i0,...,ir-1,:]*b[j0,...,js-1,:])

    In addition `a` or `b` may be scalars, in which case::

       np.inner(a,b) = a*b

    Examples
    --------
    Ordinary inner product for vectors:

    >>> a = np.array([1,2,3])
    >>> b = np.array([0,1,0])
    >>> np.inner(a, b)
    2

    A multidimensional example:

    >>> a = np.arange(24).reshape((2,3,4))
    >>> b = np.arange(4)
    >>> np.inner(a, b)
    array([[ 14,  38,  62],
           [ 86, 110, 134]])

    An example where `b` is a scalar:

    >>> np.inner(np.eye(2), 7)
    array([[ 7.,  0.],
           [ 0.,  7.]])
    
    """
    # BUILTIN
    raise NotImplementedError


def insert(arr, obj, values, axis=None):
    """Insert values along the given axis before the given indices.

    Parameters
    ----------
    arr : array_like
        Input array.
    obj : int, slice or sequence of ints
        Object that defines the index or indices before which `values` is
        inserted.
    values : array_like
        Values to insert into `arr`. If the type of `values` is different
        from that of `arr`, `values` is converted to the type of `arr`.
    axis : int, optional
        Axis along which to insert `values`.  If `axis` is None then `arr`
        is flattened first.

    Returns
    -------
    out : ndarray
        A copy of `arr` with `values` inserted.  Note that `insert`
        does not occur in-place: a new array is returned. If
        `axis` is None, `out` is a flattened array.

    See Also
    --------
    append : Append elements at the end of an array.
    delete : Delete elements from an array.

    Examples
    --------
    >>> a = np.array([[1, 1], [2, 2], [3, 3]])
    >>> a
    array([[1, 1],
           [2, 2],
           [3, 3]])
    >>> np.insert(a, 1, 5)
    array([1, 5, 1, 2, 2, 3, 3])
    >>> np.insert(a, 1, 5, axis=1)
    array([[1, 5, 1],
           [2, 5, 2],
           [3, 5, 3]])

    >>> b = a.flatten()
    >>> b
    array([1, 1, 2, 2, 3, 3])
    >>> np.insert(b, [2, 2], [5, 6])
    array([1, 1, 5, 6, 2, 2, 3, 3])

    >>> np.insert(b, slice(2, 4), [5, 6])
    array([1, 1, 5, 2, 6, 2, 3, 3])

    >>> np.insert(b, [2, 2], [7.13, False]) # type casting
    array([1, 1, 7, 0, 2, 2, 3, 3])

    >>> x = np.arange(8).reshape(2, 4)
    >>> idx = (1, 3)
    >>> np.insert(x, idx, 999, axis=1)
    array([[  0, 999,   1,   2, 999,   3],
           [  4, 999,   5,   6, 999,   7]])
    
    """
    raise NotImplementedError


def int_asbuffer():
    """
    
    """
    # BUILTIN
    raise NotImplementedError


def interp(x, xp, fp, left=None, right=None):
    """One-dimensional linear interpolation.

    Returns the one-dimensional piecewise linear interpolant to a function
    with given values at discrete data-points.

    Parameters
    ----------
    x : array_like
        The x-coordinates of the interpolated values.

    xp : 1-D sequence of floats
        The x-coordinates of the data points, must be increasing.

    fp : 1-D sequence of floats
        The y-coordinates of the data points, same length as `xp`.

    left : float, optional
        Value to return for `x < xp[0]`, default is `fp[0]`.

    right : float, optional
        Value to return for `x > xp[-1]`, defaults is `fp[-1]`.

    Returns
    -------
    y : {float, ndarray}
        The interpolated values, same shape as `x`.

    Raises
    ------
    ValueError
        If `xp` and `fp` have different length

    Notes
    -----
    Does not check that the x-coordinate sequence `xp` is increasing.
    If `xp` is not increasing, the results are nonsense.
    A simple check for increasingness is::

        np.all(np.diff(xp) > 0)


    Examples
    --------
    >>> xp = [1, 2, 3]
    >>> fp = [3, 2, 0]
    >>> np.interp(2.5, xp, fp)
    1.0
    >>> np.interp([0, 1, 1.5, 2.72, 3.14], xp, fp)
    array([ 3. ,  3. ,  2.5 ,  0.56,  0. ])
    >>> UNDEF = -99.0
    >>> np.interp(3.14, xp, fp, right=UNDEF)
    -99.0

    Plot an interpolant to the sine function:

    >>> x = np.linspace(0, 2*np.pi, 10)
    >>> y = np.sin(x)
    >>> xvals = np.linspace(0, 2*np.pi, 50)
    >>> yinterp = np.interp(xvals, x, y)
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x, y, 'o')
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.plot(xvals, yinterp, '-x')
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.show()
    
    """
    raise NotImplementedError


def intersect1d(ar1, ar2, assume_unique=False):
    """Find the intersection of two arrays.

    Return the sorted, unique values that are in both of the input arrays.

    Parameters
    ----------
    ar1, ar2 : array_like
        Input arrays.
    assume_unique : bool
        If True, the input arrays are both assumed to be unique, which
        can speed up the calculation.  Default is False.

    Returns
    -------
    out : ndarray
        Sorted 1D array of common and unique elements.

    See Also
    --------
    numpy.lib.arraysetops : Module with a number of other functions for
                            performing set operations on arrays.

    Examples
    --------
    >>> np.intersect1d([1, 3, 4, 3], [3, 1, 2, 1])
    array([1, 3])
    
    """
    raise NotImplementedError


def ipmt(rate, per, nper, pv, fv=0.0, when='end'):
    """Not implemented. Compute the payment portion for loan interest.

    Parameters
    ----------
    rate : scalar or array_like of shape(M, )
        Rate of interest as decimal (not per cent) per period
    per : scalar or array_like of shape(M, )
        Interest paid against the loan changes during the life or the loan.
        The `per` is the payment period to calculate the interest amount.
    nper : scalar or array_like of shape(M, )
        Number of compounding periods
    pv : scalar or array_like of shape(M, )
        Present value
    fv : scalar or array_like of shape(M, ), optional
        Future value
    when : {{'begin', 1}, {'end', 0}}, {string, int}, optional
        When payments are due ('begin' (1) or 'end' (0)).
        Defaults to {'end', 0}.

    Returns
    -------
    out : ndarray
        Interest portion of payment.  If all input is scalar, returns a scalar
        float.  If any input is array_like, returns interest payment for each
        input element. If multiple inputs are array_like, they all must have
        the same shape.

    See Also
    --------
    ppmt, pmt, pv

    Notes
    -----
    The total payment is made up of payment against principal plus interest.

    ``pmt = ppmt + ipmt``
    
    """
    raise NotImplementedError


def irr(values):
    """Return the Internal Rate of Return (IRR).

    This is the "average" periodically compounded rate of return
    that gives a net present value of 0.0; for a more complete explanation,
    see Notes below.

    Parameters
    ----------
    values : array_like, shape(N,)
        Input cash flows per time period.  By convention, net "deposits"
        are negative and net "withdrawals" are positive.  Thus, for example,
        at least the first element of `values`, which represents the initial
        investment, will typically be negative.

    Returns
    -------
    out : float
        Internal Rate of Return for periodic input values.

    Notes
    -----
    The IRR is perhaps best understood through an example (illustrated
    using np.irr in the Examples section below).  Suppose one invests
    100 units and then makes the following withdrawals at regular
    (fixed) intervals: 39, 59, 55, 20.  Assuming the ending value is 0,
    one's 100 unit investment yields 173 units; however, due to the
    combination of compounding and the periodic withdrawals, the
    "average" rate of return is neither simply 0.73/4 nor (1.73)^0.25-1.
    Rather, it is the solution (for :math:`r`) of the equation:

    .. math:: -100 + \frac{39}{1+r} + \frac{59}{(1+r)^2}
     + \frac{55}{(1+r)^3} + \frac{20}{(1+r)^4} = 0

    In general, for `values` :math:`= [v_0, v_1, ... v_M]`,
    irr is the solution of the equation: [G]_

    .. math:: \sum_{t=0}^M{\frac{v_t}{(1+irr)^{t}}} = 0

    References
    ----------
    .. [G] L. J. Gitman, "Principles of Managerial Finance, Brief," 3rd ed.,
       Addison-Wesley, 2003, pg. 348.

    Examples
    --------
    >>> np.irr([-100, 39, 59, 55, 20])
    0.2809484211599611

    (Compare with the Example given for numpy.lib.financial.npv)
    
    """
    raise NotImplementedError


def iscomplex(x):
    """Returns a bool array, where True if input element is complex.

    What is tested is whether the input has a non-zero imaginary part, not if
    the input type is complex.

    Parameters
    ----------
    x : array_like
        Input array.

    Returns
    -------
    out : ndarray of bools
        Output array.

    See Also
    --------
    isreal
    iscomplexobj : Return True if x is a complex type or an array of complex
                   numbers.

    Examples
    --------
    >>> np.iscomplex([1+1j, 1+0j, 4.5, 3, 2, 2j])
    array([ True, False, False, False, False,  True], dtype=bool)
    
    """
    raise NotImplementedError


def iscomplexobj(x):
    """Return True if x is a complex type or an array of complex numbers.

    The type of the input is checked, not the value. So even if the input
    has an imaginary part equal to zero, `iscomplexobj` evaluates to True
    if the data type is complex.

    Parameters
    ----------
    x : any
        The input can be of any type and shape.

    Returns
    -------
    y : bool
        The return value, True if `x` is of a complex type.

    See Also
    --------
    isrealobj, iscomplex

    Examples
    --------
    >>> np.iscomplexobj(1)
    False
    >>> np.iscomplexobj(1+0j)
    True
    >>> np.iscomplexobj([3, 1+0j, True])
    True
    
    """
    raise NotImplementedError


def isfortran(a):
    """Returns True if array is arranged in Fortran-order in memory
    and dimension > 1.

    Parameters
    ----------
    a : ndarray
        Input array.


    Examples
    --------

    np.array allows to specify whether the array is written in C-contiguous
    order (last index varies the fastest), or FORTRAN-contiguous order in
    memory (first index varies the fastest).

    >>> a = np.array([[1, 2, 3], [4, 5, 6]], order='C')
    >>> a
    array([[1, 2, 3],
           [4, 5, 6]])
    >>> np.isfortran(a)
    False

    >>> b = np.array([[1, 2, 3], [4, 5, 6]], order='FORTRAN')
    >>> b
    array([[1, 2, 3],
           [4, 5, 6]])
    >>> np.isfortran(b)
    True


    The transpose of a C-ordered array is a FORTRAN-ordered array.

    >>> a = np.array([[1, 2, 3], [4, 5, 6]], order='C')
    >>> a
    array([[1, 2, 3],
           [4, 5, 6]])
    >>> np.isfortran(a)
    False
    >>> b = a.T
    >>> b
    array([[1, 4],
           [2, 5],
           [3, 6]])
    >>> np.isfortran(b)
    True

    1-D arrays always evaluate as False.

    >>> np.isfortran(np.array([1, 2], order='FORTRAN'))
    False
    
    """
    raise NotImplementedError


def isneginf(x, y=None):
    """Test element-wise for negative infinity, return result as bool array.

    Parameters
    ----------
    x : array_like
        The input array.
    y : array_like, optional
        A boolean array with the same shape and type as `x` to store the
        result.

    Returns
    -------
    y : ndarray
        A boolean array with the same dimensions as the input.
        If second argument is not supplied then a numpy boolean array is
        returned with values True where the corresponding element of the
        input is negative infinity and values False where the element of
        the input is not negative infinity.

        If a second argument is supplied the result is stored there. If the
        type of that array is a numeric type the result is represented as
        zeros and ones, if the type is boolean then as False and True. The
        return value `y` is then a reference to that array.

    See Also
    --------
    isinf, isposinf, isnan, isfinite

    Notes
    -----
    Numpy uses the IEEE Standard for Binary Floating-Point for Arithmetic
    (IEEE 754).

    Errors result if the second argument is also supplied when x is a scalar
    input, or if first and second arguments have different shapes.

    Examples
    --------
    >>> np.isneginf(np.NINF)
    array(True, dtype=bool)
    >>> np.isneginf(np.inf)
    array(False, dtype=bool)
    >>> np.isneginf(np.PINF)
    array(False, dtype=bool)
    >>> np.isneginf([-np.inf, 0., np.inf])
    array([ True, False, False], dtype=bool)

    >>> x = np.array([-np.inf, 0., np.inf])
    >>> y = np.array([2, 2, 2])
    >>> np.isneginf(x, y)
    array([1, 0, 0])
    >>> y
    array([1, 0, 0])
    
    """
    raise NotImplementedError


def isposinf(x, y=None):
    """Test element-wise for positive infinity, return result as bool array.

    Parameters
    ----------
    x : array_like
        The input array.
    y : array_like, optional
        A boolean array with the same shape as `x` to store the result.

    Returns
    -------
    y : ndarray
        A boolean array with the same dimensions as the input.
        If second argument is not supplied then a boolean array is returned
        with values True where the corresponding element of the input is
        positive infinity and values False where the element of the input is
        not positive infinity.

        If a second argument is supplied the result is stored there. If the
        type of that array is a numeric type the result is represented as zeros
        and ones, if the type is boolean then as False and True.
        The return value `y` is then a reference to that array.

    See Also
    --------
    isinf, isneginf, isfinite, isnan

    Notes
    -----
    Numpy uses the IEEE Standard for Binary Floating-Point for Arithmetic
    (IEEE 754).

    Errors result if the second argument is also supplied when `x` is a
    scalar input, or if first and second arguments have different shapes.

    Examples
    --------
    >>> np.isposinf(np.PINF)
    array(True, dtype=bool)
    >>> np.isposinf(np.inf)
    array(True, dtype=bool)
    >>> np.isposinf(np.NINF)
    array(False, dtype=bool)
    >>> np.isposinf([-np.inf, 0., np.inf])
    array([False, False,  True], dtype=bool)

    >>> x = np.array([-np.inf, 0., np.inf])
    >>> y = np.array([2, 2, 2])
    >>> np.isposinf(x, y)
    array([0, 0, 1])
    >>> y
    array([0, 0, 1])
    
    """
    raise NotImplementedError


def isreal(x):
    """Returns a bool array, where True if input element is real.

    If element has complex type with zero complex part, the return value
    for that element is True.

    Parameters
    ----------
    x : array_like
        Input array.

    Returns
    -------
    out : ndarray, bool
        Boolean array of same shape as `x`.

    See Also
    --------
    iscomplex
    isrealobj : Return True if x is not a complex type.

    Examples
    --------
    >>> np.isreal([1+1j, 1+0j, 4.5, 3, 2, 2j])
    array([False,  True,  True,  True,  True, False], dtype=bool)
    
    """
    raise NotImplementedError


def isrealobj(x):
    """Return True if x is a not complex type or an array of complex numbers.

    The type of the input is checked, not the value. So even if the input
    has an imaginary part equal to zero, `isrealobj` evaluates to False
    if the data type is complex.

    Parameters
    ----------
    x : any
        The input can be of any type and shape.

    Returns
    -------
    y : bool
        The return value, False if `x` is of a complex type.

    See Also
    --------
    iscomplexobj, isreal

    Examples
    --------
    >>> np.isrealobj(1)
    True
    >>> np.isrealobj(1+0j)
    False
    >>> np.isrealobj([3, 1+0j, True])
    False
    
    """
    raise NotImplementedError


def isscalar(num):
    """Returns True if the type of `num` is a scalar type.

    Parameters
    ----------
    num : any
        Input argument, can be of any type and shape.

    Returns
    -------
    val : bool
        True if `num` is a scalar type, False if it is not.

    Examples
    --------
    >>> np.isscalar(3.1)
    True
    >>> np.isscalar([3.1])
    False
    >>> np.isscalar(False)
    True
    
    """
    raise NotImplementedError


def issctype(rep):
    """Determines whether the given object represents a scalar data-type.

    Parameters
    ----------
    rep : any
        If `rep` is an instance of a scalar dtype, True is returned. If not,
        False is returned.

    Returns
    -------
    out : bool
        Boolean result of check whether `rep` is a scalar dtype.

    See Also
    --------
    issubsctype, issubdtype, obj2sctype, sctype2char

    Examples
    --------
    >>> np.issctype(np.int32)
    True
    >>> np.issctype(list)
    False
    >>> np.issctype(1.1)
    False

    Strings are also a scalar type:

    >>> np.issctype(np.dtype('str'))
    True
    
    """
    raise NotImplementedError


def issubclass_(arg1, arg2):
    """Determine if a class is a subclass of a second class.

    `issubclass_` is equivalent to the Python built-in ``issubclass``,
    except that it returns False instead of raising a TypeError is one
    of the arguments is not a class.

    Parameters
    ----------
    arg1 : class
        Input class. True is returned if `arg1` is a subclass of `arg2`.
    arg2 : class or tuple of classes.
        Input class. If a tuple of classes, True is returned if `arg1` is a
        subclass of any of the tuple elements.

    Returns
    -------
    out : bool
        Whether `arg1` is a subclass of `arg2` or not.

    See Also
    --------
    issubsctype, issubdtype, issctype

    Examples
    --------
    >>> np.issubclass_(np.int32, np.int)
    True
    >>> np.issubclass_(np.int32, np.float)
    False
    
    """
    raise NotImplementedError


def issubdtype(arg1, arg2):
    """Returns True if first argument is a typecode lower/equal in type hierarchy.

    Parameters
    ----------
    arg1, arg2 : dtype_like
        dtype or string representing a typecode.

    Returns
    -------
    out : bool

    See Also
    --------
    issubsctype, issubclass_
    numpy.core.numerictypes : Overview of numpy type hierarchy.

    Examples
    --------
    >>> np.issubdtype('S1', str)
    True
    >>> np.issubdtype(np.float64, np.float32)
    False
    
    """
    raise NotImplementedError


def issubsctype(arg1, arg2):
    """Determine if the first argument is a subclass of the second argument.

    Parameters
    ----------
    arg1, arg2 : dtype or dtype specifier
        Data-types.

    Returns
    -------
    out : bool
        The result.

    See Also
    --------
    issctype, issubdtype,obj2sctype

    Examples
    --------
    >>> np.issubsctype('S8', str)
    True
    >>> np.issubsctype(np.array([1]), np.int)
    True
    >>> np.issubsctype(np.array([1]), np.float)
    False
    
    """
    raise NotImplementedError


def iterable(y):
    """Check whether or not an object can be iterated over.

    Parameters
    ----------
    y : object
      Input object.

    Returns
    -------
    b : {0, 1}
      Return 1 if the object has an iterator method or is a sequence,
      and 0 otherwise.


    Examples
    --------
    >>> np.iterable([1, 2, 3])
    1
    >>> np.iterable(2)
    0
    
    """
    raise NotImplementedError


def ix_():
    """Construct an open mesh from multiple sequences.

    This function takes N 1-D sequences and returns N outputs with N
    dimensions each, such that the shape is 1 in all but one dimension
    and the dimension with the non-unit shape value cycles through all
    N dimensions.

    Using `ix_` one can quickly construct index arrays that will index
    the cross product. ``a[np.ix_([1,3],[2,5])]`` returns the array
    ``[[a[1,2] a[1,5]], [a[3,2] a[3,5]]]``.

    Parameters
    ----------
    args : 1-D sequences

    Returns
    -------
    out : tuple of ndarrays
        N arrays with N dimensions each, with N the number of input
        sequences. Together these arrays form an open mesh.

    See Also
    --------
    ogrid, mgrid, meshgrid

    Examples
    --------
    >>> a = np.arange(10).reshape(2, 5)
    >>> a
    array([[0, 1, 2, 3, 4],
           [5, 6, 7, 8, 9]])
    >>> ixgrid = np.ix_([0,1], [2,4])
    >>> ixgrid
    (array([[0],
           [1]]), array([[2, 4]]))
    >>> ixgrid[0].shape, ixgrid[1].shape
    ((2, 1), (1, 2))
    >>> a[ixgrid]
    array([[2, 4],
           [7, 9]])
    
    """
    raise NotImplementedError


def kaiser(M, beta):
    """Return the Kaiser window.

    The Kaiser window is a taper formed by using a Bessel function.

    Parameters
    ----------
    M : int
        Number of points in the output window. If zero or less, an
        empty array is returned.
    beta : float
        Shape parameter for window.

    Returns
    -------
    out : array
        The window, normalized to one (the value one
        appears only if the number of samples is odd).

    See Also
    --------
    bartlett, blackman, hamming, hanning

    Notes
    -----
    The Kaiser window is defined as

    .. math::  w(n) = I_0\left( \beta \sqrt{1-\frac{4n^2}{(M-1)^2}}
               \right)/I_0(\beta)

    with

    .. math:: \quad -\frac{M-1}{2} \leq n \leq \frac{M-1}{2},

    where :math:`I_0` is the modified zeroth-order Bessel function.

    The Kaiser was named for Jim Kaiser, who discovered a simple approximation
    to the DPSS window based on Bessel functions.
    The Kaiser window is a very good approximation to the Digital Prolate
    Spheroidal Sequence, or Slepian window, which is the transform which
    maximizes the energy in the main lobe of the window relative to total
    energy.

    The Kaiser can approximate many other windows by varying the beta
    parameter.

    ====  =======================
    beta  Window shape
    ====  =======================
    0     Rectangular
    5     Similar to a Hamming
    6     Similar to a Hanning
    8.6   Similar to a Blackman
    ====  =======================

    A beta value of 14 is probably a good starting point. Note that as beta
    gets large, the window narrows, and so the number of samples needs to be
    large enough to sample the increasingly narrow spike, otherwise nans will
    get returned.


    Most references to the Kaiser window come from the signal processing
    literature, where it is used as one of many windowing functions for
    smoothing values.  It is also known as an apodization (which means
    "removing the foot", i.e. smoothing discontinuities at the beginning
    and end of the sampled signal) or tapering function.

    References
    ----------
    .. [1] J. F. Kaiser, "Digital Filters" - Ch 7 in "Systems analysis by
           digital computer", Editors: F.F. Kuo and J.F. Kaiser, p 218-285.
           John Wiley and Sons, New York, (1966).
    .. [2] E.R. Kanasewich, "Time Sequence Analysis in Geophysics", The
           University of Alberta Press, 1975, pp. 177-178.
    .. [3] Wikipedia, "Window function",
           http://en.wikipedia.org/wiki/Window_function

    Examples
    --------
    >>> from numpy import kaiser
    >>> kaiser(12, 14)
    array([  7.72686684e-06,   3.46009194e-03,   4.65200189e-02,
             2.29737120e-01,   5.99885316e-01,   9.45674898e-01,
             9.45674898e-01,   5.99885316e-01,   2.29737120e-01,
             4.65200189e-02,   3.46009194e-03,   7.72686684e-06])


    Plot the window and the frequency response:

    >>> from numpy import clip, log10, array, kaiser, linspace
    >>> from numpy.fft import fft, fftshift
    >>> import matplotlib.pyplot as plt

    >>> window = kaiser(51, 14)
    >>> plt.plot(window)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Kaiser window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Sample")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    >>> plt.figure()
    <matplotlib.figure.Figure object at 0x...>
    >>> A = fft(window, 2048) / 25.5
    >>> mag = abs(fftshift(A))
    >>> freq = linspace(-0.5,0.5,len(A))
    >>> response = 20*log10(mag)
    >>> response = clip(response,-100,100)
    >>> plt.plot(freq, response)
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Frequency response of Kaiser window")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Magnitude [dB]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("Normalized frequency [cycles per sample]")
    <matplotlib.text.Text object at 0x...>
    >>> plt.axis('tight')
    (-0.5, 0.5, -100.0, ...)
    >>> plt.show()
    
    """
    raise NotImplementedError


def kron(a, b):
    """Kronecker product of two arrays.

    Computes the Kronecker product, a composite array made of blocks of the
    second array scaled by the first.

    Parameters
    ----------
    a, b : array_like

    Returns
    -------
    out : ndarray

    See Also
    --------

    outer : The outer product

    Notes
    -----

    The function assumes that the number of dimenensions of `a` and `b`
    are the same, if necessary prepending the smallest with ones.
    If `a.shape = (r0,r1,..,rN)` and `b.shape = (s0,s1,...,sN)`,
    the Kronecker product has shape `(r0*s0, r1*s1, ..., rN*SN)`.
    The elements are products of elements from `a` and `b`, organized
    explicitly by::

        kron(a,b)[k0,k1,...,kN] = a[i0,i1,...,iN] * b[j0,j1,...,jN]

    where::

        kt = it * st + jt,  t = 0,...,N

    In the common 2-D case (N=1), the block structure can be visualized::

        [[ a[0,0]*b,   a[0,1]*b,  ... , a[0,-1]*b  ],
         [  ...                              ...   ],
         [ a[-1,0]*b,  a[-1,1]*b, ... , a[-1,-1]*b ]]


    Examples
    --------
    >>> np.kron([1,10,100], [5,6,7])
    array([  5,   6,   7,  50,  60,  70, 500, 600, 700])
    >>> np.kron([5,6,7], [1,10,100])
    array([  5,  50, 500,   6,  60, 600,   7,  70, 700])

    >>> np.kron(np.eye(2), np.ones((2,2)))
    array([[ 1.,  1.,  0.,  0.],
           [ 1.,  1.,  0.,  0.],
           [ 0.,  0.,  1.,  1.],
           [ 0.,  0.,  1.,  1.]])

    >>> a = np.arange(100).reshape((2,5,2,5))
    >>> b = np.arange(24).reshape((2,3,4))
    >>> c = np.kron(a,b)
    >>> c.shape
    (2, 10, 6, 20)
    >>> I = (1,3,0,2)
    >>> J = (0,2,1)
    >>> J1 = (0,) + J             # extend to ndim=4
    >>> S1 = (1,) + b.shape
    >>> K = tuple(np.array(I) * np.array(S1) + np.array(J1))
    >>> c[K] == a[I]*b[J]
    True
    
    """
    raise NotImplementedError


def lexsort(keys, axis=-1):
    """lexsort(keys, axis=-1)

    Perform an indirect sort using a sequence of keys.

    Given multiple sorting keys, which can be interpreted as columns in a
    spreadsheet, lexsort returns an array of integer indices that describes
    the sort order by multiple columns. The last key in the sequence is used
    for the primary sort order, the second-to-last key for the secondary sort
    order, and so on. The keys argument must be a sequence of objects that
    can be converted to arrays of the same shape. If a 2D array is provided
    for the keys argument, it's rows are interpreted as the sorting keys and
    sorting is according to the last row, second last row etc.

    Parameters
    ----------
    keys : (k,N) array or tuple containing k (N,)-shaped sequences
        The `k` different "columns" to be sorted.  The last column (or row if
        `keys` is a 2D array) is the primary sort key.
    axis : int, optional
        Axis to be indirectly sorted.  By default, sort over the last axis.

    Returns
    -------
    indices : (N,) ndarray of ints
        Array of indices that sort the keys along the specified axis.

    See Also
    --------
    argsort : Indirect sort.
    ndarray.sort : In-place sort.
    sort : Return a sorted copy of an array.

    Examples
    --------
    Sort names: first by surname, then by name.

    >>> surnames =    ('Hertz',    'Galilei', 'Hertz')
    >>> first_names = ('Heinrich', 'Galileo', 'Gustav')
    >>> ind = np.lexsort((first_names, surnames))
    >>> ind
    array([1, 2, 0])

    >>> [surnames[i] + ", " + first_names[i] for i in ind]
    ['Galilei, Galileo', 'Hertz, Gustav', 'Hertz, Heinrich']

    Sort two columns of numbers:

    >>> a = [1,5,1,4,3,4,4] # First column
    >>> b = [9,4,0,4,0,2,1] # Second column
    >>> ind = np.lexsort((b,a)) # Sort by a, then by b
    >>> print ind
    [2 0 4 6 5 3 1]

    >>> [(a[i],b[i]) for i in ind]
    [(1, 0), (1, 9), (3, 0), (4, 1), (4, 2), (4, 4), (5, 4)]

    Note that sorting is first according to the elements of ``a``.
    Secondary sorting is according to the elements of ``b``.

    A normal ``argsort`` would have yielded:

    >>> [(a[i],b[i]) for i in np.argsort(a)]
    [(1, 9), (1, 0), (3, 0), (4, 4), (4, 2), (4, 1), (5, 4)]

    Structured arrays are sorted lexically by ``argsort``:

    >>> x = np.array([(1,9), (5,4), (1,0), (4,4), (3,0), (4,2), (4,1)],
    ...              dtype=np.dtype([('x', int), ('y', int)]))

    >>> np.argsort(x) # or np.argsort(x, order=('x', 'y'))
    array([2, 0, 4, 6, 5, 3, 1])
    
    """
    # BUILTIN
    raise NotImplementedError


def load(file, mmap_mode=None):
    """Load a pickled, ``.npy``, or ``.npz`` binary file.

    Parameters
    ----------
    file : file-like object or string
        The file to read.  It must support ``seek()`` and ``read()`` methods.
        If the filename extension is ``.gz``, the file is first decompressed.
    mmap_mode: {None, 'r+', 'r', 'w+', 'c'}, optional
        If not None, then memory-map the file, using the given mode
        (see `numpy.memmap`).  The mode has no effect for pickled or
        zipped files.
        A memory-mapped array is stored on disk, and not directly loaded
        into memory.  However, it can be accessed and sliced like any
        ndarray.  Memory mapping is especially useful for accessing
        small fragments of large files without reading the entire file
        into memory.

    Returns
    -------
    result : array, tuple, dict, etc.
        Data stored in the file.

    Raises
    ------
    IOError
        If the input file does not exist or cannot be read.

    See Also
    --------
    save, savez, loadtxt
    memmap : Create a memory-map to an array stored in a file on disk.

    Notes
    -----
    - If the file contains pickle data, then whatever is stored in the
      pickle is returned.
    - If the file is a ``.npy`` file, then an array is returned.
    - If the file is a ``.npz`` file, then a dictionary-like object is
      returned, containing ``{filename: array}`` key-value pairs, one for
      each file in the archive.

    Examples
    --------
    Store data to disk, and load it again:

    >>> np.save('/tmp/123', np.array([[1, 2, 3], [4, 5, 6]]))
    >>> np.load('/tmp/123.npy')
    array([[1, 2, 3],
           [4, 5, 6]])

    Mem-map the stored array, and then access the second row
    directly from disk:

    >>> X = np.load('/tmp/123.npy', mmap_mode='r')
    >>> X[1, :]
    memmap([4, 5, 6])
    
    """
    raise NotImplementedError


def loads():
    """loads(string) -- Load a pickle from the given string
    
    """
    # BUILTIN
    raise NotImplementedError


def loadtxt(fname, dtype=float, comments='#', delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0):
    """Load data from a text file.

    Each row in the text file must have the same number of values.

    Parameters
    ----------
    fname : file or str
        File, filename, or generator to read.  If the filename extension is
        ``.gz`` or ``.bz2``, the file is first decompressed. Note that
        generators should return byte strings for Python 3k.
    dtype : data-type, optional
        Data-type of the resulting array; default: float.  If this is a
        record data-type, the resulting array will be 1-dimensional, and
        each row will be interpreted as an element of the array.  In this
        case, the number of columns used must match the number of fields in
        the data-type.
    comments : str, optional
        The character used to indicate the start of a comment;
        default: '#'.
    delimiter : str, optional
        The string used to separate values.  By default, this is any
        whitespace.
    converters : dict, optional
        A dictionary mapping column number to a function that will convert
        that column to a float.  E.g., if column 0 is a date string:
        ``converters = {0: datestr2num}``.  Converters can also be used to
        provide a default value for missing data (but see also `genfromtxt`):
        ``converters = {3: lambda s: float(s.strip() or 0)}``.  Default: None.
    skiprows : int, optional
        Skip the first `skiprows` lines; default: 0.
    usecols : sequence, optional
        Which columns to read, with 0 being the first.  For example,
        ``usecols = (1,4,5)`` will extract the 2nd, 5th and 6th columns.
        The default, None, results in all columns being read.
    unpack : bool, optional
        If True, the returned array is transposed, so that arguments may be
        unpacked using ``x, y, z = loadtxt(...)``.  When used with a record
        data-type, arrays are returned for each field.  Default is False.
    ndmin : int, optional
        The returned array will have at least `ndmin` dimensions.
        Otherwise mono-dimensional axes will be squeezed. 
        Legal values: 0 (default), 1 or 2.
        .. versionadded:: 1.6.0

    Returns
    -------
    out : ndarray
        Data read from the text file.

    See Also
    --------
    load, fromstring, fromregex
    genfromtxt : Load data with missing values handled as specified.
    scipy.io.loadmat : reads MATLAB data files

    Notes
    -----
    This function aims to be a fast reader for simply formatted files.  The
    `genfromtxt` function provides more sophisticated handling of, e.g.,
    lines with missing values.

    Examples
    --------
    >>> from StringIO import StringIO   # StringIO behaves like a file object
    >>> c = StringIO("0 1\n2 3")
    >>> np.loadtxt(c)
    array([[ 0.,  1.],
           [ 2.,  3.]])

    >>> d = StringIO("M 21 72\nF 35 58")
    >>> np.loadtxt(d, dtype={'names': ('gender', 'age', 'weight'),
    ...                      'formats': ('S1', 'i4', 'f4')})
    array([('M', 21, 72.0), ('F', 35, 58.0)],
          dtype=[('gender', '|S1'), ('age', '<i4'), ('weight', '<f4')])

    >>> c = StringIO("1,0,2\n3,0,4")
    >>> x, y = np.loadtxt(c, delimiter=',', usecols=(0, 2), unpack=True)
    >>> x
    array([ 1.,  3.])
    >>> y
    array([ 2.,  4.])
    
    """
    raise NotImplementedError


def lookfor(what, module=None, import_modules=True, regenerate=False, output=None):
    """Do a keyword search on docstrings.

    A list of of objects that matched the search is displayed,
    sorted by relevance. All given keywords need to be found in the
    docstring for it to be returned as a result, but the order does
    not matter.

    Parameters
    ----------
    what : str
        String containing words to look for.
    module : str or list, optional
        Name of module(s) whose docstrings to go through.
    import_modules : bool, optional
        Whether to import sub-modules in packages. Default is True.
    regenerate : bool, optional
        Whether to re-generate the docstring cache. Default is False.
    output : file-like, optional
        File-like object to write the output to. If omitted, use a pager.

    See Also
    --------
    source, info

    Notes
    -----
    Relevance is determined only roughly, by checking if the keywords occur
    in the function name, at the start of a docstring, etc.

    Examples
    --------
    >>> np.lookfor('binary representation')
    Search results for 'binary representation'
    ------------------------------------------
    numpy.binary_repr
        Return the binary representation of the input number as a string.
    numpy.core.setup_common.long_double_representation
        Given a binary dump as given by GNU od -b, look for long double
    numpy.base_repr
        Return a string representation of a number in the given base system.
    ...
    
    """
    raise NotImplementedError


def mafromtxt(fname):
    """Load ASCII data stored in a text file and return a masked array.

    For a complete description of all the input parameters, see `genfromtxt`.

    See Also
    --------
    numpy.genfromtxt : generic function to load ASCII data.
    
    """
    raise NotImplementedError


def mask_indices(n, mask_func, k=0):
    """Return the indices to access (n, n) arrays, given a masking function.

    Assume `mask_func` is a function that, for a square array a of size
    ``(n, n)`` with a possible offset argument `k`, when called as
    ``mask_func(a, k)`` returns a new array with zeros in certain locations
    (functions like `triu` or `tril` do precisely this). Then this function
    returns the indices where the non-zero values would be located.

    Parameters
    ----------
    n : int
        The returned indices will be valid to access arrays of shape (n, n).
    mask_func : callable
        A function whose call signature is similar to that of `triu`, `tril`.
        That is, ``mask_func(x, k)`` returns a boolean array, shaped like `x`.
        `k` is an optional argument to the function.
    k : scalar
        An optional argument which is passed through to `mask_func`. Functions
        like `triu`, `tril` take a second argument that is interpreted as an
        offset.

    Returns
    -------
    indices : tuple of arrays.
        The `n` arrays of indices corresponding to the locations where
        ``mask_func(np.ones((n, n)), k)`` is True.

    See Also
    --------
    triu, tril, triu_indices, tril_indices

    Notes
    -----
    .. versionadded:: 1.4.0

    Examples
    --------
    These are the indices that would allow you to access the upper triangular
    part of any 3x3 array:

    >>> iu = np.mask_indices(3, np.triu)

    For example, if `a` is a 3x3 array:

    >>> a = np.arange(9).reshape(3, 3)
    >>> a
    array([[0, 1, 2],
           [3, 4, 5],
           [6, 7, 8]])
    >>> a[iu]
    array([0, 1, 2, 4, 5, 8])

    An offset can be passed also to the masking function.  This gets us the
    indices starting on the first diagonal right of the main one:

    >>> iu1 = np.mask_indices(3, np.triu, 1)

    with which we now extract only three elements:

    >>> a[iu1]
    array([1, 2, 5])
    
    """
    raise NotImplementedError


def mat(data, dtype=None):
    """Interpret the input as a matrix.

    Unlike `matrix`, `asmatrix` does not make a copy if the input is already
    a matrix or an ndarray.  Equivalent to ``matrix(data, copy=False)``.

    Parameters
    ----------
    data : array_like
        Input data.

    Returns
    -------
    mat : matrix
        `data` interpreted as a matrix.

    Examples
    --------
    >>> x = np.array([[1, 2], [3, 4]])

    >>> m = np.asmatrix(x)

    >>> x[0,0] = 5

    >>> m
    matrix([[5, 2],
            [3, 4]])
    
    """
    raise NotImplementedError


def max(a, axis=None, out=None):
    """Return the maximum of an array or maximum along an axis.

    Parameters
    ----------
    a : array_like
        Input data.
    axis : int, optional
        Axis along which to operate.  By default flattened input is used.
    out : ndarray, optional
        Alternate output array in which to place the result.  Must be of
        the same shape and buffer length as the expected output.  See
        `doc.ufuncs` (Section "Output arguments") for more details.

    Returns
    -------
    amax : ndarray or scalar
        Maximum of `a`. If `axis` is None, the result is a scalar value.
        If `axis` is given, the result is an array of dimension
        ``a.ndim - 1``.

    See Also
    --------
    nanmax : NaN values are ignored instead of being propagated.
    fmax : same behavior as the C99 fmax function.
    argmax : indices of the maximum values.

    Notes
    -----
    NaN values are propagated, that is if at least one item is NaN, the
    corresponding max value will be NaN as well.  To ignore NaN values
    (MATLAB behavior), please use nanmax.

    Examples
    --------
    >>> a = np.arange(4).reshape((2,2))
    >>> a
    array([[0, 1],
           [2, 3]])
    >>> np.amax(a)
    3
    >>> np.amax(a, axis=0)
    array([2, 3])
    >>> np.amax(a, axis=1)
    array([1, 3])

    >>> b = np.arange(5, dtype=np.float)
    >>> b[2] = np.NaN
    >>> np.amax(b)
    nan
    >>> np.nanmax(b)
    4.0
    
    """
    raise NotImplementedError


def maximum_sctype(t):
    """Return the scalar type of highest precision of the same kind as the input.

    Parameters
    ----------
    t : dtype or dtype specifier
        The input data type. This can be a `dtype` object or an object that
        is convertible to a `dtype`.

    Returns
    -------
    out : dtype
        The highest precision data type of the same kind (`dtype.kind`) as `t`.

    See Also
    --------
    obj2sctype, mintypecode, sctype2char
    dtype

    Examples
    --------
    >>> np.maximum_sctype(np.int)
    <type 'numpy.int64'>
    >>> np.maximum_sctype(np.uint8)
    <type 'numpy.uint64'>
    >>> np.maximum_sctype(np.complex)
    <type 'numpy.complex192'>

    >>> np.maximum_sctype(str)
    <type 'numpy.string_'>

    >>> np.maximum_sctype('i2')
    <type 'numpy.int64'>
    >>> np.maximum_sctype('f4')
    <type 'numpy.float96'>
    
    """
    raise NotImplementedError


def may_share_memory(a, b):
    """Determine if two arrays can share memory

    The memory-bounds of a and b are computed.  If they overlap then
    this function returns True.  Otherwise, it returns False.

    A return of True does not necessarily mean that the two arrays
    share any element.  It just means that they *might*.

    Parameters
    ----------
    a, b : ndarray

    Returns
    -------
    out : bool

    Examples
    --------
    >>> np.may_share_memory(np.array([1,2]), np.array([5,8,9]))
    False
    
    """
    raise NotImplementedError


def mean(a, axis=None, dtype=None, out=None):
    """Compute the arithmetic mean along the specified axis.

    Returns the average of the array elements.  The average is taken over
    the flattened array by default, otherwise over the specified axis.
    `float64` intermediate and return values are used for integer inputs.

    Parameters
    ----------
    a : array_like
        Array containing numbers whose mean is desired. If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the means are computed. The default is to compute
        the mean of the flattened array.
    dtype : data-type, optional
        Type to use in computing the mean.  For integer inputs, the default
        is `float64`; for floating point inputs, it is the same as the
        input dtype.
    out : ndarray, optional
        Alternate output array in which to place the result.  The default
        is ``None``; if provided, it must have the same shape as the
        expected output, but the type will be cast if necessary.
        See `doc.ufuncs` for details.

    Returns
    -------
    m : ndarray, see dtype parameter above
        If `out=None`, returns a new array containing the mean values,
        otherwise a reference to the output array is returned.

    See Also
    --------
    average : Weighted average

    Notes
    -----
    The arithmetic mean is the sum of the elements along the axis divided
    by the number of elements.

    Note that for floating-point input, the mean is computed using the
    same precision the input has.  Depending on the input data, this can
    cause the results to be inaccurate, especially for `float32` (see
    example below).  Specifying a higher-precision accumulator using the
    `dtype` keyword can alleviate this issue.

    Examples
    --------
    >>> a = np.array([[1, 2], [3, 4]])
    >>> np.mean(a)
    2.5
    >>> np.mean(a, axis=0)
    array([ 2.,  3.])
    >>> np.mean(a, axis=1)
    array([ 1.5,  3.5])

    In single precision, `mean` can be inaccurate:

    >>> a = np.zeros((2, 512*512), dtype=np.float32)
    >>> a[0, :] = 1.0
    >>> a[1, :] = 0.1
    >>> np.mean(a)
    0.546875

    Computing the mean in float64 is more accurate:

    >>> np.mean(a, dtype=np.float64)
    0.55000000074505806
    
    """
    raise NotImplementedError


def median(a, axis=None, out=None, overwrite_input=False):
    """Compute the median along the specified axis.

    Returns the median of the array elements.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : {None, int}, optional
        Axis along which the medians are computed. The default (axis=None)
        is to compute the median along a flattened version of the array.
    out : ndarray, optional
        Alternative output array in which to place the result. It must
        have the same shape and buffer length as the expected output,
        but the type (of the output) will be cast if necessary.
    overwrite_input : {False, True}, optional
       If True, then allow use of memory of input array (a) for
       calculations. The input array will be modified by the call to
       median. This will save memory when you do not need to preserve
       the contents of the input array. Treat the input as undefined,
       but it will probably be fully or partially sorted. Default is
       False. Note that, if `overwrite_input` is True and the input
       is not already an ndarray, an error will be raised.

    Returns
    -------
    median : ndarray
        A new array holding the result (unless `out` is specified, in
        which case that array is returned instead).  If the input contains
        integers, or floats of smaller precision than 64, then the output
        data-type is float64.  Otherwise, the output data-type is the same
        as that of the input.

    See Also
    --------
    mean, percentile

    Notes
    -----
    Given a vector V of length N, the median of V is the middle value of
    a sorted copy of V, ``V_sorted`` - i.e., ``V_sorted[(N-1)/2]``, when N is
    odd.  When N is even, it is the average of the two middle values of
    ``V_sorted``.

    Examples
    --------
    >>> a = np.array([[10, 7, 4], [3, 2, 1]])
    >>> a
    array([[10,  7,  4],
           [ 3,  2,  1]])
    >>> np.median(a)
    3.5
    >>> np.median(a, axis=0)
    array([ 6.5,  4.5,  2.5])
    >>> np.median(a, axis=1)
    array([ 7.,  2.])
    >>> m = np.median(a, axis=0)
    >>> out = np.zeros_like(m)
    >>> np.median(a, axis=0, out=m)
    array([ 6.5,  4.5,  2.5])
    >>> m
    array([ 6.5,  4.5,  2.5])
    >>> b = a.copy()
    >>> np.median(b, axis=1, overwrite_input=True)
    array([ 7.,  2.])
    >>> assert not np.all(a==b)
    >>> b = a.copy()
    >>> np.median(b, axis=None, overwrite_input=True)
    3.5
    >>> assert not np.all(a==b)
    
    """
    raise NotImplementedError


def meshgrid(x, y):
    """Return coordinate matrices from two coordinate vectors.

    Parameters
    ----------
    x, y : ndarray
        Two 1-D arrays representing the x and y coordinates of a grid.

    Returns
    -------
    X, Y : ndarray
        For vectors `x`, `y` with lengths ``Nx=len(x)`` and ``Ny=len(y)``,
        return `X`, `Y` where `X` and `Y` are ``(Ny, Nx)`` shaped arrays
        with the elements of `x` and y repeated to fill the matrix along
        the first dimension for `x`, the second for `y`.

    See Also
    --------
    index_tricks.mgrid : Construct a multi-dimensional "meshgrid"
                         using indexing notation.
    index_tricks.ogrid : Construct an open multi-dimensional "meshgrid"
                         using indexing notation.

    Examples
    --------
    >>> X, Y = np.meshgrid([1,2,3], [4,5,6,7])
    >>> X
    array([[1, 2, 3],
           [1, 2, 3],
           [1, 2, 3],
           [1, 2, 3]])
    >>> Y
    array([[4, 4, 4],
           [5, 5, 5],
           [6, 6, 6],
           [7, 7, 7]])

    `meshgrid` is very useful to evaluate functions on a grid.

    >>> x = np.arange(-5, 5, 0.1)
    >>> y = np.arange(-5, 5, 0.1)
    >>> xx, yy = np.meshgrid(x, y)
    >>> z = np.sin(xx**2+yy**2)/(xx**2+yy**2)
    
    """
    raise NotImplementedError


def min(a, axis=None, out=None):
    """Return the minimum of an array or minimum along an axis.

    Parameters
    ----------
    a : array_like
        Input data.
    axis : int, optional
        Axis along which to operate.  By default a flattened input is used.
    out : ndarray, optional
        Alternative output array in which to place the result.  Must
        be of the same shape and buffer length as the expected output.
        See `doc.ufuncs` (Section "Output arguments") for more details.

    Returns
    -------
    amin : ndarray
        A new array or a scalar array with the result.

    See Also
    --------
    nanmin: nan values are ignored instead of being propagated
    fmin: same behavior as the C99 fmin function
    argmin: Return the indices of the minimum values.

    amax, nanmax, fmax

    Notes
    -----
    NaN values are propagated, that is if at least one item is nan, the
    corresponding min value will be nan as well. To ignore NaN values (matlab
    behavior), please use nanmin.

    Examples
    --------
    >>> a = np.arange(4).reshape((2,2))
    >>> a
    array([[0, 1],
           [2, 3]])
    >>> np.amin(a)           # Minimum of the flattened array
    0
    >>> np.amin(a, axis=0)         # Minima along the first axis
    array([0, 1])
    >>> np.amin(a, axis=1)         # Minima along the second axis
    array([0, 2])

    >>> b = np.arange(5, dtype=np.float)
    >>> b[2] = np.NaN
    >>> np.amin(b)
    nan
    >>> np.nanmin(b)
    0.0
    
    """
    raise NotImplementedError


def min_scalar_type(a):
    """min_scalar_type(a)

    For scalar ``a``, returns the data type with the smallest size
    and smallest scalar kind which can hold its value.  For non-scalar
    array ``a``, returns the vector's dtype unmodified.

    Floating point values are not demoted to integers,
    and complex values are not demoted to floats.

    Parameters
    ----------
    a : scalar or array_like
        The value whose minimal data type is to be found.

    Returns
    -------
    out : dtype
        The minimal data type.

    Notes
    -----
    .. versionadded:: 1.6.0

    See Also
    --------
    result_type, promote_types, dtype, can_cast

    Examples
    --------
    >>> np.min_scalar_type(10)
    dtype('uint8')

    >>> np.min_scalar_type(-260)
    dtype('int16')

    >>> np.min_scalar_type(3.1)
    dtype('float16')

    >>> np.min_scalar_type(1e50)
    dtype('float64')

    >>> np.min_scalar_type(np.arange(4,dtype='f8'))
    dtype('float64')
    
    """
    # BUILTIN
    raise NotImplementedError


def mintypecode(typechars, typeset='GDFgdf', default='d'):
    """Return the character for the minimum-size type to which given types can
    be safely cast.

    The returned type character must represent the smallest size dtype such
    that an array of the returned type can handle the data from an array of
    all types in `typechars` (or if `typechars` is an array, then its
    dtype.char).

    Parameters
    ----------
    typechars : list of str or array_like
        If a list of strings, each string should represent a dtype.
        If array_like, the character representation of the array dtype is used.
    typeset : str or list of str, optional
        The set of characters that the returned character is chosen from.
        The default set is 'GDFgdf'.
    default : str, optional
        The default character, this is returned if none of the characters in
        `typechars` matches a character in `typeset`.

    Returns
    -------
    typechar : str
        The character representing the minimum-size type that was found.

    See Also
    --------
    dtype, sctype2char, maximum_sctype

    Examples
    --------
    >>> np.mintypecode(['d', 'f', 'S'])
    'd'
    >>> x = np.array([1.1, 2-3.j])
    >>> np.mintypecode(x)
    'D'

    >>> np.mintypecode('abceh', default='G')
    'G'
    
    """
    raise NotImplementedError


def mirr(values, finance_rate, reinvest_rate):
    """Modified internal rate of return.

    Parameters
    ----------
    values : array_like
        Cash flows (must contain at least one positive and one negative value)
        or nan is returned.  The first value is considered a sunk cost at time zero.
    finance_rate : scalar
        Interest rate paid on the cash flows
    reinvest_rate : scalar
        Interest rate received on the cash flows upon reinvestment

    Returns
    -------
    out : float
        Modified internal rate of return
    
    """
    raise NotImplementedError


def msort(a):
    """Return a copy of an array sorted along the first axis.

    Parameters
    ----------
    a : array_like
        Array to be sorted.

    Returns
    -------
    sorted_array : ndarray
        Array of the same type and shape as `a`.

    See Also
    --------
    sort

    Notes
    -----
    ``np.msort(a)`` is equivalent to  ``np.sort(a, axis=0)``.
    
    """
    raise NotImplementedError


def nan_to_num(x):
    """Replace nan with zero and inf with finite numbers.

    Returns an array or scalar replacing Not a Number (NaN) with zero,
    (positive) infinity with a very large number and negative infinity
    with a very small (or negative) number.

    Parameters
    ----------
    x : array_like
        Input data.

    Returns
    -------
    out : ndarray, float
        Array with the same shape as `x` and dtype of the element in `x`  with
        the greatest precision. NaN is replaced by zero, and infinity
        (-infinity) is replaced by the largest (smallest or most negative)
        floating point value that fits in the output dtype. All finite numbers
        are upcast to the output dtype (default float64).

    See Also
    --------
    isinf : Shows which elements are negative or negative infinity.
    isneginf : Shows which elements are negative infinity.
    isposinf : Shows which elements are positive infinity.
    isnan : Shows which elements are Not a Number (NaN).
    isfinite : Shows which elements are finite (not NaN, not infinity)

    Notes
    -----
    Numpy uses the IEEE Standard for Binary Floating-Point for Arithmetic
    (IEEE 754). This means that Not a Number is not equivalent to infinity.


    Examples
    --------
    >>> np.set_printoptions(precision=8)
    >>> x = np.array([np.inf, -np.inf, np.nan, -128, 128])
    >>> np.nan_to_num(x)
    array([  1.79769313e+308,  -1.79769313e+308,   0.00000000e+000,
            -1.28000000e+002,   1.28000000e+002])
    
    """
    raise NotImplementedError


def nanargmax(a, axis=None):
    """Return indices of the maximum values over an axis, ignoring NaNs.

    Parameters
    ----------
    a : array_like
        Input data.
    axis : int, optional
        Axis along which to operate.  By default flattened input is used.

    Returns
    -------
    index_array : ndarray
        An array of indices or a single index value.

    See Also
    --------
    argmax, nanargmin

    Examples
    --------
    >>> a = np.array([[np.nan, 4], [2, 3]])
    >>> np.argmax(a)
    0
    >>> np.nanargmax(a)
    1
    >>> np.nanargmax(a, axis=0)
    array([1, 0])
    >>> np.nanargmax(a, axis=1)
    array([1, 1])
    
    """
    raise NotImplementedError


def nanargmin(a, axis=None):
    """Return indices of the minimum values over an axis, ignoring NaNs.

    Parameters
    ----------
    a : array_like
        Input data.
    axis : int, optional
        Axis along which to operate.  By default flattened input is used.

    Returns
    -------
    index_array : ndarray
        An array of indices or a single index value.

    See Also
    --------
    argmin, nanargmax

    Examples
    --------
    >>> a = np.array([[np.nan, 4], [2, 3]])
    >>> np.argmin(a)
    0
    >>> np.nanargmin(a)
    2
    >>> np.nanargmin(a, axis=0)
    array([1, 1])
    >>> np.nanargmin(a, axis=1)
    array([1, 0])
    
    """
    raise NotImplementedError


def nanmax(a, axis=None):
    """Return the maximum of an array or maximum along an axis ignoring any NaNs.

    Parameters
    ----------
    a : array_like
        Array containing numbers whose maximum is desired. If `a` is not
        an array, a conversion is attempted.
    axis : int, optional
        Axis along which the maximum is computed. The default is to compute
        the maximum of the flattened array.

    Returns
    -------
    nanmax : ndarray
        An array with the same shape as `a`, with the specified axis removed.
        If `a` is a 0-d array, or if axis is None, a ndarray scalar is
        returned.  The the same dtype as `a` is returned.

    See Also
    --------
    numpy.amax : Maximum across array including any Not a Numbers.
    numpy.nanmin : Minimum across array ignoring any Not a Numbers.
    isnan : Shows which elements are Not a Number (NaN).
    isfinite: Shows which elements are not: Not a Number, positive and
             negative infinity

    Notes
    -----
    Numpy uses the IEEE Standard for Binary Floating-Point for Arithmetic
    (IEEE 754). This means that Not a Number is not equivalent to infinity.
    Positive infinity is treated as a very large number and negative infinity
    is treated as a very small (i.e. negative) number.

    If the input has a integer type the function is equivalent to np.max.

    Examples
    --------
    >>> a = np.array([[1, 2], [3, np.nan]])
    >>> np.nanmax(a)
    3.0
    >>> np.nanmax(a, axis=0)
    array([ 3.,  2.])
    >>> np.nanmax(a, axis=1)
    array([ 2.,  3.])

    When positive infinity and negative infinity are present:

    >>> np.nanmax([1, 2, np.nan, np.NINF])
    2.0
    >>> np.nanmax([1, 2, np.nan, np.inf])
    inf
    
    """
    raise NotImplementedError


def nanmin(a, axis=None):
    """Return the minimum of an array or minimum along an axis ignoring any NaNs.

    Parameters
    ----------
    a : array_like
        Array containing numbers whose minimum is desired.
    axis : int, optional
        Axis along which the minimum is computed.The default is to compute
        the minimum of the flattened array.

    Returns
    -------
    nanmin : ndarray
        A new array or a scalar array with the result.

    See Also
    --------
    numpy.amin : Minimum across array including any Not a Numbers.
    numpy.nanmax : Maximum across array ignoring any Not a Numbers.
    isnan : Shows which elements are Not a Number (NaN).
    isfinite: Shows which elements are not: Not a Number, positive and
             negative infinity

    Notes
    -----
    Numpy uses the IEEE Standard for Binary Floating-Point for Arithmetic
    (IEEE 754). This means that Not a Number is not equivalent to infinity.
    Positive infinity is treated as a very large number and negative infinity
    is treated as a very small (i.e. negative) number.

    If the input has a integer type the function is equivalent to np.min.


    Examples
    --------
    >>> a = np.array([[1, 2], [3, np.nan]])
    >>> np.nanmin(a)
    1.0
    >>> np.nanmin(a, axis=0)
    array([ 1.,  2.])
    >>> np.nanmin(a, axis=1)
    array([ 1.,  3.])

    When positive infinity and negative infinity are present:

    >>> np.nanmin([1, 2, np.nan, np.inf])
    1.0
    >>> np.nanmin([1, 2, np.nan, np.NINF])
    -inf
    
    """
    raise NotImplementedError


def nansum(a, axis=None):
    """Return the sum of array elements over a given axis treating
    Not a Numbers (NaNs) as zero.

    Parameters
    ----------
    a : array_like
        Array containing numbers whose sum is desired. If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the sum is computed. The default is to compute
        the sum of the flattened array.

    Returns
    -------
    y : ndarray
        An array with the same shape as a, with the specified axis removed.
        If a is a 0-d array, or if axis is None, a scalar is returned with
        the same dtype as `a`.

    See Also
    --------
    numpy.sum : Sum across array including Not a Numbers.
    isnan : Shows which elements are Not a Number (NaN).
    isfinite: Shows which elements are not: Not a Number, positive and
             negative infinity

    Notes
    -----
    Numpy uses the IEEE Standard for Binary Floating-Point for Arithmetic
    (IEEE 754). This means that Not a Number is not equivalent to infinity.
    If positive or negative infinity are present the result is positive or
    negative infinity. But if both positive and negative infinity are present,
    the result is Not A Number (NaN).

    Arithmetic is modular when using integer types (all elements of `a` must
    be finite i.e. no elements that are NaNs, positive infinity and negative
    infinity because NaNs are floating point types), and no error is raised
    on overflow.


    Examples
    --------
    >>> np.nansum(1)
    1
    >>> np.nansum([1])
    1
    >>> np.nansum([1, np.nan])
    1.0
    >>> a = np.array([[1, 1], [1, np.nan]])
    >>> np.nansum(a)
    3.0
    >>> np.nansum(a, axis=0)
    array([ 2.,  1.])

    When positive infinity and negative infinity are present

    >>> np.nansum([1, np.nan, np.inf])
    inf
    >>> np.nansum([1, np.nan, np.NINF])
    -inf
    >>> np.nansum([1, np.nan, np.inf, np.NINF])
    nan
    
    """
    raise NotImplementedError


def ndfromtxt(fname):
    """Load ASCII data stored in a file and return it as a single array.

    Complete description of all the optional input parameters is available in
    the docstring of the `genfromtxt` function.

    See Also
    --------
    numpy.genfromtxt : generic function.
    
    """
    raise NotImplementedError


def ndim(a):
    """Return the number of dimensions of an array.

    Parameters
    ----------
    a : array_like
        Input array.  If it is not already an ndarray, a conversion is
        attempted.

    Returns
    -------
    number_of_dimensions : int
        The number of dimensions in `a`.  Scalars are zero-dimensional.

    See Also
    --------
    ndarray.ndim : equivalent method
    shape : dimensions of array
    ndarray.shape : dimensions of array

    Examples
    --------
    >>> np.ndim([[1,2,3],[4,5,6]])
    2
    >>> np.ndim(np.array([[1,2,3],[4,5,6]]))
    2
    >>> np.ndim(1)
    0
    
    """
    raise NotImplementedError


def nested_iters():
    """
    
    """
    # BUILTIN
    raise NotImplementedError


def newbuffer(size):
    """newbuffer(size)

    Return a new uninitialized buffer object of size bytes
    
    """
    # BUILTIN
    raise NotImplementedError


def nonzero(a):
    """Return the indices of the elements that are non-zero.

    Returns a tuple of arrays, one for each dimension of `a`, containing
    the indices of the non-zero elements in that dimension. The
    corresponding non-zero values can be obtained with::

        a[nonzero(a)]

    To group the indices by element, rather than dimension, use::

        transpose(nonzero(a))

    The result of this is always a 2-D array, with a row for
    each non-zero element.

    Parameters
    ----------
    a : array_like
        Input array.

    Returns
    -------
    tuple_of_arrays : tuple
        Indices of elements that are non-zero.

    See Also
    --------
    flatnonzero :
        Return indices that are non-zero in the flattened version of the input
        array.
    ndarray.nonzero :
        Equivalent ndarray method.
    count_nonzero :
        Counts the number of non-zero elements in the input array.

    Examples
    --------
    >>> x = np.eye(3)
    >>> x
    array([[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]])
    >>> np.nonzero(x)
    (array([0, 1, 2]), array([0, 1, 2]))

    >>> x[np.nonzero(x)]
    array([ 1.,  1.,  1.])
    >>> np.transpose(np.nonzero(x))
    array([[0, 0],
           [1, 1],
           [2, 2]])

    A common use for ``nonzero`` is to find the indices of an array, where
    a condition is True.  Given an array `a`, the condition `a` > 3 is a
    boolean array and since False is interpreted as 0, np.nonzero(a > 3)
    yields the indices of the `a` where the condition is true.

    >>> a = np.array([[1,2,3],[4,5,6],[7,8,9]])
    >>> a > 3
    array([[False, False, False],
           [ True,  True,  True],
           [ True,  True,  True]], dtype=bool)
    >>> np.nonzero(a > 3)
    (array([1, 1, 1, 2, 2, 2]), array([0, 1, 2, 0, 1, 2]))

    The ``nonzero`` method of the boolean array can also be called.

    >>> (a > 3).nonzero()
    (array([1, 1, 1, 2, 2, 2]), array([0, 1, 2, 0, 1, 2]))
    
    """
    raise NotImplementedError


def nper(rate, pmt, pv, fv=0, when='end'):
    """Compute the number of periodic payments.

    Parameters
    ----------
    rate : array_like
        Rate of interest (per period)
    pmt : array_like
        Payment
    pv : array_like
        Present value
    fv : array_like, optional
        Future value
    when : {{'begin', 1}, {'end', 0}}, {string, int}, optional
        When payments are due ('begin' (1) or 'end' (0))

    Notes
    -----
    The number of periods ``nper`` is computed by solving the equation::

     fv + pv*(1+rate)**nper + pmt*(1+rate*when)/rate*((1+rate)**nper-1) = 0

    but if ``rate = 0`` then::

     fv + pv + pmt*nper = 0

    Examples
    --------
    If you only had $150/month to pay towards the loan, how long would it take
    to pay-off a loan of $8,000 at 7% annual interest?

    >>> np.nper(0.07/12, -150, 8000)
    64.073348770661852

    So, over 64 months would be required to pay off the loan.

    The same analysis could be done with several different interest rates
    and/or payments and/or total amounts to produce an entire table.

    >>> np.nper(*(np.ogrid[0.07/12: 0.08/12: 0.01/12,
    ...                    -150   : -99     : 50    ,
    ...                    8000   : 9001    : 1000]))
    array([[[  64.07334877,   74.06368256],
            [ 108.07548412,  127.99022654]],
           [[  66.12443902,   76.87897353],
            [ 114.70165583,  137.90124779]]])
    
    """
    raise NotImplementedError


def npv(rate, values):
    """Returns the NPV (Net Present Value) of a cash flow series.

    Parameters
    ----------
    rate : scalar
        The discount rate.
    values : array_like, shape(M, )
        The values of the time series of cash flows.  The (fixed) time
        interval between cash flow "events" must be the same as that
        for which `rate` is given (i.e., if `rate` is per year, then
        precisely a year is understood to elapse between each cash flow
        event).  By convention, investments or "deposits" are negative,
        income or "withdrawals" are positive; `values` must begin with
        the initial investment, thus `values[0]` will typically be
        negative.

    Returns
    -------
    out : float
        The NPV of the input cash flow series `values` at the discount `rate`.

    Notes
    -----
    Returns the result of: [G]_

    .. math :: \sum_{t=0}^M{\frac{values_t}{(1+rate)^{t}}}

    References
    ----------
    .. [G] L. J. Gitman, "Principles of Managerial Finance, Brief," 3rd ed.,
       Addison-Wesley, 2003, pg. 346.

    Examples
    --------
    >>> np.npv(0.281,[-100, 39, 59, 55, 20])
    -0.0066187288356340801

    (Compare with the Example given for numpy.lib.financial.irr)
    
    """
    raise NotImplementedError


def obj2sctype(rep, default=None):
    """Return the scalar dtype or NumPy equivalent of Python type of an object.

    Parameters
    ----------
    rep : any
        The object of which the type is returned.
    default : any, optional
        If given, this is returned for objects whose types can not be
        determined. If not given, None is returned for those objects.

    Returns
    -------
    dtype : dtype or Python type
        The data type of `rep`.

    See Also
    --------
    sctype2char, issctype, issubsctype, issubdtype, maximum_sctype

    Examples
    --------
    >>> np.obj2sctype(np.int32)
    <type 'numpy.int32'>
    >>> np.obj2sctype(np.array([1., 2.]))
    <type 'numpy.float64'>
    >>> np.obj2sctype(np.array([1.j]))
    <type 'numpy.complex128'>

    >>> np.obj2sctype(dict)
    <type 'numpy.object_'>
    >>> np.obj2sctype('string')
    <type 'numpy.string_'>

    >>> np.obj2sctype(1, default=list)
    <type 'list'>
    
    """
    raise NotImplementedError


def outer(a, b):
    """Compute the outer product of two vectors.

    Given two vectors, ``a = [a0, a1, ..., aM]`` and
    ``b = [b0, b1, ..., bN]``,
    the outer product [1]_ is::

      [[a0*b0  a0*b1 ... a0*bN ]
       [a1*b0    .
       [ ...          .
       [aM*b0            aM*bN ]]

    Parameters
    ----------
    a, b : array_like, shape (M,), (N,)
        First and second input vectors.  Inputs are flattened if they
        are not already 1-dimensional.

    Returns
    -------
    out : ndarray, shape (M, N)
        ``out[i, j] = a[i] * b[j]``

    See also
    --------
    inner, einsum

    References
    ----------
    .. [1] : G. H. Golub and C. F. van Loan, *Matrix Computations*, 3rd
             ed., Baltimore, MD, Johns Hopkins University Press, 1996,
             pg. 8.

    Examples
    --------
    Make a (*very* coarse) grid for computing a Mandelbrot set:

    >>> rl = np.outer(np.ones((5,)), np.linspace(-2, 2, 5))
    >>> rl
    array([[-2., -1.,  0.,  1.,  2.],
           [-2., -1.,  0.,  1.,  2.],
           [-2., -1.,  0.,  1.,  2.],
           [-2., -1.,  0.,  1.,  2.],
           [-2., -1.,  0.,  1.,  2.]])
    >>> im = np.outer(1j*np.linspace(2, -2, 5), np.ones((5,)))
    >>> im
    array([[ 0.+2.j,  0.+2.j,  0.+2.j,  0.+2.j,  0.+2.j],
           [ 0.+1.j,  0.+1.j,  0.+1.j,  0.+1.j,  0.+1.j],
           [ 0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j],
           [ 0.-1.j,  0.-1.j,  0.-1.j,  0.-1.j,  0.-1.j],
           [ 0.-2.j,  0.-2.j,  0.-2.j,  0.-2.j,  0.-2.j]])
    >>> grid = rl + im
    >>> grid
    array([[-2.+2.j, -1.+2.j,  0.+2.j,  1.+2.j,  2.+2.j],
           [-2.+1.j, -1.+1.j,  0.+1.j,  1.+1.j,  2.+1.j],
           [-2.+0.j, -1.+0.j,  0.+0.j,  1.+0.j,  2.+0.j],
           [-2.-1.j, -1.-1.j,  0.-1.j,  1.-1.j,  2.-1.j],
           [-2.-2.j, -1.-2.j,  0.-2.j,  1.-2.j,  2.-2.j]])

    An example using a "vector" of letters:

    >>> x = np.array(['a', 'b', 'c'], dtype=object)
    >>> np.outer(x, [1, 2, 3])
    array([[a, aa, aaa],
           [b, bb, bbb],
           [c, cc, ccc]], dtype=object)
    
    """
    raise NotImplementedError


def packbits(myarray, axis=None):
    """packbits(myarray, axis=None)

    Packs the elements of a binary-valued array into bits in a uint8 array.

    The result is padded to full bytes by inserting zero bits at the end.

    Parameters
    ----------
    myarray : array_like
        An integer type array whose elements should be packed to bits.
    axis : int, optional
        The dimension over which bit-packing is done.
        ``None`` implies packing the flattened array.

    Returns
    -------
    packed : ndarray
        Array of type uint8 whose elements represent bits corresponding to the
        logical (0 or nonzero) value of the input elements. The shape of
        `packed` has the same number of dimensions as the input (unless `axis`
        is None, in which case the output is 1-D).

    See Also
    --------
    unpackbits: Unpacks elements of a uint8 array into a binary-valued output
                array.

    Examples
    --------
    >>> a = np.array([[[1,0,1],
    ...                [0,1,0]],
    ...               [[1,1,0],
    ...                [0,0,1]]])
    >>> b = np.packbits(a, axis=-1)
    >>> b
    array([[[160],[64]],[[192],[32]]], dtype=uint8)

    Note that in binary 160 = 1010 0000, 64 = 0100 0000, 192 = 1100 0000,
    and 32 = 0010 0000.
    
    """
    # BUILTIN
    raise NotImplementedError


def percentile(a, q, axis=None, out=None, overwrite_input=False):
    """Compute the qth percentile of the data along the specified axis.

    Returns the qth percentile of the array elements.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    q : float in range of [0,100] (or sequence of floats)
        Percentile to compute which must be between 0 and 100 inclusive.
    axis : int, optional
        Axis along which the percentiles are computed. The default (None)
        is to compute the median along a flattened version of the array.
    out : ndarray, optional
        Alternative output array in which to place the result. It must
        have the same shape and buffer length as the expected output,
        but the type (of the output) will be cast if necessary.
    overwrite_input : bool, optional
       If True, then allow use of memory of input array `a` for
       calculations. The input array will be modified by the call to
       median. This will save memory when you do not need to preserve
       the contents of the input array. Treat the input as undefined,
       but it will probably be fully or partially sorted.
       Default is False. Note that, if `overwrite_input` is True and the
       input is not already an array, an error will be raised.

    Returns
    -------
    pcntile : ndarray
        A new array holding the result (unless `out` is specified, in
        which case that array is returned instead).  If the input contains
        integers, or floats of smaller precision than 64, then the output
        data-type is float64.  Otherwise, the output data-type is the same
        as that of the input.

    See Also
    --------
    mean, median

    Notes
    -----
    Given a vector V of length N, the qth percentile of V is the qth ranked
    value in a sorted copy of V.  A weighted average of the two nearest
    neighbors is used if the normalized ranking does not match q exactly.
    The same as the median if ``q=0.5``, the same as the minimum if ``q=0``
    and the same as the maximum if ``q=1``.

    Examples
    --------
    >>> a = np.array([[10, 7, 4], [3, 2, 1]])
    >>> a
    array([[10,  7,  4],
           [ 3,  2,  1]])
    >>> np.percentile(a, 50)
    3.5
    >>> np.percentile(a, 0.5, axis=0)
    array([ 6.5,  4.5,  2.5])
    >>> np.percentile(a, 50, axis=1)
    array([ 7.,  2.])

    >>> m = np.percentile(a, 50, axis=0)
    >>> out = np.zeros_like(m)
    >>> np.percentile(a, 50, axis=0, out=m)
    array([ 6.5,  4.5,  2.5])
    >>> m
    array([ 6.5,  4.5,  2.5])

    >>> b = a.copy()
    >>> np.percentile(b, 50, axis=1, overwrite_input=True)
    array([ 7.,  2.])
    >>> assert not np.all(a==b)
    >>> b = a.copy()
    >>> np.percentile(b, 50, axis=None, overwrite_input=True)
    3.5
    
    """
    raise NotImplementedError


def piecewise(x, condlist, funclist):
    """Evaluate a piecewise-defined function.

    Given a set of conditions and corresponding functions, evaluate each
    function on the input data wherever its condition is true.

    Parameters
    ----------
    x : ndarray
        The input domain.
    condlist : list of bool arrays
        Each boolean array corresponds to a function in `funclist`.  Wherever
        `condlist[i]` is True, `funclist[i](x)` is used as the output value.

        Each boolean array in `condlist` selects a piece of `x`,
        and should therefore be of the same shape as `x`.

        The length of `condlist` must correspond to that of `funclist`.
        If one extra function is given, i.e. if
        ``len(funclist) - len(condlist) == 1``, then that extra function
        is the default value, used wherever all conditions are false.
    funclist : list of callables, f(x,*args,**kw), or scalars
        Each function is evaluated over `x` wherever its corresponding
        condition is True.  It should take an array as input and give an array
        or a scalar value as output.  If, instead of a callable,
        a scalar is provided then a constant function (``lambda x: scalar``) is
        assumed.
    args : tuple, optional
        Any further arguments given to `piecewise` are passed to the functions
        upon execution, i.e., if called ``piecewise(..., ..., 1, 'a')``, then
        each function is called as ``f(x, 1, 'a')``.
    kw : dict, optional
        Keyword arguments used in calling `piecewise` are passed to the
        functions upon execution, i.e., if called
        ``piecewise(..., ..., lambda=1)``, then each function is called as
        ``f(x, lambda=1)``.

    Returns
    -------
    out : ndarray
        The output is the same shape and type as x and is found by
        calling the functions in `funclist` on the appropriate portions of `x`,
        as defined by the boolean arrays in `condlist`.  Portions not covered
        by any condition have undefined values.


    See Also
    --------
    choose, select, where

    Notes
    -----
    This is similar to choose or select, except that functions are
    evaluated on elements of `x` that satisfy the corresponding condition from
    `condlist`.

    The result is::

            |--
            |funclist[0](x[condlist[0]])
      out = |funclist[1](x[condlist[1]])
            |...
            |funclist[n2](x[condlist[n2]])
            |--

    Examples
    --------
    Define the sigma function, which is -1 for ``x < 0`` and +1 for ``x >= 0``.

    >>> x = np.arange(6) - 2.5
    >>> np.piecewise(x, [x < 0, x >= 0], [-1, 1])
    array([-1., -1., -1.,  1.,  1.,  1.])

    Define the absolute value, which is ``-x`` for ``x <0`` and ``x`` for
    ``x >= 0``.

    >>> np.piecewise(x, [x < 0, x >= 0], [lambda x: -x, lambda x: x])
    array([ 2.5,  1.5,  0.5,  0.5,  1.5,  2.5])
    
    """
    raise NotImplementedError


def pkgload():
    """Load one or more packages into parent package top-level namespace.

       This function is intended to shorten the need to import many
       subpackages, say of scipy, constantly with statements such as

         import scipy.linalg, scipy.fftpack, scipy.etc...

       Instead, you can say:

         import scipy
         scipy.pkgload('linalg','fftpack',...)

       or

         scipy.pkgload()

       to load all of them in one call.

       If a name which doesn't exist in scipy's namespace is
       given, a warning is shown.

       Parameters
       ----------
        *packages : arg-tuple
             the names (one or more strings) of all the modules one
             wishes to load into the top-level namespace.
        verbose= : integer
             verbosity level [default: -1].
             verbose=-1 will suspend also warnings.
        force= : bool
             when True, force reloading loaded packages [default: False].
        postpone= : bool
             when True, don't load packages [default: False]
    
    """
    raise NotImplementedError


def place(arr, mask, vals):
    """Change elements of an array based on conditional and input values.

    Similar to ``np.putmask(arr, mask, vals)``, the difference is that `place`
    uses the first N elements of `vals`, where N is the number of True values
    in `mask`, while `putmask` uses the elements where `mask` is True.

    Note that `extract` does the exact opposite of `place`.

    Parameters
    ----------
    arr : array_like
        Array to put data into.
    mask : array_like
        Boolean mask array. Must have the same size as `a`.
    vals : 1-D sequence
        Values to put into `a`. Only the first N elements are used, where
        N is the number of True values in `mask`. If `vals` is smaller
        than N it will be repeated.

    See Also
    --------
    putmask, put, take, extract

    Examples
    --------
    >>> arr = np.arange(6).reshape(2, 3)
    >>> np.place(arr, arr>2, [44, 55])
    >>> arr
    array([[ 0,  1,  2],
           [44, 55, 44]])
    
    """
    raise NotImplementedError


def pmt(rate, nper, pv, fv=0, when='end'):
    """Compute the payment against loan principal plus interest.

    Given:
     * a present value, `pv` (e.g., an amount borrowed)
     * a future value, `fv` (e.g., 0)
     * an interest `rate` compounded once per period, of which
       there are
     * `nper` total
     * and (optional) specification of whether payment is made
       at the beginning (`when` = {'begin', 1}) or the end
       (`when` = {'end', 0}) of each period

    Return:
       the (fixed) periodic payment.

    Parameters
    ----------
    rate : array_like
        Rate of interest (per period)
    nper : array_like
        Number of compounding periods
    pv : array_like
        Present value
    fv : array_like (optional)
        Future value (default = 0)
    when : {{'begin', 1}, {'end', 0}}, {string, int}
        When payments are due ('begin' (1) or 'end' (0))

    Returns
    -------
    out : ndarray
        Payment against loan plus interest.  If all input is scalar, returns a
        scalar float.  If any input is array_like, returns payment for each
        input element. If multiple inputs are array_like, they all must have
        the same shape.

    Notes
    -----
    The payment is computed by solving the equation::

     fv +
     pv*(1 + rate)**nper +
     pmt*(1 + rate*when)/rate*((1 + rate)**nper - 1) == 0

    or, when ``rate == 0``::

      fv + pv + pmt * nper == 0

    for ``pmt``.

    Note that computing a monthly mortgage payment is only
    one use for this function.  For example, pmt returns the
    periodic deposit one must make to achieve a specified
    future balance given an initial deposit, a fixed,
    periodically compounded interest rate, and the total
    number of periods.

    References
    ----------
    .. [WRW] Wheeler, D. A., E. Rathke, and R. Weir (Eds.) (2009, May).
       Open Document Format for Office Applications (OpenDocument)v1.2,
       Part 2: Recalculated Formula (OpenFormula) Format - Annotated Version,
       Pre-Draft 12. Organization for the Advancement of Structured Information
       Standards (OASIS). Billerica, MA, USA. [ODT Document].
       Available:
       http://www.oasis-open.org/committees/documents.php
       ?wg_abbrev=office-formulaOpenDocument-formula-20090508.odt

    Examples
    --------
    What is the monthly payment needed to pay off a $200,000 loan in 15
    years at an annual interest rate of 7.5%?

    >>> np.pmt(0.075/12, 12*15, 200000)
    -1854.0247200054619

    In order to pay-off (i.e., have a future-value of 0) the $200,000 obtained
    today, a monthly payment of $1,854.02 would be required.  Note that this
    example illustrates usage of `fv` having a default value of 0.
    
    """
    raise NotImplementedError


def poly(seq_of_zeros):
    """Find the coefficients of a polynomial with the given sequence of roots.

    Returns the coefficients of the polynomial whose leading coefficient
    is one for the given sequence of zeros (multiple roots must be included
    in the sequence as many times as their multiplicity; see Examples).
    A square matrix (or array, which will be treated as a matrix) can also
    be given, in which case the coefficients of the characteristic polynomial
    of the matrix are returned.

    Parameters
    ----------
    seq_of_zeros : array_like, shape (N,) or (N, N)
        A sequence of polynomial roots, or a square array or matrix object.

    Returns
    -------
    c : ndarray
        1D array of polynomial coefficients from highest to lowest degree:

        ``c[0] * x**(N) + c[1] * x**(N-1) + ... + c[N-1] * x + c[N]``
        where c[0] always equals 1.

    Raises
    ------
    ValueError
        If input is the wrong shape (the input must be a 1-D or square
        2-D array).

    See Also
    --------
    polyval : Evaluate a polynomial at a point.
    roots : Return the roots of a polynomial.
    polyfit : Least squares polynomial fit.
    poly1d : A one-dimensional polynomial class.

    Notes
    -----
    Specifying the roots of a polynomial still leaves one degree of
    freedom, typically represented by an undetermined leading
    coefficient. [1]_ In the case of this function, that coefficient -
    the first one in the returned array - is always taken as one. (If
    for some reason you have one other point, the only automatic way
    presently to leverage that information is to use ``polyfit``.)

    The characteristic polynomial, :math:`p_a(t)`, of an `n`-by-`n`
    matrix **A** is given by

        :math:`p_a(t) = \mathrm{det}(t\, \mathbf{I} - \mathbf{A})`,

    where **I** is the `n`-by-`n` identity matrix. [2]_

    References
    ----------
    .. [1] M. Sullivan and M. Sullivan, III, "Algebra and Trignometry,
       Enhanced With Graphing Utilities," Prentice-Hall, pg. 318, 1996.

    .. [2] G. Strang, "Linear Algebra and Its Applications, 2nd Edition,"
       Academic Press, pg. 182, 1980.

    Examples
    --------
    Given a sequence of a polynomial's zeros:

    >>> np.poly((0, 0, 0)) # Multiple root example
    array([1, 0, 0, 0])
    
    The line above represents z**3 + 0*z**2 + 0*z + 0.

    >>> np.poly((-1./2, 0, 1./2))
    array([ 1.  ,  0.  , -0.25,  0.  ])
    
    The line above represents z**3 - z/4

    >>> np.poly((np.random.random(1.)[0], 0, np.random.random(1.)[0]))
    array([ 1.        , -0.77086955,  0.08618131,  0.        ]) #random

    Given a square array object:

    >>> P = np.array([[0, 1./3], [-1./2, 0]])
    >>> np.poly(P)
    array([ 1.        ,  0.        ,  0.16666667])

    Or a square matrix object:

    >>> np.poly(np.matrix(P))
    array([ 1.        ,  0.        ,  0.16666667])

    Note how in all cases the leading coefficient is always 1.
    
    """
    raise NotImplementedError


def polyadd(a1, a2):
    """Find the sum of two polynomials.

    Returns the polynomial resulting from the sum of two input polynomials.
    Each input must be either a poly1d object or a 1D sequence of polynomial
    coefficients, from highest to lowest degree.

    Parameters
    ----------
    a1, a2 : array_like or poly1d object
        Input polynomials.

    Returns
    -------
    out : ndarray or poly1d object
        The sum of the inputs. If either input is a poly1d object, then the
        output is also a poly1d object. Otherwise, it is a 1D array of
        polynomial coefficients from highest to lowest degree.

    See Also
    --------
    poly1d : A one-dimensional polynomial class.
    poly, polyadd, polyder, polydiv, polyfit, polyint, polysub, polyval

    Examples
    --------
    >>> np.polyadd([1, 2], [9, 5, 4])
    array([9, 6, 6])

    Using poly1d objects:

    >>> p1 = np.poly1d([1, 2])
    >>> p2 = np.poly1d([9, 5, 4])
    >>> print p1
    1 x + 2
    >>> print p2
       2
    9 x + 5 x + 4
    >>> print np.polyadd(p1, p2)
       2
    9 x + 6 x + 6
    
    """
    raise NotImplementedError


def polyder(p, m=1):
    """Return the derivative of the specified order of a polynomial.

    Parameters
    ----------
    p : poly1d or sequence
        Polynomial to differentiate.
        A sequence is interpreted as polynomial coefficients, see `poly1d`.
    m : int, optional
        Order of differentiation (default: 1)

    Returns
    -------
    der : poly1d
        A new polynomial representing the derivative.

    See Also
    --------
    polyint : Anti-derivative of a polynomial.
    poly1d : Class for one-dimensional polynomials.

    Examples
    --------
    The derivative of the polynomial :math:`x^3 + x^2 + x^1 + 1` is:

    >>> p = np.poly1d([1,1,1,1])
    >>> p2 = np.polyder(p)
    >>> p2
    poly1d([3, 2, 1])

    which evaluates to:

    >>> p2(2.)
    17.0

    We can verify this, approximating the derivative with
    ``(f(x + h) - f(x))/h``:

    >>> (p(2. + 0.001) - p(2.)) / 0.001
    17.007000999997857

    The fourth-order derivative of a 3rd-order polynomial is zero:

    >>> np.polyder(p, 2)
    poly1d([6, 2])
    >>> np.polyder(p, 3)
    poly1d([6])
    >>> np.polyder(p, 4)
    poly1d([ 0.])
    
    """
    raise NotImplementedError


def polydiv(u, v):
    """Returns the quotient and remainder of polynomial division.

    The input arrays are the coefficients (including any coefficients
    equal to zero) of the "numerator" (dividend) and "denominator"
    (divisor) polynomials, respectively.

    Parameters
    ----------
    u : array_like or poly1d
        Dividend polynomial's coefficients.

    v : array_like or poly1d
        Divisor polynomial's coefficients.

    Returns
    -------
    q : ndarray
        Coefficients, including those equal to zero, of the quotient.
    r : ndarray
        Coefficients, including those equal to zero, of the remainder.

    See Also
    --------
    poly, polyadd, polyder, polydiv, polyfit, polyint, polymul, polysub,
    polyval

    Notes
    -----
    Both `u` and `v` must be 0-d or 1-d (ndim = 0 or 1), but `u.ndim` need
    not equal `v.ndim`. In other words, all four possible combinations -
    ``u.ndim = v.ndim = 0``, ``u.ndim = v.ndim = 1``,
    ``u.ndim = 1, v.ndim = 0``, and ``u.ndim = 0, v.ndim = 1`` - work.

    Examples
    --------
    .. math:: \frac{3x^2 + 5x + 2}{2x + 1} = 1.5x + 1.75, remainder 0.25

    >>> x = np.array([3.0, 5.0, 2.0])
    >>> y = np.array([2.0, 1.0])
    >>> np.polydiv(x, y)
    (array([ 1.5 ,  1.75]), array([ 0.25]))
    
    """
    raise NotImplementedError


def polyfit(x, y, deg, rcond=None, full=False):
    """Least squares polynomial fit.

    Fit a polynomial ``p(x) = p[0] * x**deg + ... + p[deg]`` of degree `deg`
    to points `(x, y)`. Returns a vector of coefficients `p` that minimises
    the squared error.

    Parameters
    ----------
    x : array_like, shape (M,)
        x-coordinates of the M sample points ``(x[i], y[i])``.
    y : array_like, shape (M,) or (M, K)
        y-coordinates of the sample points. Several data sets of sample
        points sharing the same x-coordinates can be fitted at once by
        passing in a 2D-array that contains one dataset per column.
    deg : int
        Degree of the fitting polynomial
    rcond : float, optional
        Relative condition number of the fit. Singular values smaller than this
        relative to the largest singular value will be ignored. The default
        value is len(x)*eps, where eps is the relative precision of the float
        type, about 2e-16 in most cases.
    full : bool, optional
        Switch determining nature of return value. When it is
        False (the default) just the coefficients are returned, when True
        diagnostic information from the singular value decomposition is also
        returned.

    Returns
    -------
    p : ndarray, shape (M,) or (M, K)
        Polynomial coefficients, highest power first.
        If `y` was 2-D, the coefficients for `k`-th data set are in ``p[:,k]``.

    residuals, rank, singular_values, rcond : present only if `full` = True
        Residuals of the least-squares fit, the effective rank of the scaled
        Vandermonde coefficient matrix, its singular values, and the specified
        value of `rcond`. For more details, see `linalg.lstsq`.

    Warns
    -----
    RankWarning
        The rank of the coefficient matrix in the least-squares fit is
        deficient. The warning is only raised if `full` = False.

        The warnings can be turned off by

        >>> import warnings
        >>> warnings.simplefilter('ignore', np.RankWarning)

    See Also
    --------
    polyval : Computes polynomial values.
    linalg.lstsq : Computes a least-squares fit.
    scipy.interpolate.UnivariateSpline : Computes spline fits.

    Notes
    -----
    The solution minimizes the squared error

    .. math ::
        E = \sum_{j=0}^k |p(x_j) - y_j|^2

    in the equations::

        x[0]**n * p[n] + ... + x[0] * p[1] + p[0] = y[0]
        x[1]**n * p[n] + ... + x[1] * p[1] + p[0] = y[1]
        ...
        x[k]**n * p[n] + ... + x[k] * p[1] + p[0] = y[k]

    The coefficient matrix of the coefficients `p` is a Vandermonde matrix.

    `polyfit` issues a `RankWarning` when the least-squares fit is badly
    conditioned. This implies that the best fit is not well-defined due
    to numerical error. The results may be improved by lowering the polynomial
    degree or by replacing `x` by `x` - `x`.mean(). The `rcond` parameter
    can also be set to a value smaller than its default, but the resulting
    fit may be spurious: including contributions from the small singular
    values can add numerical noise to the result.

    Note that fitting polynomial coefficients is inherently badly conditioned
    when the degree of the polynomial is large or the interval of sample points
    is badly centered. The quality of the fit should always be checked in these
    cases. When polynomial fits are not satisfactory, splines may be a good
    alternative.

    References
    ----------
    .. [1] Wikipedia, "Curve fitting",
           http://en.wikipedia.org/wiki/Curve_fitting
    .. [2] Wikipedia, "Polynomial interpolation",
           http://en.wikipedia.org/wiki/Polynomial_interpolation

    Examples
    --------
    >>> x = np.array([0.0, 1.0, 2.0, 3.0,  4.0,  5.0])
    >>> y = np.array([0.0, 0.8, 0.9, 0.1, -0.8, -1.0])
    >>> z = np.polyfit(x, y, 3)
    >>> z
    array([ 0.08703704, -0.81349206,  1.69312169, -0.03968254])

    It is convenient to use `poly1d` objects for dealing with polynomials:

    >>> p = np.poly1d(z)
    >>> p(0.5)
    0.6143849206349179
    >>> p(3.5)
    -0.34732142857143039
    >>> p(10)
    22.579365079365115

    High-order polynomials may oscillate wildly:

    >>> p30 = np.poly1d(np.polyfit(x, y, 30))
    /... RankWarning: Polyfit may be poorly conditioned...
    >>> p30(4)
    -0.80000000000000204
    >>> p30(5)
    -0.99999999999999445
    >>> p30(4.5)
    -0.10547061179440398

    Illustration:

    >>> import matplotlib.pyplot as plt
    >>> xp = np.linspace(-2, 6, 100)
    >>> plt.plot(x, y, '.', xp, p(xp), '-', xp, p30(xp), '--')
    [<matplotlib.lines.Line2D object at 0x...>, <matplotlib.lines.Line2D object at 0x...>, <matplotlib.lines.Line2D object at 0x...>]
    >>> plt.ylim(-2,2)
    (-2, 2)
    >>> plt.show()
    
    """
    raise NotImplementedError


def polyint(p, m=1, k=None):
    """Return an antiderivative (indefinite integral) of a polynomial.

    The returned order `m` antiderivative `P` of polynomial `p` satisfies
    :math:`\frac{d^m}{dx^m}P(x) = p(x)` and is defined up to `m - 1`
    integration constants `k`. The constants determine the low-order
    polynomial part

    .. math:: \frac{k_{m-1}}{0!} x^0 + \ldots + \frac{k_0}{(m-1)!}x^{m-1}

    of `P` so that :math:`P^{(j)}(0) = k_{m-j-1}`.

    Parameters
    ----------
    p : {array_like, poly1d}
        Polynomial to differentiate.
        A sequence is interpreted as polynomial coefficients, see `poly1d`.
    m : int, optional
        Order of the antiderivative. (Default: 1)
    k : {None, list of `m` scalars, scalar}, optional
        Integration constants. They are given in the order of integration:
        those corresponding to highest-order terms come first.

        If ``None`` (default), all constants are assumed to be zero.
        If `m = 1`, a single scalar can be given instead of a list.

    See Also
    --------
    polyder : derivative of a polynomial
    poly1d.integ : equivalent method

    Examples
    --------
    The defining property of the antiderivative:

    >>> p = np.poly1d([1,1,1])
    >>> P = np.polyint(p)
    >>> P
    poly1d([ 0.33333333,  0.5       ,  1.        ,  0.        ])
    >>> np.polyder(P) == p
    True

    The integration constants default to zero, but can be specified:

    >>> P = np.polyint(p, 3)
    >>> P(0)
    0.0
    >>> np.polyder(P)(0)
    0.0
    >>> np.polyder(P, 2)(0)
    0.0
    >>> P = np.polyint(p, 3, k=[6,5,3])
    >>> P
    poly1d([ 0.01666667,  0.04166667,  0.16666667,  3. ,  5. ,  3. ])

    Note that 3 = 6 / 2!, and that the constants are given in the order of
    integrations. Constant of the highest-order polynomial term comes first:

    >>> np.polyder(P, 2)(0)
    6.0
    >>> np.polyder(P, 1)(0)
    5.0
    >>> P(0)
    3.0
    
    """
    raise NotImplementedError


def polymul(a1, a2):
    """Find the product of two polynomials.

    Finds the polynomial resulting from the multiplication of the two input
    polynomials. Each input must be either a poly1d object or a 1D sequence
    of polynomial coefficients, from highest to lowest degree.

    Parameters
    ----------
    a1, a2 : array_like or poly1d object
        Input polynomials.

    Returns
    -------
    out : ndarray or poly1d object
        The polynomial resulting from the multiplication of the inputs. If
        either inputs is a poly1d object, then the output is also a poly1d
        object. Otherwise, it is a 1D array of polynomial coefficients from
        highest to lowest degree.

    See Also
    --------
    poly1d : A one-dimensional polynomial class.
    poly, polyadd, polyder, polydiv, polyfit, polyint, polysub,
    polyval

    Examples
    --------
    >>> np.polymul([1, 2, 3], [9, 5, 1])
    array([ 9, 23, 38, 17,  3])

    Using poly1d objects:

    >>> p1 = np.poly1d([1, 2, 3])
    >>> p2 = np.poly1d([9, 5, 1])
    >>> print p1
       2
    1 x + 2 x + 3
    >>> print p2
       2
    9 x + 5 x + 1
    >>> print np.polymul(p1, p2)
       4      3      2
    9 x + 23 x + 38 x + 17 x + 3
    
    """
    raise NotImplementedError


def polysub(a1, a2):
    """Difference (subtraction) of two polynomials.

    Given two polynomials `a1` and `a2`, returns ``a1 - a2``.
    `a1` and `a2` can be either array_like sequences of the polynomials'
    coefficients (including coefficients equal to zero), or `poly1d` objects.

    Parameters
    ----------
    a1, a2 : array_like or poly1d
        Minuend and subtrahend polynomials, respectively.

    Returns
    -------
    out : ndarray or poly1d
        Array or `poly1d` object of the difference polynomial's coefficients.

    See Also
    --------
    polyval, polydiv, polymul, polyadd

    Examples
    --------
    .. math:: (2 x^2 + 10 x - 2) - (3 x^2 + 10 x -4) = (-x^2 + 2)

    >>> np.polysub([2, 10, -2], [3, 10, -4])
    array([-1,  0,  2])
    
    """
    raise NotImplementedError


def polyval(p, x):
    """Evaluate a polynomial at specific values.

    If `p` is of length N, this function returns the value:

        ``p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1]``

    If `x` is a sequence, then `p(x)` is returned for each element of `x`.
    If `x` is another polynomial then the composite polynomial `p(x(t))`
    is returned.

    Parameters
    ----------
    p : array_like or poly1d object
       1D array of polynomial coefficients (including coefficients equal
       to zero) from highest degree to the constant term, or an
       instance of poly1d.
    x : array_like or poly1d object
       A number, a 1D array of numbers, or an instance of poly1d, "at"
       which to evaluate `p`.

    Returns
    -------
    values : ndarray or poly1d
       If `x` is a poly1d instance, the result is the composition of the two
       polynomials, i.e., `x` is "substituted" in `p` and the simplified
       result is returned. In addition, the type of `x` - array_like or
       poly1d - governs the type of the output: `x` array_like => `values`
       array_like, `x` a poly1d object => `values` is also.

    See Also
    --------
    poly1d: A polynomial class.

    Notes
    -----
    Horner's scheme [1]_ is used to evaluate the polynomial. Even so,
    for polynomials of high degree the values may be inaccurate due to
    rounding errors. Use carefully.

    References
    ----------
    .. [1] I. N. Bronshtein, K. A. Semendyayev, and K. A. Hirsch (Eng.
       trans. Ed.), *Handbook of Mathematics*, New York, Van Nostrand
       Reinhold Co., 1985, pg. 720.

    Examples
    --------
    >>> np.polyval([3,0,1], 5)  # 3 * 5**2 + 0 * 5**1 + 1
    76
    >>> np.polyval([3,0,1], np.poly1d(5))
    poly1d([ 76.])
    >>> np.polyval(np.poly1d([3,0,1]), 5)
    76
    >>> np.polyval(np.poly1d([3,0,1]), np.poly1d(5))
    poly1d([ 76.])
    
    """
    raise NotImplementedError


def ppmt(rate, per, nper, pv, fv=0.0, when='end'):
    """Not implemented. Compute the payment against loan principal.

    Parameters
    ----------
    rate : array_like
        Rate of interest (per period)
    per : array_like, int
        Amount paid against the loan changes.  The `per` is the period of
        interest.
    nper : array_like
        Number of compounding periods
    pv : array_like
        Present value
    fv : array_like, optional
        Future value
    when : {{'begin', 1}, {'end', 0}}, {string, int}
        When payments are due ('begin' (1) or 'end' (0))

    See Also
    --------
    pmt, pv, ipmt
    
    """
    raise NotImplementedError


def prod(a, axis=None, dtype=None, out=None):
    """Return the product of array elements over a given axis.

    Parameters
    ----------
    a : array_like
        Input data.
    axis : int, optional
        Axis over which the product is taken.  By default, the product
        of all elements is calculated.
    dtype : data-type, optional
        The data-type of the returned array, as well as of the accumulator
        in which the elements are multiplied.  By default, if `a` is of
        integer type, `dtype` is the default platform integer. (Note: if
        the type of `a` is unsigned, then so is `dtype`.)  Otherwise,
        the dtype is the same as that of `a`.
    out : ndarray, optional
        Alternative output array in which to place the result. It must have
        the same shape as the expected output, but the type of the
        output values will be cast if necessary.

    Returns
    -------
    product_along_axis : ndarray, see `dtype` parameter above.
        An array shaped as `a` but with the specified axis removed.
        Returns a reference to `out` if specified.

    See Also
    --------
    ndarray.prod : equivalent method
    numpy.doc.ufuncs : Section "Output arguments"

    Notes
    -----
    Arithmetic is modular when using integer types, and no error is
    raised on overflow.  That means that, on a 32-bit platform:

    >>> x = np.array([536870910, 536870910, 536870910, 536870910])
    >>> np.prod(x) #random
    16

    Examples
    --------
    By default, calculate the product of all elements:

    >>> np.prod([1.,2.])
    2.0

    Even when the input array is two-dimensional:

    >>> np.prod([[1.,2.],[3.,4.]])
    24.0

    But we can also specify the axis over which to multiply:

    >>> np.prod([[1.,2.],[3.,4.]], axis=1)
    array([  2.,  12.])

    If the type of `x` is unsigned, then the output type is
    the unsigned platform integer:

    >>> x = np.array([1, 2, 3], dtype=np.uint8)
    >>> np.prod(x).dtype == np.uint
    True

    If `x` is of a signed integer type, then the output type
    is the default platform integer:

    >>> x = np.array([1, 2, 3], dtype=np.int8)
    >>> np.prod(x).dtype == np.int
    True
    
    """
    raise NotImplementedError


def product(a, axis=None, dtype=None, out=None):
    """Return the product of array elements over a given axis.

    See Also
    --------
    prod : equivalent function; see for details.
    
    """
    raise NotImplementedError


def promote_types(type1, type2):
    """promote_types(type1, type2)

    Returns the data type with the smallest size and smallest scalar
    kind to which both ``type1`` and ``type2`` may be safely cast.
    The returned data type is always in native byte order.

    This function is symmetric and associative.

    Parameters
    ----------
    type1 : dtype or dtype specifier
        First data type.
    type2 : dtype or dtype specifier
        Second data type.

    Returns
    -------
    out : dtype
        The promoted data type.

    Notes
    -----
    .. versionadded:: 1.6.0

    See Also
    --------
    result_type, dtype, can_cast

    Examples
    --------
    >>> np.promote_types('f4', 'f8')
    dtype('float64')

    >>> np.promote_types('i8', 'f4')
    dtype('float64')

    >>> np.promote_types('>i8', '<c8')
    dtype('complex128')

    >>> np.promote_types('i1', 'S8')
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: invalid type promotion
    
    """
    # BUILTIN
    raise NotImplementedError


def ptp(a, axis=None, out=None):
    """Range of values (maximum - minimum) along an axis.

    The name of the function comes from the acronym for 'peak to peak'.

    Parameters
    ----------
    a : array_like
        Input values.
    axis : int, optional
        Axis along which to find the peaks.  By default, flatten the
        array.
    out : array_like
        Alternative output array in which to place the result. It must
        have the same shape and buffer length as the expected output,
        but the type of the output values will be cast if necessary.

    Returns
    -------
    ptp : ndarray
        A new array holding the result, unless `out` was
        specified, in which case a reference to `out` is returned.

    Examples
    --------
    >>> x = np.arange(4).reshape((2,2))
    >>> x
    array([[0, 1],
           [2, 3]])

    >>> np.ptp(x, axis=0)
    array([2, 2])

    >>> np.ptp(x, axis=1)
    array([1, 1])
    
    """
    raise NotImplementedError


def put(a, ind, v, mode='raise'):
    """Replaces specified elements of an array with given values.

    The indexing works on the flattened target array. `put` is roughly
    equivalent to:

    ::

        a.flat[ind] = v

    Parameters
    ----------
    a : ndarray
        Target array.
    ind : array_like
        Target indices, interpreted as integers.
    v : array_like
        Values to place in `a` at target indices. If `v` is shorter than
        `ind` it will be repeated as necessary.
    mode : {'raise', 'wrap', 'clip'}, optional
        Specifies how out-of-bounds indices will behave.

        * 'raise' -- raise an error (default)
        * 'wrap' -- wrap around
        * 'clip' -- clip to the range

        'clip' mode means that all indices that are too large are replaced
        by the index that addresses the last element along that axis. Note
        that this disables indexing with negative numbers.

    See Also
    --------
    putmask, place

    Examples
    --------
    >>> a = np.arange(5)
    >>> np.put(a, [0, 2], [-44, -55])
    >>> a
    array([-44,   1, -55,   3,   4])

    >>> a = np.arange(5)
    >>> np.put(a, 22, -5, mode='clip')
    >>> a
    array([ 0,  1,  2,  3, -5])
    
    """
    raise NotImplementedError


def putmask(a, mask, values):
    """putmask(a, mask, values)

    Changes elements of an array based on conditional and input values.

    Sets ``a.flat[n] = values[n]`` for each n where ``mask.flat[n]==True``.

    If `values` is not the same size as `a` and `mask` then it will repeat.
    This gives behavior different from ``a[mask] = values``.

    Parameters
    ----------
    a : array_like
        Target array.
    mask : array_like
        Boolean mask array. It has to be the same shape as `a`.
    values : array_like
        Values to put into `a` where `mask` is True. If `values` is smaller
        than `a` it will be repeated.

    See Also
    --------
    place, put, take

    Examples
    --------
    >>> x = np.arange(6).reshape(2, 3)
    >>> np.putmask(x, x>2, x**2)
    >>> x
    array([[ 0,  1,  2],
           [ 9, 16, 25]])

    If `values` is smaller than `a` it is repeated:

    >>> x = np.arange(5)
    >>> np.putmask(x, x>1, [-33, -44])
    >>> x
    array([  0,   1, -33, -44, -33])
    
    """
    # BUILTIN
    raise NotImplementedError


def pv(rate, nper, pmt, fv=0.0, when='end'):
    """Compute the present value.

    Given:
     * a future value, `fv`
     * an interest `rate` compounded once per period, of which
       there are
     * `nper` total
     * a (fixed) payment, `pmt`, paid either
     * at the beginning (`when` = {'begin', 1}) or the end
       (`when` = {'end', 0}) of each period

    Return:
       the value now

    Parameters
    ----------
    rate : array_like
        Rate of interest (per period)
    nper : array_like
        Number of compounding periods
    pmt : array_like
        Payment
    fv : array_like, optional
        Future value
    when : {{'begin', 1}, {'end', 0}}, {string, int}, optional
        When payments are due ('begin' (1) or 'end' (0))

    Returns
    -------
    out : ndarray, float
        Present value of a series of payments or investments.

    Notes
    -----
    The present value is computed by solving the equation::

     fv +
     pv*(1 + rate)**nper +
     pmt*(1 + rate*when)/rate*((1 + rate)**nper - 1) = 0

    or, when ``rate = 0``::

     fv + pv + pmt * nper = 0

    for `pv`, which is then returned.

    References
    ----------
    .. [WRW] Wheeler, D. A., E. Rathke, and R. Weir (Eds.) (2009, May).
       Open Document Format for Office Applications (OpenDocument)v1.2,
       Part 2: Recalculated Formula (OpenFormula) Format - Annotated Version,
       Pre-Draft 12. Organization for the Advancement of Structured Information
       Standards (OASIS). Billerica, MA, USA. [ODT Document].
       Available:
       http://www.oasis-open.org/committees/documents.php?wg_abbrev=office-formula
       OpenDocument-formula-20090508.odt

    Examples
    --------
    What is the present value (e.g., the initial investment)
    of an investment that needs to total $15692.93
    after 10 years of saving $100 every month?  Assume the
    interest rate is 5% (annually) compounded monthly.

    >>> np.pv(0.05/12, 10*12, -100, 15692.93)
    -100.00067131625819

    By convention, the negative sign represents cash flow out
    (i.e., money not available today).  Thus, to end up with
    $15,692.93 in 10 years saving $100 a month at 5% annual
    interest, one's initial deposit should also be $100.

    If any input is array_like, ``pv`` returns an array of equal shape.
    Let's compare different interest rates in the example above:

    >>> a = np.array((0.05, 0.04, 0.03))/12
    >>> np.pv(a, 10*12, -100, 15692.93)
    array([ -100.00067132,  -649.26771385, -1273.78633713])

    So, to end up with the same $15692.93 under the same $100 per month
    "savings plan," for annual interest rates of 4% and 3%, one would
    need initial investments of $649.27 and $1273.79, respectively.
    
    """
    raise NotImplementedError


def rank(a):
    """Return the number of dimensions of an array.

    If `a` is not already an array, a conversion is attempted.
    Scalars are zero dimensional.

    Parameters
    ----------
    a : array_like
        Array whose number of dimensions is desired. If `a` is not an array,
        a conversion is attempted.

    Returns
    -------
    number_of_dimensions : int
        The number of dimensions in the array.

    See Also
    --------
    ndim : equivalent function
    ndarray.ndim : equivalent property
    shape : dimensions of array
    ndarray.shape : dimensions of array

    Notes
    -----
    In the old Numeric package, `rank` was the term used for the number of
    dimensions, but in Numpy `ndim` is used instead.

    Examples
    --------
    >>> np.rank([1,2,3])
    1
    >>> np.rank(np.array([[1,2,3],[4,5,6]]))
    2
    >>> np.rank(1)
    0
    
    """
    raise NotImplementedError


def rate(nper, pmt, pv, fv, when='end', guess=0.1, tol=1e-06, maxiter=100):
    """Compute the rate of interest per period.

    Parameters
    ----------
    nper : array_like
        Number of compounding periods
    pmt : array_like
        Payment
    pv : array_like
        Present value
    fv : array_like
        Future value
    when : {{'begin', 1}, {'end', 0}}, {string, int}, optional
        When payments are due ('begin' (1) or 'end' (0))
    guess : float, optional
        Starting guess for solving the rate of interest
    tol : float, optional
        Required tolerance for the solution
    maxiter : int, optional
        Maximum iterations in finding the solution

    Notes
    -----
    The rate of interest is computed by iteratively solving the
    (non-linear) equation::

     fv + pv*(1+rate)**nper + pmt*(1+rate*when)/rate * ((1+rate)**nper - 1) = 0

    for ``rate``.

    References
    ----------
    Wheeler, D. A., E. Rathke, and R. Weir (Eds.) (2009, May). Open Document
    Format for Office Applications (OpenDocument)v1.2, Part 2: Recalculated
    Formula (OpenFormula) Format - Annotated Version, Pre-Draft 12.
    Organization for the Advancement of Structured Information Standards
    (OASIS). Billerica, MA, USA. [ODT Document]. Available:
    http://www.oasis-open.org/committees/documents.php?wg_abbrev=office-formula
    OpenDocument-formula-20090508.odt
    
    """
    raise NotImplementedError


def ravel(a, order='C'):
    """Return a flattened array.

    A 1-D array, containing the elements of the input, is returned.  A copy is
    made only if needed.

    Parameters
    ----------
    a : array_like
        Input array.  The elements in ``a`` are read in the order specified by
        `order`, and packed as a 1-D array.
    order : {'C','F', 'A', 'K'}, optional
        The elements of ``a`` are read in this order. 'C' means to view
        the elements in C (row-major) order. 'F' means to view the elements
        in Fortran (column-major) order. 'A' means to view the elements
        in 'F' order if a is Fortran contiguous, 'C' order otherwise.
        'K' means to view the elements in the order they occur in memory,
        except for reversing the data when strides are negative.
        By default, 'C' order is used.

    Returns
    -------
    1d_array : ndarray
        Output of the same dtype as `a`, and of shape ``(a.size(),)``.

    See Also
    --------
    ndarray.flat : 1-D iterator over an array.
    ndarray.flatten : 1-D array copy of the elements of an array
                      in row-major order.

    Notes
    -----
    In row-major order, the row index varies the slowest, and the column
    index the quickest.  This can be generalized to multiple dimensions,
    where row-major order implies that the index along the first axis
    varies slowest, and the index along the last quickest.  The opposite holds
    for Fortran-, or column-major, mode.

    Examples
    --------
    It is equivalent to ``reshape(-1, order=order)``.

    >>> x = np.array([[1, 2, 3], [4, 5, 6]])
    >>> print np.ravel(x)
    [1 2 3 4 5 6]

    >>> print x.reshape(-1)
    [1 2 3 4 5 6]

    >>> print np.ravel(x, order='F')
    [1 4 2 5 3 6]

    When ``order`` is 'A', it will preserve the array's 'C' or 'F' ordering:

    >>> print np.ravel(x.T)
    [1 4 2 5 3 6]
    >>> print np.ravel(x.T, order='A')
    [1 2 3 4 5 6]

    When ``order`` is 'K', it will preserve orderings that are neither 'C'
    nor 'F', but won't reverse axes:

    >>> a = np.arange(3)[::-1]; a
    array([2, 1, 0])
    >>> a.ravel(order='C')
    array([2, 1, 0])
    >>> a.ravel(order='K')
    array([2, 1, 0])

    >>> a = np.arange(12).reshape(2,3,2).swapaxes(1,2); a
    array([[[ 0,  2,  4],
            [ 1,  3,  5]],
           [[ 6,  8, 10],
            [ 7,  9, 11]]])
    >>> a.ravel(order='C')
    array([ 0,  2,  4,  1,  3,  5,  6,  8, 10,  7,  9, 11])
    >>> a.ravel(order='K')
    array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11])
    
    """
    raise NotImplementedError


def ravel_multi_index(multi_index, dims, mode='raise', order='C'):
    """ravel_multi_index(multi_index, dims, mode='raise', order='C')

    Converts a tuple of index arrays into an array of flat
    indices, applying boundary modes to the multi-index.

    Parameters
    ----------
    multi_index : tuple of array_like
        A tuple of integer arrays, one array for each dimension.
    dims : tuple of ints
        The shape of array into which the indices from ``multi_index`` apply.
    mode : {'raise', 'wrap', 'clip'}, optional
        Specifies how out-of-bounds indices are handled.  Can specify
        either one mode or a tuple of modes, one mode per index.

        * 'raise' -- raise an error (default)
        * 'wrap' -- wrap around
        * 'clip' -- clip to the range

        In 'clip' mode, a negative index which would normally
        wrap will clip to 0 instead.
    order : {'C', 'F'}, optional
        Determines whether the multi-index should be viewed as indexing in
        C (row-major) order or FORTRAN (column-major) order.

    Returns
    -------
    raveled_indices : ndarray
        An array of indices into the flattened version of an array
        of dimensions ``dims``.

    See Also
    --------
    unravel_index

    Notes
    -----
    .. versionadded:: 1.6.0

    Examples
    --------
    >>> arr = np.array([[3,6,6],[4,5,1]])
    >>> np.ravel_multi_index(arr, (7,6))
    array([22, 41, 37])
    >>> np.ravel_multi_index(arr, (7,6), order='F')
    array([31, 41, 13])
    >>> np.ravel_multi_index(arr, (4,6), mode='clip')
    array([22, 23, 19])
    >>> np.ravel_multi_index(arr, (4,4), mode=('clip','wrap'))
    array([12, 13, 13])

    >>> np.ravel_multi_index((3,1,4,1), (6,7,8,9))
    1621
    
    """
    # BUILTIN
    raise NotImplementedError


def real(val):
    """Return the real part of the elements of the array.

    Parameters
    ----------
    val : array_like
        Input array.

    Returns
    -------
    out : ndarray
        Output array. If `val` is real, the type of `val` is used for the
        output.  If `val` has complex elements, the returned type is float.

    See Also
    --------
    real_if_close, imag, angle

    Examples
    --------
    >>> a = np.array([1+2j, 3+4j, 5+6j])
    >>> a.real
    array([ 1.,  3.,  5.])
    >>> a.real = 9
    >>> a
    array([ 9.+2.j,  9.+4.j,  9.+6.j])
    >>> a.real = np.array([9, 8, 7])
    >>> a
    array([ 9.+2.j,  8.+4.j,  7.+6.j])
    
    """
    raise NotImplementedError


def real_if_close(a, tol=100):
    """If complex input returns a real array if complex parts are close to zero.

    "Close to zero" is defined as `tol` * (machine epsilon of the type for
    `a`).

    Parameters
    ----------
    a : array_like
        Input array.
    tol : float
        Tolerance in machine epsilons for the complex part of the elements
        in the array.

    Returns
    -------
    out : ndarray
        If `a` is real, the type of `a` is used for the output.  If `a`
        has complex elements, the returned type is float.

    See Also
    --------
    real, imag, angle

    Notes
    -----
    Machine epsilon varies from machine to machine and between data types
    but Python floats on most platforms have a machine epsilon equal to
    2.2204460492503131e-16.  You can use 'np.finfo(np.float).eps' to print
    out the machine epsilon for floats.

    Examples
    --------
    >>> np.finfo(np.float).eps
    2.2204460492503131e-16

    >>> np.real_if_close([2.1 + 4e-14j], tol=1000)
    array([ 2.1])
    >>> np.real_if_close([2.1 + 4e-13j], tol=1000)
    array([ 2.1 +4.00000000e-13j])
    
    """
    raise NotImplementedError


def recfromcsv(fname):
    """Load ASCII data stored in a comma-separated file.

    The returned array is a record array (if ``usemask=False``, see
    `recarray`) or a masked record array (if ``usemask=True``,
    see `ma.mrecords.MaskedRecords`).

    For a complete description of all the input parameters, see `genfromtxt`.

    See Also
    --------
    numpy.genfromtxt : generic function to load ASCII data.
    
    """
    raise NotImplementedError


def recfromtxt(fname):
    """Load ASCII data from a file and return it in a record array.

    If ``usemask=False`` a standard `recarray` is returned,
    if ``usemask=True`` a MaskedRecords array is returned.

    Complete description of all the optional input parameters is available in
    the docstring of the `genfromtxt` function.

    See Also
    --------
    numpy.genfromtxt : generic function

    Notes
    -----
    By default, `dtype` is None, which means that the data-type of the output
    array will be determined from the data.
    
    """
    raise NotImplementedError


def repeat(a, repeats, axis=None):
    """Repeat elements of an array.

    Parameters
    ----------
    a : array_like
        Input array.
    repeats : {int, array of ints}
        The number of repetitions for each element.  `repeats` is broadcasted
        to fit the shape of the given axis.
    axis : int, optional
        The axis along which to repeat values.  By default, use the
        flattened input array, and return a flat output array.

    Returns
    -------
    repeated_array : ndarray
        Output array which has the same shape as `a`, except along
        the given axis.

    See Also
    --------
    tile : Tile an array.

    Examples
    --------
    >>> x = np.array([[1,2],[3,4]])
    >>> np.repeat(x, 2)
    array([1, 1, 2, 2, 3, 3, 4, 4])
    >>> np.repeat(x, 3, axis=1)
    array([[1, 1, 1, 2, 2, 2],
           [3, 3, 3, 4, 4, 4]])
    >>> np.repeat(x, [1, 2], axis=0)
    array([[1, 2],
           [3, 4],
           [3, 4]])
    
    """
    raise NotImplementedError


def require(a, dtype=None, requirements=None):
    """Return an ndarray of the provided type that satisfies requirements.

    This function is useful to be sure that an array with the correct flags
    is returned for passing to compiled code (perhaps through ctypes).

    Parameters
    ----------
    a : array_like
       The object to be converted to a type-and-requirement-satisfying array.
    dtype : data-type
       The required data-type, the default data-type is float64).
    requirements : str or list of str
       The requirements list can be any of the following

       * 'F_CONTIGUOUS' ('F') - ensure a Fortran-contiguous array
       * 'C_CONTIGUOUS' ('C') - ensure a C-contiguous array
       * 'ALIGNED' ('A')      - ensure a data-type aligned array
       * 'WRITEABLE' ('W')    - ensure a writable array
       * 'OWNDATA' ('O')      - ensure an array that owns its own data

    See Also
    --------
    asarray : Convert input to an ndarray.
    asanyarray : Convert to an ndarray, but pass through ndarray subclasses.
    ascontiguousarray : Convert input to a contiguous array.
    asfortranarray : Convert input to an ndarray with column-major
                     memory order.
    ndarray.flags : Information about the memory layout of the array.

    Notes
    -----
    The returned array will be guaranteed to have the listed requirements
    by making a copy if needed.

    Examples
    --------
    >>> x = np.arange(6).reshape(2,3)
    >>> x.flags
      C_CONTIGUOUS : True
      F_CONTIGUOUS : False
      OWNDATA : False
      WRITEABLE : True
      ALIGNED : True
      UPDATEIFCOPY : False

    >>> y = np.require(x, dtype=np.float32, requirements=['A', 'O', 'W', 'F'])
    >>> y.flags
      C_CONTIGUOUS : False
      F_CONTIGUOUS : True
      OWNDATA : True
      WRITEABLE : True
      ALIGNED : True
      UPDATEIFCOPY : False
    
    """
    raise NotImplementedError


def reshape(a, newshape, order='C'):
    """Gives a new shape to an array without changing its data.

    Parameters
    ----------
    a : array_like
        Array to be reshaped.
    newshape : int or tuple of ints
        The new shape should be compatible with the original shape. If
        an integer, then the result will be a 1-D array of that length.
        One shape dimension can be -1. In this case, the value is inferred
        from the length of the array and remaining dimensions.
    order : {'C', 'F', 'A'}, optional
        Determines whether the array data should be viewed as in C
        (row-major) order, FORTRAN (column-major) order, or the C/FORTRAN
        order should be preserved.

    Returns
    -------
    reshaped_array : ndarray
        This will be a new view object if possible; otherwise, it will
        be a copy.


    See Also
    --------
    ndarray.reshape : Equivalent method.

    Notes
    -----

    It is not always possible to change the shape of an array without
    copying the data. If you want an error to be raise if the data is copied,
    you should assign the new shape to the shape attribute of the array::

     >>> a = np.zeros((10, 2))
     # A transpose make the array non-contiguous
     >>> b = a.T
     # Taking a view makes it possible to modify the shape without modiying the
     # initial object.
     >>> c = b.view()
     >>> c.shape = (20)
     AttributeError: incompatible shape for a non-contiguous array


    Examples
    --------
    >>> a = np.array([[1,2,3], [4,5,6]])
    >>> np.reshape(a, 6)
    array([1, 2, 3, 4, 5, 6])
    >>> np.reshape(a, 6, order='F')
    array([1, 4, 2, 5, 3, 6])

    >>> np.reshape(a, (3,-1))       # the unspecified value is inferred to be 2
    array([[1, 2],
           [3, 4],
           [5, 6]])
    
    """
    raise NotImplementedError


def resize(a, new_shape):
    """Return a new array with the specified shape.

    If the new array is larger than the original array, then the new
    array is filled with repeated copies of `a`.  Note that this behavior
    is different from a.resize(new_shape) which fills with zeros instead
    of repeated copies of `a`.

    Parameters
    ----------
    a : array_like
        Array to be resized.

    new_shape : int or tuple of int
        Shape of resized array.

    Returns
    -------
    reshaped_array : ndarray
        The new array is formed from the data in the old array, repeated
        if necessary to fill out the required number of elements.  The
        data are repeated in the order that they are stored in memory.

    See Also
    --------
    ndarray.resize : resize an array in-place.

    Examples
    --------
    >>> a=np.array([[0,1],[2,3]])
    >>> np.resize(a,(1,4))
    array([[0, 1, 2, 3]])
    >>> np.resize(a,(2,4))
    array([[0, 1, 2, 3],
           [0, 1, 2, 3]])
    
    """
    raise NotImplementedError


def restoredot():
    """Restore `dot`, `vdot`, and `innerproduct` to the default non-BLAS
    implementations.

    Typically, the user will only need to call this when troubleshooting and
    installation problem, reproducing the conditions of a build without an
    accelerated BLAS, or when being very careful about benchmarking linear
    algebra operations.

    See Also
    --------
    alterdot : `restoredot` undoes the effects of `alterdot`.
    
    """
    # BUILTIN
    raise NotImplementedError


def result_type(*arrays_and_dtypes):
    """result_type(*arrays_and_dtypes)

    Returns the type that results from applying the NumPy
    type promotion rules to the arguments.

    Type promotion in NumPy works similarly to the rules in languages
    like C++, with some slight differences.  When both scalars and
    arrays are used, the array's type takes precedence and the actual value
    of the scalar is taken into account.

    For example, calculating 3*a, where a is an array of 32-bit floats,
    intuitively should result in a 32-bit float output.  If the 3 is a
    32-bit integer, the NumPy rules indicate it can't convert losslessly
    into a 32-bit float, so a 64-bit float should be the result type.
    By examining the value of the constant, '3', we see that it fits in
    an 8-bit integer, which can be cast losslessly into the 32-bit float.

    Parameters
    ----------
    arrays_and_dtypes : list of arrays and dtypes
        The operands of some operation whose result type is needed.

    Returns
    -------
    out : dtype
        The result type.

    See also
    --------
    dtype, promote_types, min_scalar_type, can_cast

    Notes
    -----
    .. versionadded:: 1.6.0

    The specific algorithm used is as follows.

    Categories are determined by first checking which of boolean,
    integer (int/uint), or floating point (float/complex) the maximum
    kind of all the arrays and the scalars are.
    
    If there are only scalars or the maximum category of the scalars
    is higher than the maximum category of the arrays,
    the data types are combined with :func:`promote_types`
    to produce the return value.

    Otherwise, `min_scalar_type` is called on each array, and
    the resulting data types are all combined with :func:`promote_types`
    to produce the return value.

    The set of int values is not a subset of the uint values for types
    with the same number of bits, something not reflected in
    :func:`min_scalar_type`, but handled as a special case in `result_type`.

    Examples
    --------
    >>> np.result_type(3, np.arange(7, dtype='i1'))
    dtype('int8')

    >>> np.result_type('i4', 'c8')
    dtype('complex128')

    >>> np.result_type(3.0, -2)
    dtype('float64')
    
    """
    # BUILTIN
    raise NotImplementedError


def roll(a, shift, axis=None):
    """Roll array elements along a given axis.

    Elements that roll beyond the last position are re-introduced at
    the first.

    Parameters
    ----------
    a : array_like
        Input array.
    shift : int
        The number of places by which elements are shifted.
    axis : int, optional
        The axis along which elements are shifted.  By default, the array
        is flattened before shifting, after which the original
        shape is restored.

    Returns
    -------
    res : ndarray
        Output array, with the same shape as `a`.

    See Also
    --------
    rollaxis : Roll the specified axis backwards, until it lies in a
               given position.

    Examples
    --------
    >>> x = np.arange(10)
    >>> np.roll(x, 2)
    array([8, 9, 0, 1, 2, 3, 4, 5, 6, 7])

    >>> x2 = np.reshape(x, (2,5))
    >>> x2
    array([[0, 1, 2, 3, 4],
           [5, 6, 7, 8, 9]])
    >>> np.roll(x2, 1)
    array([[9, 0, 1, 2, 3],
           [4, 5, 6, 7, 8]])
    >>> np.roll(x2, 1, axis=0)
    array([[5, 6, 7, 8, 9],
           [0, 1, 2, 3, 4]])
    >>> np.roll(x2, 1, axis=1)
    array([[4, 0, 1, 2, 3],
           [9, 5, 6, 7, 8]])
    
    """
    raise NotImplementedError


def rollaxis(a, axis, start=0):
    """Roll the specified axis backwards, until it lies in a given position.

    Parameters
    ----------
    a : ndarray
        Input array.
    axis : int
        The axis to roll backwards.  The positions of the other axes do not
        change relative to one another.
    start : int, optional
        The axis is rolled until it lies before this position.  The default,
        0, results in a "complete" roll.

    Returns
    -------
    res : ndarray
        Output array.

    See Also
    --------
    roll : Roll the elements of an array by a number of positions along a
        given axis.

    Examples
    --------
    >>> a = np.ones((3,4,5,6))
    >>> np.rollaxis(a, 3, 1).shape
    (3, 6, 4, 5)
    >>> np.rollaxis(a, 2).shape
    (5, 3, 4, 6)
    >>> np.rollaxis(a, 1, 4).shape
    (3, 5, 6, 4)
    
    """
    raise NotImplementedError


def roots(p):
    """Return the roots of a polynomial with coefficients given in p.

    The values in the rank-1 array `p` are coefficients of a polynomial.
    If the length of `p` is n+1 then the polynomial is described by::

      p[0] * x**n + p[1] * x**(n-1) + ... + p[n-1]*x + p[n]

    Parameters
    ----------
    p : array_like
        Rank-1 array of polynomial coefficients.

    Returns
    -------
    out : ndarray
        An array containing the complex roots of the polynomial.

    Raises
    ------
    ValueError :
        When `p` cannot be converted to a rank-1 array.

    See also
    --------
    poly : Find the coefficients of a polynomial with a given sequence
           of roots.
    polyval : Evaluate a polynomial at a point.
    polyfit : Least squares polynomial fit.
    poly1d : A one-dimensional polynomial class.

    Notes
    -----
    The algorithm relies on computing the eigenvalues of the
    companion matrix [1]_.

    References
    ----------
    .. [1] R. A. Horn & C. R. Johnson, *Matrix Analysis*.  Cambridge, UK:
        Cambridge University Press, 1999, pp. 146-7.

    Examples
    --------
    >>> coeff = [3.2, 2, 1]
    >>> np.roots(coeff)
    array([-0.3125+0.46351241j, -0.3125-0.46351241j])
    
    """
    raise NotImplementedError


def rot90(m, k=1):
    """Rotate an array by 90 degrees in the counter-clockwise direction.

    The first two dimensions are rotated; therefore, the array must be at
    least 2-D.

    Parameters
    ----------
    m : array_like
        Array of two or more dimensions.
    k : integer
        Number of times the array is rotated by 90 degrees.

    Returns
    -------
    y : ndarray
        Rotated array.

    See Also
    --------
    fliplr : Flip an array horizontally.
    flipud : Flip an array vertically.

    Examples
    --------
    >>> m = np.array([[1,2],[3,4]], int)
    >>> m
    array([[1, 2],
           [3, 4]])
    >>> np.rot90(m)
    array([[2, 4],
           [1, 3]])
    >>> np.rot90(m, 2)
    array([[4, 3],
           [2, 1]])
    
    """
    raise NotImplementedError


def round(a, decimals=0, out=None):
    """Round an array to the given number of decimals.

    Refer to `around` for full documentation.

    See Also
    --------
    around : equivalent function
    
    """
    raise NotImplementedError


def round_(a, decimals=0, out=None):
    """Round an array to the given number of decimals.

    Refer to `around` for full documentation.

    See Also
    --------
    around : equivalent function
    
    """
    raise NotImplementedError


def row_stack(tup):
    """Stack arrays in sequence vertically (row wise).

    Take a sequence of arrays and stack them vertically to make a single
    array. Rebuild arrays divided by `vsplit`.

    Parameters
    ----------
    tup : sequence of ndarrays
        Tuple containing arrays to be stacked. The arrays must have the same
        shape along all but the first axis.

    Returns
    -------
    stacked : ndarray
        The array formed by stacking the given arrays.

    See Also
    --------
    hstack : Stack arrays in sequence horizontally (column wise).
    dstack : Stack arrays in sequence depth wise (along third dimension).
    concatenate : Join a sequence of arrays together.
    vsplit : Split array into a list of multiple sub-arrays vertically.


    Notes
    -----
    Equivalent to ``np.concatenate(tup, axis=0)``

    Examples
    --------
    >>> a = np.array([1, 2, 3])
    >>> b = np.array([2, 3, 4])
    >>> np.vstack((a,b))
    array([[1, 2, 3],
           [2, 3, 4]])

    >>> a = np.array([[1], [2], [3]])
    >>> b = np.array([[2], [3], [4]])
    >>> np.vstack((a,b))
    array([[1],
           [2],
           [3],
           [2],
           [3],
           [4]])
    
    """
    raise NotImplementedError


def safe_eval(source):
    """Protected string evaluation.

    Evaluate a string containing a Python literal expression without
    allowing the execution of arbitrary non-literal code.

    Parameters
    ----------
    source : str
        The string to evaluate.

    Returns
    -------
    obj : object
       The result of evaluating `source`.

    Raises
    ------
    SyntaxError
        If the code has invalid Python syntax, or if it contains non-literal
        code.

    Examples
    --------
    >>> np.safe_eval('1')
    1
    >>> np.safe_eval('[1, 2, 3]')
    [1, 2, 3]
    >>> np.safe_eval('{"foo": ("bar", 10.0)}')
    {'foo': ('bar', 10.0)}

    >>> np.safe_eval('import os')
    Traceback (most recent call last):
      ...
    SyntaxError: invalid syntax

    >>> np.safe_eval('open("/home/user/.ssh/id_dsa").read()')
    Traceback (most recent call last):
      ...
    SyntaxError: Unsupported source construct: compiler.ast.CallFunc
    
    """
    raise NotImplementedError


def save(file, arr):
    """Save an array to a binary file in NumPy ``.npy`` format.

    Parameters
    ----------
    file : file or str
        File or filename to which the data is saved.  If file is a file-object,
        then the filename is unchanged.  If file is a string, a ``.npy``
        extension will be appended to the file name if it does not already
        have one.
    arr : array_like
        Array data to be saved.

    See Also
    --------
    savez : Save several arrays into a ``.npz`` archive
    savetxt, load

    Notes
    -----
    For a description of the ``.npy`` format, see `format`.

    Examples
    --------
    >>> from tempfile import TemporaryFile
    >>> outfile = TemporaryFile()

    >>> x = np.arange(10)
    >>> np.save(outfile, x)

    >>> outfile.seek(0) # Only needed here to simulate closing & reopening file
    >>> np.load(outfile)
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    
    """
    raise NotImplementedError


def savetxt(fname, X, fmt='%.18e', delimiter=' ', newline='\n'):
    """Save an array to a text file.

    Parameters
    ----------
    fname : filename or file handle
        If the filename ends in ``.gz``, the file is automatically saved in
        compressed gzip format.  `loadtxt` understands gzipped files
        transparently.
    X : array_like
        Data to be saved to a text file.
    fmt : str or sequence of strs
        A single format (%10.5f), a sequence of formats, or a
        multi-format string, e.g. 'Iteration %d -- %10.5f', in which
        case `delimiter` is ignored.
    delimiter : str
        Character separating columns.
    newline : str
        .. versionadded:: 1.5.0

        Character separating lines.


    See Also
    --------
    save : Save an array to a binary file in NumPy ``.npy`` format
    savez : Save several arrays into a ``.npz`` compressed archive

    Notes
    -----
    Further explanation of the `fmt` parameter
    (``%[flag]width[.precision]specifier``):

    flags:
        ``-`` : left justify

        ``+`` : Forces to preceed result with + or -.

        ``0`` : Left pad the number with zeros instead of space (see width).

    width:
        Minimum number of characters to be printed. The value is not truncated
        if it has more characters.

    precision:
        - For integer specifiers (eg. ``d,i,o,x``), the minimum number of
          digits.
        - For ``e, E`` and ``f`` specifiers, the number of digits to print
          after the decimal point.
        - For ``g`` and ``G``, the maximum number of significant digits.
        - For ``s``, the maximum number of characters.

    specifiers:
        ``c`` : character

        ``d`` or ``i`` : signed decimal integer

        ``e`` or ``E`` : scientific notation with ``e`` or ``E``.

        ``f`` : decimal floating point

        ``g,G`` : use the shorter of ``e,E`` or ``f``

        ``o`` : signed octal

        ``s`` : string of characters

        ``u`` : unsigned decimal integer

        ``x,X`` : unsigned hexadecimal integer

    This explanation of ``fmt`` is not complete, for an exhaustive
    specification see [1]_.

    References
    ----------
    .. [1] `Format Specification Mini-Language
           <http://docs.python.org/library/string.html#
           format-specification-mini-language>`_, Python Documentation.

    Examples
    --------
    >>> x = y = z = np.arange(0.0,5.0,1.0)
    >>> np.savetxt('test.out', x, delimiter=',')   # X is an array
    >>> np.savetxt('test.out', (x,y,z))   # x,y,z equal sized 1D arrays
    >>> np.savetxt('test.out', x, fmt='%1.4e')   # use exponential notation
    
    """
    raise NotImplementedError


def savez(file):
    """Save several arrays into a single file in uncompressed ``.npz`` format.

    If arguments are passed in with no keywords, the corresponding variable
    names, in the .npz file, are 'arr_0', 'arr_1', etc. If keyword arguments
    are given, the corresponding variable names, in the ``.npz`` file will
    match the keyword names.

    Parameters
    ----------
    file : str or file
        Either the file name (string) or an open file (file-like object)
        where the data will be saved. If file is a string, the ``.npz``
        extension will be appended to the file name if it is not already there.
    *args : Arguments, optional
        Arrays to save to the file. Since it is not possible for Python to
        know the names of the arrays outside `savez`, the arrays will be saved
        with names "arr_0", "arr_1", and so on. These arguments can be any
        expression.
    **kwds : Keyword arguments, optional
        Arrays to save to the file. Arrays will be saved in the file with the
        keyword names.

    Returns
    -------
    None

    See Also
    --------
    save : Save a single array to a binary file in NumPy format.
    savetxt : Save an array to a file as plain text.

    Notes
    -----
    The ``.npz`` file format is a zipped archive of files named after the
    variables they contain.  The archive is not compressed and each file
    in the archive contains one variable in ``.npy`` format. For a
    description of the ``.npy`` format, see `format`.

    When opening the saved ``.npz`` file with `load` a `NpzFile` object is
    returned. This is a dictionary-like object which can be queried for
    its list of arrays (with the ``.files`` attribute), and for the arrays
    themselves.

    Examples
    --------
    >>> from tempfile import TemporaryFile
    >>> outfile = TemporaryFile()
    >>> x = np.arange(10)
    >>> y = np.sin(x)

    Using `savez` with *args, the arrays are saved with default names.

    >>> np.savez(outfile, x, y)
    >>> outfile.seek(0) # Only needed here to simulate closing & reopening file
    >>> npzfile = np.load(outfile)
    >>> npzfile.files
    ['arr_1', 'arr_0']
    >>> npzfile['arr_0']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

    Using `savez` with **kwds, the arrays are saved with the keyword names.

    >>> outfile = TemporaryFile()
    >>> np.savez(outfile, x=x, y=y)
    >>> outfile.seek(0)
    >>> npzfile = np.load(outfile)
    >>> npzfile.files
    ['y', 'x']
    >>> npzfile['x']
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

    See Also
    --------
    numpy.savez_compressed : Save several arrays into a compressed .npz file format
    
    """
    raise NotImplementedError


def savez_compressed(file):
    """Save several arrays into a single file in compressed ``.npz`` format.

    If keyword arguments are given, then filenames are taken from the keywords.
    If arguments are passed in with no keywords, then stored file names are
    arr_0, arr_1, etc.

    Parameters
    ----------
    file : str
        File name of .npz file.
    args : Arguments
        Function arguments.
    kwds : Keyword arguments
        Keywords.

    See Also
    --------
    numpy.savez : Save several arrays into an uncompressed .npz file format
    
    """
    raise NotImplementedError


def sctype2char(sctype):
    """Return the string representation of a scalar dtype.

    Parameters
    ----------
    sctype : scalar dtype or object
        If a scalar dtype, the corresponding string character is
        returned. If an object, `sctype2char` tries to infer its scalar type
        and then return the corresponding string character.

    Returns
    -------
    typechar : str
        The string character corresponding to the scalar type.

    Raises
    ------
    ValueError
        If `sctype` is an object for which the type can not be inferred.

    See Also
    --------
    obj2sctype, issctype, issubsctype, mintypecode

    Examples
    --------
    >>> for sctype in [np.int32, np.float, np.complex, np.string_, np.ndarray]:
    ...     print np.sctype2char(sctype)
    l
    d
    D
    S
    O

    >>> x = np.array([1., 2-1.j])
    >>> np.sctype2char(x)
    'D'
    >>> np.sctype2char(list)
    'O'
    
    """
    raise NotImplementedError


def searchsorted(a, v, side='left'):
    """Find indices where elements should be inserted to maintain order.

    Find the indices into a sorted array `a` such that, if the corresponding
    elements in `v` were inserted before the indices, the order of `a` would
    be preserved.

    Parameters
    ----------
    a : 1-D array_like
        Input array, sorted in ascending order.
    v : array_like
        Values to insert into `a`.
    side : {'left', 'right'}, optional
        If 'left', the index of the first suitable location found is given.  If
        'right', return the last such index.  If there is no suitable
        index, return either 0 or N (where N is the length of `a`).

    Returns
    -------
    indices : array of ints
        Array of insertion points with the same shape as `v`.

    See Also
    --------
    sort : Return a sorted copy of an array.
    histogram : Produce histogram from 1-D data.

    Notes
    -----
    Binary search is used to find the required insertion points.

    As of Numpy 1.4.0 `searchsorted` works with real/complex arrays containing
    `nan` values. The enhanced sort order is documented in `sort`.

    Examples
    --------
    >>> np.searchsorted([1,2,3,4,5], 3)
    2
    >>> np.searchsorted([1,2,3,4,5], 3, side='right')
    3
    >>> np.searchsorted([1,2,3,4,5], [-10, 10, 2, 3])
    array([0, 5, 1, 2])
    
    """
    raise NotImplementedError


def select(condlist, choicelist, default=0):
    """Return an array drawn from elements in choicelist, depending on conditions.

    Parameters
    ----------
    condlist : list of bool ndarrays
        The list of conditions which determine from which array in `choicelist`
        the output elements are taken. When multiple conditions are satisfied,
        the first one encountered in `condlist` is used.
    choicelist : list of ndarrays
        The list of arrays from which the output elements are taken. It has
        to be of the same length as `condlist`.
    default : scalar, optional
        The element inserted in `output` when all conditions evaluate to False.

    Returns
    -------
    output : ndarray
        The output at position m is the m-th element of the array in
        `choicelist` where the m-th element of the corresponding array in
        `condlist` is True.

    See Also
    --------
    where : Return elements from one of two arrays depending on condition.
    take, choose, compress, diag, diagonal

    Examples
    --------
    >>> x = np.arange(10)
    >>> condlist = [x<3, x>5]
    >>> choicelist = [x, x**2]
    >>> np.select(condlist, choicelist)
    array([ 0,  1,  2,  0,  0,  0, 36, 49, 64, 81])
    
    """
    raise NotImplementedError


def set_numeric_ops():
    """set_numeric_ops(op1=func1, op2=func2, ...)

    Set numerical operators for array objects.

    Parameters
    ----------
    op1, op2, ... : callable
        Each ``op = func`` pair describes an operator to be replaced.
        For example, ``add = lambda x, y: np.add(x, y) % 5`` would replace
        addition by modulus 5 addition.

    Returns
    -------
    saved_ops : list of callables
        A list of all operators, stored before making replacements.

    Notes
    -----
    .. WARNING::
       Use with care!  Incorrect usage may lead to memory errors.

    A function replacing an operator cannot make use of that operator.
    For example, when replacing add, you may not use ``+``.  Instead,
    directly call ufuncs.

    Examples
    --------
    >>> def add_mod5(x, y):
    ...     return np.add(x, y) % 5
    ...
    >>> old_funcs = np.set_numeric_ops(add=add_mod5)

    >>> x = np.arange(12).reshape((3, 4))
    >>> x + x
    array([[0, 2, 4, 1],
           [3, 0, 2, 4],
           [1, 3, 0, 2]])

    >>> ignore = np.set_numeric_ops(**old_funcs) # restore operators
    
    """
    # BUILTIN
    raise NotImplementedError


def set_printoptions(precision=None, threshold=None, edgeitems=None, linewidth=None, suppress=None, nanstr=None, infstr=None):
    """Set printing options.

    These options determine the way floating point numbers, arrays and
    other NumPy objects are displayed.

    Parameters
    ----------
    precision : int, optional
        Number of digits of precision for floating point output (default 8).
    threshold : int, optional
        Total number of array elements which trigger summarization
        rather than full repr (default 1000).
    edgeitems : int, optional
        Number of array items in summary at beginning and end of
        each dimension (default 3).
    linewidth : int, optional
        The number of characters per line for the purpose of inserting
        line breaks (default 75).
    suppress : bool, optional
        Whether or not suppress printing of small floating point values
        using scientific notation (default False).
    nanstr : str, optional
        String representation of floating point not-a-number (default nan).
    infstr : str, optional
        String representation of floating point infinity (default inf).

    See Also
    --------
    get_printoptions, set_string_function

    Examples
    --------
    Floating point precision can be set:

    >>> np.set_printoptions(precision=4)
    >>> print np.array([1.123456789])
    [ 1.1235]

    Long arrays can be summarised:

    >>> np.set_printoptions(threshold=5)
    >>> print np.arange(10)
    [0 1 2 ..., 7 8 9]

    Small results can be suppressed:

    >>> eps = np.finfo(float).eps
    >>> x = np.arange(4.)
    >>> x**2 - (x + eps)**2
    array([ -4.9304e-32,  -4.4409e-16,   0.0000e+00,   0.0000e+00])
    >>> np.set_printoptions(suppress=True)
    >>> x**2 - (x + eps)**2
    array([-0., -0.,  0.,  0.])

    To put back the default options, you can use:

    >>> np.set_printoptions(edgeitems=3,infstr='Inf',
    ... linewidth=75, nanstr='NaN', precision=8,
    ... suppress=False, threshold=1000)
    
    """
    raise NotImplementedError


def set_string_function(f, repr=True):
    """Set a Python function to be used when pretty printing arrays.

    Parameters
    ----------
    f : function or None
        Function to be used to pretty print arrays. The function should expect
        a single array argument and return a string of the representation of
        the array. If None, the function is reset to the default NumPy function
        to print arrays.
    repr : bool, optional
        If True (default), the function for pretty printing (``__repr__``)
        is set, if False the function that returns the default string
        representation (``__str__``) is set.

    See Also
    --------
    set_printoptions, get_printoptions

    Examples
    --------
    >>> def pprint(arr):
    ...     return 'HA! - What are you going to do now?'
    ...
    >>> np.set_string_function(pprint)
    >>> a = np.arange(10)
    >>> a
    HA! - What are you going to do now?
    >>> print a
    [0 1 2 3 4 5 6 7 8 9]

    We can reset the function to the default:

    >>> np.set_string_function(None)
    >>> a
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

    `repr` affects either pretty printing or normal string representation.
    Note that ``__repr__`` is still affected by setting ``__str__``
    because the width of each array element in the returned string becomes
    equal to the length of the result of ``__str__()``.

    >>> x = np.arange(4)
    >>> np.set_string_function(lambda x:'random', repr=False)
    >>> x.__str__()
    'random'
    >>> x.__repr__()
    'array([     0,      1,      2,      3])'
    
    """
    raise NotImplementedError


def setbufsize(size):
    """Set the size of the buffer used in ufuncs.

    Parameters
    ----------
    size : int
        Size of buffer.
    
    """
    raise NotImplementedError


def setdiff1d(ar1, ar2, assume_unique=False):
    """Find the set difference of two arrays.

    Return the sorted, unique values in `ar1` that are not in `ar2`.

    Parameters
    ----------
    ar1 : array_like
        Input array.
    ar2 : array_like
        Input comparison array.
    assume_unique : bool
        If True, the input arrays are both assumed to be unique, which
        can speed up the calculation.  Default is False.

    Returns
    -------
    difference : ndarray
        Sorted 1D array of values in `ar1` that are not in `ar2`.

    See Also
    --------
    numpy.lib.arraysetops : Module with a number of other functions for
                            performing set operations on arrays.

    Examples
    --------
    >>> a = np.array([1, 2, 3, 2, 4, 1])
    >>> b = np.array([3, 4, 5, 6])
    >>> np.setdiff1d(a, b)
    array([1, 2])
    
    """
    raise NotImplementedError


def seterr(all=None, divide=None, over=None, under=None, invalid=None):
    """Set how floating-point errors are handled.

    Note that operations on integer scalar types (such as `int16`) are
    handled like floating point, and are affected by these settings.

    Parameters
    ----------
    all : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
        Set treatment for all types of floating-point errors at once:

        - ignore: Take no action when the exception occurs.
        - warn: Print a `RuntimeWarning` (via the Python `warnings` module).
        - raise: Raise a `FloatingPointError`.
        - call: Call a function specified using the `seterrcall` function.
        - print: Print a warning directly to ``stdout``.
        - log: Record error in a Log object specified by `seterrcall`.

        The default is not to change the current behavior.
    divide : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
        Treatment for division by zero.
    over : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
        Treatment for floating-point overflow.
    under : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
        Treatment for floating-point underflow.
    invalid : {'ignore', 'warn', 'raise', 'call', 'print', 'log'}, optional
        Treatment for invalid floating-point operation.

    Returns
    -------
    old_settings : dict
        Dictionary containing the old settings.

    See also
    --------
    seterrcall : Set a callback function for the 'call' mode.
    geterr, geterrcall

    Notes
    -----
    The floating-point exceptions are defined in the IEEE 754 standard [1]:

    - Division by zero: infinite result obtained from finite numbers.
    - Overflow: result too large to be expressed.
    - Underflow: result so close to zero that some precision
      was lost.
    - Invalid operation: result is not an expressible number, typically
      indicates that a NaN was produced.

    .. [1] http://en.wikipedia.org/wiki/IEEE_754

    Examples
    --------
    >>> old_settings = np.seterr(all='ignore')  #seterr to known value
    >>> np.seterr(over='raise')
    {'over': 'ignore', 'divide': 'ignore', 'invalid': 'ignore',
     'under': 'ignore'}
    >>> np.seterr(all='ignore')  # reset to default
    {'over': 'raise', 'divide': 'ignore', 'invalid': 'ignore', 'under': 'ignore'}

    >>> np.int16(32000) * np.int16(3)
    30464
    >>> old_settings = np.seterr(all='warn', over='raise')
    >>> np.int16(32000) * np.int16(3)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    FloatingPointError: overflow encountered in short_scalars

    >>> old_settings = np.seterr(all='print')
    >>> np.geterr()
    {'over': 'print', 'divide': 'print', 'invalid': 'print', 'under': 'print'}
    >>> np.int16(32000) * np.int16(3)
    Warning: overflow encountered in short_scalars
    30464
    
    """
    raise NotImplementedError


def seterrcall(func):
    """Set the floating-point error callback function or log object.

    There are two ways to capture floating-point error messages.  The first
    is to set the error-handler to 'call', using `seterr`.  Then, set
    the function to call using this function.

    The second is to set the error-handler to 'log', using `seterr`.
    Floating-point errors then trigger a call to the 'write' method of
    the provided object.

    Parameters
    ----------
    func : callable f(err, flag) or object with write method
        Function to call upon floating-point errors ('call'-mode) or
        object whose 'write' method is used to log such message ('log'-mode).

        The call function takes two arguments. The first is the
        type of error (one of "divide", "over", "under", or "invalid"),
        and the second is the status flag.  The flag is a byte, whose
        least-significant bits indicate the status::

          [0 0 0 0 invalid over under invalid]

        In other words, ``flags = divide + 2*over + 4*under + 8*invalid``.

        If an object is provided, its write method should take one argument,
        a string.

    Returns
    -------
    h : callable, log instance or None
        The old error handler.

    See Also
    --------
    seterr, geterr, geterrcall

    Examples
    --------
    Callback upon error:

    >>> def err_handler(type, flag):
    ...     print "Floating point error (%s), with flag %s" % (type, flag)
    ...

    >>> saved_handler = np.seterrcall(err_handler)
    >>> save_err = np.seterr(all='call')

    >>> np.array([1, 2, 3]) / 0.0
    Floating point error (divide by zero), with flag 1
    array([ Inf,  Inf,  Inf])

    >>> np.seterrcall(saved_handler)
    <function err_handler at 0x...>
    >>> np.seterr(**save_err)
    {'over': 'call', 'divide': 'call', 'invalid': 'call', 'under': 'call'}

    Log error message:

    >>> class Log(object):
    ...     def write(self, msg):
    ...         print "LOG: %s" % msg
    ...

    >>> log = Log()
    >>> saved_handler = np.seterrcall(log)
    >>> save_err = np.seterr(all='log')

    >>> np.array([1, 2, 3]) / 0.0
    LOG: Warning: divide by zero encountered in divide
    <BLANKLINE>
    array([ Inf,  Inf,  Inf])

    >>> np.seterrcall(saved_handler)
    <__main__.Log object at 0x...>
    >>> np.seterr(**save_err)
    {'over': 'log', 'divide': 'log', 'invalid': 'log', 'under': 'log'}
    
    """
    raise NotImplementedError


def seterrobj(errobj):
    """seterrobj(errobj)

    Set the object that defines floating-point error handling.

    The error object contains all information that defines the error handling
    behavior in Numpy. `seterrobj` is used internally by the other
    functions that set error handling behavior (`seterr`, `seterrcall`).

    Parameters
    ----------
    errobj : list
        The error object, a list containing three elements:
        [internal numpy buffer size, error mask, error callback function].

        The error mask is a single integer that holds the treatment information
        on all four floating point errors. The information for each error type
        is contained in three bits of the integer. If we print it in base 8, we
        can see what treatment is set for "invalid", "under", "over", and
        "divide" (in that order). The printed string can be interpreted with

        * 0 : 'ignore'
        * 1 : 'warn'
        * 2 : 'raise'
        * 3 : 'call'
        * 4 : 'print'
        * 5 : 'log'

    See Also
    --------
    geterrobj, seterr, geterr, seterrcall, geterrcall
    getbufsize, setbufsize

    Notes
    -----
    For complete documentation of the types of floating-point exceptions and
    treatment options, see `seterr`.

    Examples
    --------
    >>> old_errobj = np.geterrobj()  # first get the defaults
    >>> old_errobj
    [10000, 0, None]

    >>> def err_handler(type, flag):
    ...     print "Floating point error (%s), with flag %s" % (type, flag)
    ...
    >>> new_errobj = [20000, 12, err_handler]
    >>> np.seterrobj(new_errobj)
    >>> np.base_repr(12, 8)  # int for divide=4 ('print') and over=1 ('warn')
    '14'
    >>> np.geterr()
    {'over': 'warn', 'divide': 'print', 'invalid': 'ignore', 'under': 'ignore'}
    >>> np.geterrcall() is err_handler
    True
    
    """
    # BUILTIN
    raise NotImplementedError


def setxor1d(ar1, ar2, assume_unique=False):
    """Find the set exclusive-or of two arrays.

    Return the sorted, unique values that are in only one (not both) of the
    input arrays.

    Parameters
    ----------
    ar1, ar2 : array_like
        Input arrays.
    assume_unique : bool
        If True, the input arrays are both assumed to be unique, which
        can speed up the calculation.  Default is False.

    Returns
    -------
    xor : ndarray
        Sorted 1D array of unique values that are in only one of the input
        arrays.

    Examples
    --------
    >>> a = np.array([1, 2, 3, 2, 4])
    >>> b = np.array([2, 3, 5, 7, 5])
    >>> np.setxor1d(a,b)
    array([1, 4, 5, 7])
    
    """
    raise NotImplementedError


def show_config():
    """
    
    """
    raise NotImplementedError


def sinc(x):
    """Return the sinc function.

    The sinc function is :math:`\sin(\pi x)/(\pi x)`.

    Parameters
    ----------
    x : ndarray
        Array (possibly multi-dimensional) of values for which to to
        calculate ``sinc(x)``.

    Returns
    -------
    out : ndarray
        ``sinc(x)``, which has the same shape as the input.

    Notes
    -----
    ``sinc(0)`` is the limit value 1.

    The name sinc is short for "sine cardinal" or "sinus cardinalis".

    The sinc function is used in various signal processing applications,
    including in anti-aliasing, in the construction of a
    Lanczos resampling filter, and in interpolation.

    For bandlimited interpolation of discrete-time signals, the ideal
    interpolation kernel is proportional to the sinc function.

    References
    ----------
    .. [1] Weisstein, Eric W. "Sinc Function." From MathWorld--A Wolfram Web
           Resource. http://mathworld.wolfram.com/SincFunction.html
    .. [2] Wikipedia, "Sinc function",
           http://en.wikipedia.org/wiki/Sinc_function

    Examples
    --------
    >>> x = np.arange(-20., 21.)/5.
    >>> np.sinc(x)
    array([ -3.89804309e-17,  -4.92362781e-02,  -8.40918587e-02,
            -8.90384387e-02,  -5.84680802e-02,   3.89804309e-17,
             6.68206631e-02,   1.16434881e-01,   1.26137788e-01,
             8.50444803e-02,  -3.89804309e-17,  -1.03943254e-01,
            -1.89206682e-01,  -2.16236208e-01,  -1.55914881e-01,
             3.89804309e-17,   2.33872321e-01,   5.04551152e-01,
             7.56826729e-01,   9.35489284e-01,   1.00000000e+00,
             9.35489284e-01,   7.56826729e-01,   5.04551152e-01,
             2.33872321e-01,   3.89804309e-17,  -1.55914881e-01,
            -2.16236208e-01,  -1.89206682e-01,  -1.03943254e-01,
            -3.89804309e-17,   8.50444803e-02,   1.26137788e-01,
             1.16434881e-01,   6.68206631e-02,   3.89804309e-17,
            -5.84680802e-02,  -8.90384387e-02,  -8.40918587e-02,
            -4.92362781e-02,  -3.89804309e-17])

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x, np.sinc(x))
    [<matplotlib.lines.Line2D object at 0x...>]
    >>> plt.title("Sinc Function")
    <matplotlib.text.Text object at 0x...>
    >>> plt.ylabel("Amplitude")
    <matplotlib.text.Text object at 0x...>
    >>> plt.xlabel("X")
    <matplotlib.text.Text object at 0x...>
    >>> plt.show()

    It works in 2-D as well:

    >>> x = np.arange(-200., 201.)/50.
    >>> xx = np.outer(x, x)
    >>> plt.imshow(np.sinc(xx))
    <matplotlib.image.AxesImage object at 0x...>
    
    """
    raise NotImplementedError


def size(a, axis=None):
    """Return the number of elements along a given axis.

    Parameters
    ----------
    a : array_like
        Input data.
    axis : int, optional
        Axis along which the elements are counted.  By default, give
        the total number of elements.

    Returns
    -------
    element_count : int
        Number of elements along the specified axis.

    See Also
    --------
    shape : dimensions of array
    ndarray.shape : dimensions of array
    ndarray.size : number of elements in array

    Examples
    --------
    >>> a = np.array([[1,2,3],[4,5,6]])
    >>> np.size(a)
    6
    >>> np.size(a,1)
    3
    >>> np.size(a,0)
    2
    
    """
    raise NotImplementedError


def sometrue(a, axis=None, out=None):
    """Check whether some values are true.

    Refer to `any` for full documentation.

    See Also
    --------
    any : equivalent function
    
    """
    raise NotImplementedError


def sort(a, axis=-1, kind='quicksort', order=None):
    """Return a sorted copy of an array.

    Parameters
    ----------
    a : array_like
        Array to be sorted.
    axis : int or None, optional
        Axis along which to sort. If None, the array is flattened before
        sorting. The default is -1, which sorts along the last axis.
    kind : {'quicksort', 'mergesort', 'heapsort'}, optional
        Sorting algorithm. Default is 'quicksort'.
    order : list, optional
        When `a` is a structured array, this argument specifies which fields
        to compare first, second, and so on.  This list does not need to
        include all of the fields.

    Returns
    -------
    sorted_array : ndarray
        Array of the same type and shape as `a`.

    See Also
    --------
    ndarray.sort : Method to sort an array in-place.
    argsort : Indirect sort.
    lexsort : Indirect stable sort on multiple keys.
    searchsorted : Find elements in a sorted array.

    Notes
    -----
    The various sorting algorithms are characterized by their average speed,
    worst case performance, work space size, and whether they are stable. A
    stable sort keeps items with the same key in the same relative
    order. The three available algorithms have the following
    properties:

    =========== ======= ============= ============ =======
       kind      speed   worst case    work space  stable
    =========== ======= ============= ============ =======
    'quicksort'    1     O(n^2)            0          no
    'mergesort'    2     O(n*log(n))      ~n/2        yes
    'heapsort'     3     O(n*log(n))       0          no
    =========== ======= ============= ============ =======

    All the sort algorithms make temporary copies of the data when
    sorting along any but the last axis.  Consequently, sorting along
    the last axis is faster and uses less space than sorting along
    any other axis.

    The sort order for complex numbers is lexicographic. If both the real
    and imaginary parts are non-nan then the order is determined by the
    real parts except when they are equal, in which case the order is
    determined by the imaginary parts.

    Previous to numpy 1.4.0 sorting real and complex arrays containing nan
    values led to undefined behaviour. In numpy versions >= 1.4.0 nan
    values are sorted to the end. The extended sort order is:

      * Real: [R, nan]
      * Complex: [R + Rj, R + nanj, nan + Rj, nan + nanj]

    where R is a non-nan real value. Complex values with the same nan
    placements are sorted according to the non-nan part if it exists.
    Non-nan values are sorted as before.

    Examples
    --------
    >>> a = np.array([[1,4],[3,1]])
    >>> np.sort(a)                # sort along the last axis
    array([[1, 4],
           [1, 3]])
    >>> np.sort(a, axis=None)     # sort the flattened array
    array([1, 1, 3, 4])
    >>> np.sort(a, axis=0)        # sort along the first axis
    array([[1, 1],
           [3, 4]])

    Use the `order` keyword to specify a field to use when sorting a
    structured array:

    >>> dtype = [('name', 'S10'), ('height', float), ('age', int)]
    >>> values = [('Arthur', 1.8, 41), ('Lancelot', 1.9, 38),
    ...           ('Galahad', 1.7, 38)]
    >>> a = np.array(values, dtype=dtype)       # create a structured array
    >>> np.sort(a, order='height')                        # doctest: +SKIP
    array([('Galahad', 1.7, 38), ('Arthur', 1.8, 41),
           ('Lancelot', 1.8999999999999999, 38)],
          dtype=[('name', '|S10'), ('height', '<f8'), ('age', '<i4')])

    Sort by age, then height if ages are equal:

    >>> np.sort(a, order=['age', 'height'])               # doctest: +SKIP
    array([('Galahad', 1.7, 38), ('Lancelot', 1.8999999999999999, 38),
           ('Arthur', 1.8, 41)],
          dtype=[('name', '|S10'), ('height', '<f8'), ('age', '<i4')])
    
    """
    raise NotImplementedError


def sort_complex(a):
    """Sort a complex array using the real part first, then the imaginary part.

    Parameters
    ----------
    a : array_like
        Input array

    Returns
    -------
    out : complex ndarray
        Always returns a sorted complex array.

    Examples
    --------
    >>> np.sort_complex([5, 3, 6, 2, 1])
    array([ 1.+0.j,  2.+0.j,  3.+0.j,  5.+0.j,  6.+0.j])

    >>> np.sort_complex([1 + 2j, 2 - 1j, 3 - 2j, 3 - 3j, 3 + 5j])
    array([ 1.+2.j,  2.-1.j,  3.-3.j,  3.-2.j,  3.+5.j])
    
    """
    raise NotImplementedError


def source(object, output=sys.stdout):
    """Print or write to a file the source code for a Numpy object.

    The source code is only returned for objects written in Python. Many
    functions and classes are defined in C and will therefore not return
    useful information.

    Parameters
    ----------
    object : numpy object
        Input object. This can be any object (function, class, module, ...).
    output : file object, optional
        If `output` not supplied then source code is printed to screen
        (sys.stdout).  File object must be created with either write 'w' or
        append 'a' modes.

    See Also
    --------
    lookfor, info

    Examples
    --------
    >>> np.source(np.interp)                        #doctest: +SKIP
    In file: /usr/lib/python2.6/dist-packages/numpy/lib/function_base.py
    def interp(x, xp, fp, left=None, right=None):
        '''.... (full docstring printed)'''
        if isinstance(x, (float, int, number)):
            return compiled_interp([x], xp, fp, left, right).item()
        else:
            return compiled_interp(x, xp, fp, left, right)

    The source code is only returned for objects written in Python.

    >>> np.source(np.array)                         #doctest: +SKIP
    Not available for this object.
    
    """
    raise NotImplementedError


def split(ary, indices_or_sections, axis=0):
    """Split an array into multiple sub-arrays of equal size.

    Parameters
    ----------
    ary : ndarray
        Array to be divided into sub-arrays.
    indices_or_sections : int or 1-D array
        If `indices_or_sections` is an integer, N, the array will be divided
        into N equal arrays along `axis`.  If such a split is not possible,
        an error is raised.

        If `indices_or_sections` is a 1-D array of sorted integers, the entries
        indicate where along `axis` the array is split.  For example,
        ``[2, 3]`` would, for ``axis=0``, result in

          - ary[:2]
          - ary[2:3]
          - ary[3:]

        If an index exceeds the dimension of the array along `axis`,
        an empty sub-array is returned correspondingly.
    axis : int, optional
        The axis along which to split, default is 0.

    Returns
    -------
    sub-arrays : list of ndarrays
        A list of sub-arrays.

    Raises
    ------
    ValueError
        If `indices_or_sections` is given as an integer, but
        a split does not result in equal division.

    See Also
    --------
    array_split : Split an array into multiple sub-arrays of equal or
                  near-equal size.  Does not raise an exception if
                  an equal division cannot be made.
    hsplit : Split array into multiple sub-arrays horizontally (column-wise).
    vsplit : Split array into multiple sub-arrays vertically (row wise).
    dsplit : Split array into multiple sub-arrays along the 3rd axis (depth).
    concatenate : Join arrays together.
    hstack : Stack arrays in sequence horizontally (column wise).
    vstack : Stack arrays in sequence vertically (row wise).
    dstack : Stack arrays in sequence depth wise (along third dimension).

    Examples
    --------
    >>> x = np.arange(9.0)
    >>> np.split(x, 3)
    [array([ 0.,  1.,  2.]), array([ 3.,  4.,  5.]), array([ 6.,  7.,  8.])]

    >>> x = np.arange(8.0)
    >>> np.split(x, [3, 5, 6, 10])
    [array([ 0.,  1.,  2.]),
     array([ 3.,  4.]),
     array([ 5.]),
     array([ 6.,  7.]),
     array([], dtype=float64)]
    
    """
    raise NotImplementedError


def squeeze(a):
    """Remove single-dimensional entries from the shape of an array.

    Parameters
    ----------
    a : array_like
        Input data.

    Returns
    -------
    squeezed : ndarray
        The input array, but with with all dimensions of length 1
        removed.  Whenever possible, a view on `a` is returned.

    Examples
    --------
    >>> x = np.array([[[0], [1], [2]]])
    >>> x.shape
    (1, 3, 1)
    >>> np.squeeze(x).shape
    (3,)
    
    """
    raise NotImplementedError


def std(a, axis=None, dtype=None, out=None, ddof=0):
    """Compute the standard deviation along the specified axis.

    Returns the standard deviation, a measure of the spread of a distribution,
    of the array elements. The standard deviation is computed for the
    flattened array by default, otherwise over the specified axis.

    Parameters
    ----------
    a : array_like
        Calculate the standard deviation of these values.
    axis : int, optional
        Axis along which the standard deviation is computed. The default is
        to compute the standard deviation of the flattened array.
    dtype : dtype, optional
        Type to use in computing the standard deviation. For arrays of
        integer type the default is float64, for arrays of float types it is
        the same as the array type.
    out : ndarray, optional
        Alternative output array in which to place the result. It must have
        the same shape as the expected output but the type (of the calculated
        values) will be cast if necessary.
    ddof : int, optional
        Means Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.
        By default `ddof` is zero.

    Returns
    -------
    standard_deviation : ndarray, see dtype parameter above.
        If `out` is None, return a new array containing the standard deviation,
        otherwise return a reference to the output array.

    See Also
    --------
    var, mean
    numpy.doc.ufuncs : Section "Output arguments"

    Notes
    -----
    The standard deviation is the square root of the average of the squared
    deviations from the mean, i.e., ``std = sqrt(mean(abs(x - x.mean())**2))``.

    The average squared deviation is normally calculated as ``x.sum() / N``, where
    ``N = len(x)``.  If, however, `ddof` is specified, the divisor ``N - ddof``
    is used instead. In standard statistical practice, ``ddof=1`` provides an
    unbiased estimator of the variance of the infinite population. ``ddof=0``
    provides a maximum likelihood estimate of the variance for normally
    distributed variables. The standard deviation computed in this function
    is the square root of the estimated variance, so even with ``ddof=1``, it
    will not be an unbiased estimate of the standard deviation per se.

    Note that, for complex numbers, `std` takes the absolute
    value before squaring, so that the result is always real and nonnegative.

    For floating-point input, the *std* is computed using the same
    precision the input has. Depending on the input data, this can cause
    the results to be inaccurate, especially for float32 (see example below).
    Specifying a higher-accuracy accumulator using the `dtype` keyword can
    alleviate this issue.

    Examples
    --------
    >>> a = np.array([[1, 2], [3, 4]])
    >>> np.std(a)
    1.1180339887498949
    >>> np.std(a, axis=0)
    array([ 1.,  1.])
    >>> np.std(a, axis=1)
    array([ 0.5,  0.5])

    In single precision, std() can be inaccurate:

    >>> a = np.zeros((2,512*512), dtype=np.float32)
    >>> a[0,:] = 1.0
    >>> a[1,:] = 0.1
    >>> np.std(a)
    0.45172946707416706

    Computing the standard deviation in float64 is more accurate:

    >>> np.std(a, dtype=np.float64)
    0.44999999925552653
    
    """
    raise NotImplementedError


def sum(a, axis=None, dtype=None, out=None):
    """Sum of array elements over a given axis.

    Parameters
    ----------
    a : array_like
        Elements to sum.
    axis : integer, optional
        Axis over which the sum is taken. By default `axis` is None,
        and all elements are summed.
    dtype : dtype, optional
        The type of the returned array and of the accumulator in which
        the elements are summed.  By default, the dtype of `a` is used.
        An exception is when `a` has an integer type with less precision
        than the default platform integer.  In that case, the default
        platform integer is used instead.
    out : ndarray, optional
        Array into which the output is placed.  By default, a new array is
        created.  If `out` is given, it must be of the appropriate shape
        (the shape of `a` with `axis` removed, i.e.,
        ``numpy.delete(a.shape, axis)``).  Its type is preserved. See
        `doc.ufuncs` (Section "Output arguments") for more details.

    Returns
    -------
    sum_along_axis : ndarray
        An array with the same shape as `a`, with the specified
        axis removed.   If `a` is a 0-d array, or if `axis` is None, a scalar
        is returned.  If an output array is specified, a reference to
        `out` is returned.

    See Also
    --------
    ndarray.sum : Equivalent method.

    cumsum : Cumulative sum of array elements.

    trapz : Integration of array values using the composite trapezoidal rule.

    mean, average

    Notes
    -----
    Arithmetic is modular when using integer types, and no error is
    raised on overflow.

    Examples
    --------
    >>> np.sum([0.5, 1.5])
    2.0
    >>> np.sum([0.5, 0.7, 0.2, 1.5], dtype=np.int32)
    1
    >>> np.sum([[0, 1], [0, 5]])
    6
    >>> np.sum([[0, 1], [0, 5]], axis=0)
    array([0, 6])
    >>> np.sum([[0, 1], [0, 5]], axis=1)
    array([1, 5])

    If the accumulator is too small, overflow occurs:

    >>> np.ones(128, dtype=np.int8).sum(dtype=np.int8)
    -128
    
    """
    raise NotImplementedError


def swapaxes(a, axis1, axis2):
    """Interchange two axes of an array.

    Parameters
    ----------
    a : array_like
        Input array.
    axis1 : int
        First axis.
    axis2 : int
        Second axis.

    Returns
    -------
    a_swapped : ndarray
        If `a` is an ndarray, then a view of `a` is returned; otherwise
        a new array is created.

    Examples
    --------
    >>> x = np.array([[1,2,3]])
    >>> np.swapaxes(x,0,1)
    array([[1],
           [2],
           [3]])

    >>> x = np.array([[[0,1],[2,3]],[[4,5],[6,7]]])
    >>> x
    array([[[0, 1],
            [2, 3]],
           [[4, 5],
            [6, 7]]])

    >>> np.swapaxes(x,0,2)
    array([[[0, 4],
            [2, 6]],
           [[1, 5],
            [3, 7]]])
    
    """
    raise NotImplementedError


def take(a, indices, axis=None, out=None, mode='raise'):
    """Take elements from an array along an axis.

    This function does the same thing as "fancy" indexing (indexing arrays
    using arrays); however, it can be easier to use if you need elements
    along a given axis.

    Parameters
    ----------
    a : array_like
        The source array.
    indices : array_like
        The indices of the values to extract.
    axis : int, optional
        The axis over which to select values. By default, the flattened
        input array is used.
    out : ndarray, optional
        If provided, the result will be placed in this array. It should
        be of the appropriate shape and dtype.
    mode : {'raise', 'wrap', 'clip'}, optional
        Specifies how out-of-bounds indices will behave.

        * 'raise' -- raise an error (default)
        * 'wrap' -- wrap around
        * 'clip' -- clip to the range

        'clip' mode means that all indices that are too large are replaced
        by the index that addresses the last element along that axis. Note
        that this disables indexing with negative numbers.

    Returns
    -------
    subarray : ndarray
        The returned array has the same type as `a`.

    See Also
    --------
    ndarray.take : equivalent method

    Examples
    --------
    >>> a = [4, 3, 5, 7, 6, 8]
    >>> indices = [0, 1, 4]
    >>> np.take(a, indices)
    array([4, 3, 6])

    In this example if `a` is an ndarray, "fancy" indexing can be used.

    >>> a = np.array(a)
    >>> a[indices]
    array([4, 3, 6])
    
    """
    raise NotImplementedError


def tensordot(a, b, axes=2):
    """Compute tensor dot product along specified axes for arrays >= 1-D.

    Given two tensors (arrays of dimension greater than or equal to one),
    ``a`` and ``b``, and an array_like object containing two array_like
    objects, ``(a_axes, b_axes)``, sum the products of ``a``'s and ``b``'s
    elements (components) over the axes specified by ``a_axes`` and
    ``b_axes``. The third argument can be a single non-negative
    integer_like scalar, ``N``; if it is such, then the last ``N``
    dimensions of ``a`` and the first ``N`` dimensions of ``b`` are summed
    over.

    Parameters
    ----------
    a, b : array_like, len(shape) >= 1
        Tensors to "dot".

    axes : variable type

    * integer_like scalar
      Number of axes to sum over (applies to both arrays); or

    * array_like, shape = (2,), both elements array_like
      Axes to be summed over, first sequence applying to ``a``, second
      to ``b``.

    See Also
    --------
    dot, einsum

    Notes
    -----
    When there is more than one axis to sum over - and they are not the last
    (first) axes of ``a`` (``b``) - the argument ``axes`` should consist of
    two sequences of the same length, with the first axis to sum over given
    first in both sequences, the second axis second, and so forth.

    Examples
    --------
    A "traditional" example:

    >>> a = np.arange(60.).reshape(3,4,5)
    >>> b = np.arange(24.).reshape(4,3,2)
    >>> c = np.tensordot(a,b, axes=([1,0],[0,1]))
    >>> c.shape
    (5, 2)
    >>> c
    array([[ 4400.,  4730.],
           [ 4532.,  4874.],
           [ 4664.,  5018.],
           [ 4796.,  5162.],
           [ 4928.,  5306.]])
    >>> # A slower but equivalent way of computing the same...
    >>> d = np.zeros((5,2))
    >>> for i in range(5):
    ...   for j in range(2):
    ...     for k in range(3):
    ...       for n in range(4):
    ...         d[i,j] += a[k,n,i] * b[n,k,j]
    >>> c == d
    array([[ True,  True],
           [ True,  True],
           [ True,  True],
           [ True,  True],
           [ True,  True]], dtype=bool)

    An extended example taking advantage of the overloading of + and \*:

    >>> a = np.array(range(1, 9))
    >>> a.shape = (2, 2, 2)
    >>> A = np.array(('a', 'b', 'c', 'd'), dtype=object)
    >>> A.shape = (2, 2)
    >>> a; A
    array([[[1, 2],
            [3, 4]],
           [[5, 6],
            [7, 8]]])
    array([[a, b],
           [c, d]], dtype=object)

    >>> np.tensordot(a, A) # third argument default is 2
    array([abbcccdddd, aaaaabbbbbbcccccccdddddddd], dtype=object)

    >>> np.tensordot(a, A, 1)
    array([[[acc, bdd],
            [aaacccc, bbbdddd]],
           [[aaaaacccccc, bbbbbdddddd],
            [aaaaaaacccccccc, bbbbbbbdddddddd]]], dtype=object)

    >>> np.tensordot(a, A, 0) # "Left for reader" (result too long to incl.)
    array([[[[[a, b],
              [c, d]],
              ...

    >>> np.tensordot(a, A, (0, 1))
    array([[[abbbbb, cddddd],
            [aabbbbbb, ccdddddd]],
           [[aaabbbbbbb, cccddddddd],
            [aaaabbbbbbbb, ccccdddddddd]]], dtype=object)

    >>> np.tensordot(a, A, (2, 1))
    array([[[abb, cdd],
            [aaabbbb, cccdddd]],
           [[aaaaabbbbbb, cccccdddddd],
            [aaaaaaabbbbbbbb, cccccccdddddddd]]], dtype=object)

    >>> np.tensordot(a, A, ((0, 1), (0, 1)))
    array([abbbcccccddddddd, aabbbbccccccdddddddd], dtype=object)

    >>> np.tensordot(a, A, ((2, 1), (1, 0)))
    array([acccbbdddd, aaaaacccccccbbbbbbdddddddd], dtype=object)
    
    """
    raise NotImplementedError


def tile(A, reps):
    """Construct an array by repeating A the number of times given by reps.

    If `reps` has length ``d``, the result will have dimension of
    ``max(d, A.ndim)``.

    If ``A.ndim < d``, `A` is promoted to be d-dimensional by prepending new
    axes. So a shape (3,) array is promoted to (1, 3) for 2-D replication,
    or shape (1, 1, 3) for 3-D replication. If this is not the desired
    behavior, promote `A` to d-dimensions manually before calling this
    function.

    If ``A.ndim > d``, `reps` is promoted to `A`.ndim by pre-pending 1's to it.
    Thus for an `A` of shape (2, 3, 4, 5), a `reps` of (2, 2) is treated as
    (1, 1, 2, 2).

    Parameters
    ----------
    A : array_like
        The input array.
    reps : array_like
        The number of repetitions of `A` along each axis.

    Returns
    -------
    c : ndarray
        The tiled output array.

    See Also
    --------
    repeat : Repeat elements of an array.

    Examples
    --------
    >>> a = np.array([0, 1, 2])
    >>> np.tile(a, 2)
    array([0, 1, 2, 0, 1, 2])
    >>> np.tile(a, (2, 2))
    array([[0, 1, 2, 0, 1, 2],
           [0, 1, 2, 0, 1, 2]])
    >>> np.tile(a, (2, 1, 2))
    array([[[0, 1, 2, 0, 1, 2]],
           [[0, 1, 2, 0, 1, 2]]])

    >>> b = np.array([[1, 2], [3, 4]])
    >>> np.tile(b, 2)
    array([[1, 2, 1, 2],
           [3, 4, 3, 4]])
    >>> np.tile(b, (2, 1))
    array([[1, 2],
           [3, 4],
           [1, 2],
           [3, 4]])
    
    """
    raise NotImplementedError


def trace(a, offset=0, axis1=0, axis2=1, dtype=None, out=None):
    """Return the sum along diagonals of the array.

    If `a` is 2-D, the sum along its diagonal with the given offset
    is returned, i.e., the sum of elements ``a[i,i+offset]`` for all i.

    If `a` has more than two dimensions, then the axes specified by axis1 and
    axis2 are used to determine the 2-D sub-arrays whose traces are returned.
    The shape of the resulting array is the same as that of `a` with `axis1`
    and `axis2` removed.

    Parameters
    ----------
    a : array_like
        Input array, from which the diagonals are taken.
    offset : int, optional
        Offset of the diagonal from the main diagonal. Can be both positive
        and negative. Defaults to 0.
    axis1, axis2 : int, optional
        Axes to be used as the first and second axis of the 2-D sub-arrays
        from which the diagonals should be taken. Defaults are the first two
        axes of `a`.
    dtype : dtype, optional
        Determines the data-type of the returned array and of the accumulator
        where the elements are summed. If dtype has the value None and `a` is
        of integer type of precision less than the default integer
        precision, then the default integer precision is used. Otherwise,
        the precision is the same as that of `a`.
    out : ndarray, optional
        Array into which the output is placed. Its type is preserved and
        it must be of the right shape to hold the output.

    Returns
    -------
    sum_along_diagonals : ndarray
        If `a` is 2-D, the sum along the diagonal is returned.  If `a` has
        larger dimensions, then an array of sums along diagonals is returned.

    See Also
    --------
    diag, diagonal, diagflat

    Examples
    --------
    >>> np.trace(np.eye(3))
    3.0
    >>> a = np.arange(8).reshape((2,2,2))
    >>> np.trace(a)
    array([6, 8])

    >>> a = np.arange(24).reshape((2,2,2,3))
    >>> np.trace(a).shape
    (2, 3)
    
    """
    raise NotImplementedError


def transpose(a, axes=None):
    """Permute the dimensions of an array.

    Parameters
    ----------
    a : array_like
        Input array.
    axes : list of ints, optional
        By default, reverse the dimensions, otherwise permute the axes
        according to the values given.

    Returns
    -------
    p : ndarray
        `a` with its axes permuted.  A view is returned whenever
        possible.

    See Also
    --------
    rollaxis

    Examples
    --------
    >>> x = np.arange(4).reshape((2,2))
    >>> x
    array([[0, 1],
           [2, 3]])

    >>> np.transpose(x)
    array([[0, 2],
           [1, 3]])

    >>> x = np.ones((1, 2, 3))
    >>> np.transpose(x, (1, 0, 2)).shape
    (2, 1, 3)
    
    """
    raise NotImplementedError


def trapz(y, x=None, dx=1.0, axis=-1):
    """Integrate along the given axis using the composite trapezoidal rule.

    Integrate `y` (`x`) along given axis.

    Parameters
    ----------
    y : array_like
        Input array to integrate.
    x : array_like, optional
        If `x` is None, then spacing between all `y` elements is `dx`.
    dx : scalar, optional
        If `x` is None, spacing given by `dx` is assumed. Default is 1.
    axis : int, optional
        Specify the axis.

    Returns
    -------
    out : float
        Definite integral as approximated by trapezoidal rule.

    See Also
    --------
    sum, cumsum

    Notes
    -----
    Image [2]_ illustrates trapezoidal rule -- y-axis locations of points will
    be taken from `y` array, by default x-axis distances between points will be
    1.0, alternatively they can be provided with `x` array or with `dx` scalar.
    Return value will be equal to combined area under the red lines.


    References
    ----------
    .. [1] Wikipedia page: http://en.wikipedia.org/wiki/Trapezoidal_rule

    .. [2] Illustration image:
           http://en.wikipedia.org/wiki/File:Composite_trapezoidal_rule_illustration.png

    Examples
    --------
    >>> np.trapz([1,2,3])
    4.0
    >>> np.trapz([1,2,3], x=[4,6,8])
    8.0
    >>> np.trapz([1,2,3], dx=2)
    8.0
    >>> a = np.arange(6).reshape(2, 3)
    >>> a
    array([[0, 1, 2],
           [3, 4, 5]])
    >>> np.trapz(a, axis=0)
    array([ 1.5,  2.5,  3.5])
    >>> np.trapz(a, axis=1)
    array([ 2.,  8.])
    
    """
    raise NotImplementedError


def tri(N, M=None, k=0, dtype=float):
    """An array with ones at and below the given diagonal and zeros elsewhere.

    Parameters
    ----------
    N : int
        Number of rows in the array.
    M : int, optional
        Number of columns in the array.
        By default, `M` is taken equal to `N`.
    k : int, optional
        The sub-diagonal at and below which the array is filled.
        `k` = 0 is the main diagonal, while `k` < 0 is below it,
        and `k` > 0 is above.  The default is 0.
    dtype : dtype, optional
        Data type of the returned array.  The default is float.

    Returns
    -------
    T : ndarray of shape (N, M)
        Array with its lower triangle filled with ones and zero elsewhere;
        in other words ``T[i,j] == 1`` for ``i <= j + k``, 0 otherwise.

    Examples
    --------
    >>> np.tri(3, 5, 2, dtype=int)
    array([[1, 1, 1, 0, 0],
           [1, 1, 1, 1, 0],
           [1, 1, 1, 1, 1]])

    >>> np.tri(3, 5, -1)
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 1.,  0.,  0.,  0.,  0.],
           [ 1.,  1.,  0.,  0.,  0.]])
    
    """
    raise NotImplementedError


def tril(m, k=0):
    """Lower triangle of an array.

    Return a copy of an array with elements above the `k`-th diagonal zeroed.

    Parameters
    ----------
    m : array_like, shape (M, N)
        Input array.
    k : int, optional
        Diagonal above which to zero elements.  `k = 0` (the default) is the
        main diagonal, `k < 0` is below it and `k > 0` is above.

    Returns
    -------
    L : ndarray, shape (M, N)
        Lower triangle of `m`, of same shape and data-type as `m`.

    See Also
    --------
    triu : same thing, only for the upper triangle

    Examples
    --------
    >>> np.tril([[1,2,3],[4,5,6],[7,8,9],[10,11,12]], -1)
    array([[ 0,  0,  0],
           [ 4,  0,  0],
           [ 7,  8,  0],
           [10, 11, 12]])
    
    """
    raise NotImplementedError


def tril_indices(n, k=0):
    """Return the indices for the lower-triangle of an (n, n) array.

    Parameters
    ----------
    n : int
        The row dimension of the square arrays for which the returned
        indices will be valid.
    k : int, optional
        Diagonal offset (see `tril` for details).

    Returns
    -------
    inds : tuple of arrays
        The indices for the triangle. The returned tuple contains two arrays,
        each with the indices along one dimension of the array.

    See also
    --------
    triu_indices : similar function, for upper-triangular.
    mask_indices : generic function accepting an arbitrary mask function.
    tril, triu

    Notes
    -----
    .. versionadded:: 1.4.0

    Examples
    --------
    Compute two different sets of indices to access 4x4 arrays, one for the
    lower triangular part starting at the main diagonal, and one starting two
    diagonals further right:

    >>> il1 = np.tril_indices(4)
    >>> il2 = np.tril_indices(4, 2)

    Here is how they can be used with a sample array:

    >>> a = np.arange(16).reshape(4, 4)
    >>> a
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11],
           [12, 13, 14, 15]])

    Both for indexing:

    >>> a[il1]
    array([ 0,  4,  5,  8,  9, 10, 12, 13, 14, 15])

    And for assigning values:

    >>> a[il1] = -1
    >>> a
    array([[-1,  1,  2,  3],
           [-1, -1,  6,  7],
           [-1, -1, -1, 11],
           [-1, -1, -1, -1]])

    These cover almost the whole array (two diagonals right of the main one):

    >>> a[il2] = -10
    >>> a
    array([[-10, -10, -10,   3],
           [-10, -10, -10, -10],
           [-10, -10, -10, -10],
           [-10, -10, -10, -10]])
    
    """
    raise NotImplementedError


def tril_indices_from(arr, k=0):
    """Return the indices for the lower-triangle of arr.

    See `tril_indices` for full details.

    Parameters
    ----------
    arr : array_like
        The indices will be valid for square arrays whose dimensions are
        the same as arr.
    k : int, optional
        Diagonal offset (see `tril` for details).

    See Also
    --------
    tril_indices, tril

    Notes
    -----
    .. versionadded:: 1.4.0
    
    """
    raise NotImplementedError


def trim_zeros(filt, trim='fb'):
    """Trim the leading and/or trailing zeros from a 1-D array or sequence.

    Parameters
    ----------
    filt : 1-D array or sequence
        Input array.
    trim : str, optional
        A string with 'f' representing trim from front and 'b' to trim from
        back. Default is 'fb', trim zeros from both front and back of the
        array.

    Returns
    -------
    trimmed : 1-D array or sequence
        The result of trimming the input. The input data type is preserved.

    Examples
    --------
    >>> a = np.array((0, 0, 0, 1, 2, 3, 0, 2, 1, 0))
    >>> np.trim_zeros(a)
    array([1, 2, 3, 0, 2, 1])

    >>> np.trim_zeros(a, 'b')
    array([0, 0, 0, 1, 2, 3, 0, 2, 1])

    The input data type is preserved, list/tuple in means list/tuple out.

    >>> np.trim_zeros([0, 1, 2, 0])
    [1, 2]
    
    """
    raise NotImplementedError


def triu(m, k=0):
    """Upper triangle of an array.

    Return a copy of a matrix with the elements below the `k`-th diagonal
    zeroed.

    Please refer to the documentation for `tril` for further details.

    See Also
    --------
    tril : lower triangle of an array

    Examples
    --------
    >>> np.triu([[1,2,3],[4,5,6],[7,8,9],[10,11,12]], -1)
    array([[ 1,  2,  3],
           [ 4,  5,  6],
           [ 0,  8,  9],
           [ 0,  0, 12]])
    
    """
    raise NotImplementedError


def triu_indices(n, k=0):
    """Return the indices for the upper-triangle of an (n, n) array.

    Parameters
    ----------
    n : int
        The size of the arrays for which the returned indices will
        be valid.
    k : int, optional
        Diagonal offset (see `triu` for details).

    Returns
    -------
    inds : tuple of arrays
        The indices for the triangle. The returned tuple contains two arrays,
        each with the indices along one dimension of the array.

    See also
    --------
    tril_indices : similar function, for lower-triangular.
    mask_indices : generic function accepting an arbitrary mask function.
    triu, tril

    Notes
    -----
    .. versionadded:: 1.4.0

    Examples
    --------
    Compute two different sets of indices to access 4x4 arrays, one for the
    upper triangular part starting at the main diagonal, and one starting two
    diagonals further right:

    >>> iu1 = np.triu_indices(4)
    >>> iu2 = np.triu_indices(4, 2)

    Here is how they can be used with a sample array:

    >>> a = np.arange(16).reshape(4, 4)
    >>> a
    array([[ 0,  1,  2,  3],
           [ 4,  5,  6,  7],
           [ 8,  9, 10, 11],
           [12, 13, 14, 15]])

    Both for indexing:

    >>> a[iu1]
    array([ 0,  1,  2,  3,  5,  6,  7, 10, 11, 15])

    And for assigning values:

    >>> a[iu1] = -1
    >>> a
    array([[-1, -1, -1, -1],
           [ 4, -1, -1, -1],
           [ 8,  9, -1, -1],
           [12, 13, 14, -1]])

    These cover only a small part of the whole array (two diagonals right
    of the main one):

    >>> a[iu2] = -10
    >>> a
    array([[ -1,  -1, -10, -10],
           [  4,  -1,  -1, -10],
           [  8,   9,  -1,  -1],
           [ 12,  13,  14,  -1]])
    
    """
    raise NotImplementedError


def triu_indices_from(arr, k=0):
    """Return the indices for the upper-triangle of an (n, n) array.

    See `triu_indices` for full details.

    Parameters
    ----------
    arr : array_like
        The indices will be valid for square arrays whose dimensions are
        the same as arr.
    k : int, optional
      Diagonal offset (see `triu` for details).

    See Also
    --------
    triu_indices, triu

    Notes
    -----
    .. versionadded:: 1.4.0
    
    """
    raise NotImplementedError


def typename(char_):
    """Return a description for the given data type code.

    Parameters
    ----------
    char : str
        Data type code.

    Returns
    -------
    out : str
        Description of the input data type code.

    See Also
    --------
    dtype, typecodes

    Examples
    --------
    >>> typechars = ['S1', '?', 'B', 'D', 'G', 'F', 'I', 'H', 'L', 'O', 'Q',
    ...              'S', 'U', 'V', 'b', 'd', 'g', 'f', 'i', 'h', 'l', 'q']
    >>> for typechar in typechars:
    ...     print typechar, ' : ', np.typename(typechar)
    ...
    S1  :  character
    ?  :  bool
    B  :  unsigned char
    D  :  complex double precision
    G  :  complex long double precision
    F  :  complex single precision
    I  :  unsigned integer
    H  :  unsigned short
    L  :  unsigned long integer
    O  :  object
    Q  :  unsigned long long integer
    S  :  string
    U  :  unicode
    V  :  void
    b  :  signed char
    d  :  double precision
    g  :  long precision
    f  :  single precision
    i  :  integer
    h  :  short
    l  :  long integer
    q  :  long long integer
    
    """
    raise NotImplementedError


def union1d(ar1, ar2):
    """Find the union of two arrays.

    Return the unique, sorted array of values that are in either of the two
    input arrays.

    Parameters
    ----------
    ar1, ar2 : array_like
        Input arrays. They are flattened if they are not already 1D.

    Returns
    -------
    union : ndarray
        Unique, sorted union of the input arrays.

    See Also
    --------
    numpy.lib.arraysetops : Module with a number of other functions for
                            performing set operations on arrays.

    Examples
    --------
    >>> np.union1d([-1, 0, 1], [-2, 0, 2])
    array([-2, -1,  0,  1,  2])
    
    """
    raise NotImplementedError


def unique(ar, return_index=False, return_inverse=False):
    """Find the unique elements of an array.

    Returns the sorted unique elements of an array. There are two optional
    outputs in addition to the unique elements: the indices of the input array
    that give the unique values, and the indices of the unique array that
    reconstruct the input array.

    Parameters
    ----------
    ar : array_like
        Input array. This will be flattened if it is not already 1-D.
    return_index : bool, optional
        If True, also return the indices of `ar` that result in the unique
        array.
    return_inverse : bool, optional
        If True, also return the indices of the unique array that can be used
        to reconstruct `ar`.

    Returns
    -------
    unique : ndarray
        The sorted unique values.
    unique_indices : ndarray, optional
        The indices of the unique values in the (flattened) original array.
        Only provided if `return_index` is True.
    unique_inverse : ndarray, optional
        The indices to reconstruct the (flattened) original array from the
        unique array. Only provided if `return_inverse` is True.

    See Also
    --------
    numpy.lib.arraysetops : Module with a number of other functions for
                            performing set operations on arrays.

    Examples
    --------
    >>> np.unique([1, 1, 2, 2, 3, 3])
    array([1, 2, 3])
    >>> a = np.array([[1, 1], [2, 3]])
    >>> np.unique(a)
    array([1, 2, 3])

    Return the indices of the original array that give the unique values:

    >>> a = np.array(['a', 'b', 'b', 'c', 'a'])
    >>> u, indices = np.unique(a, return_index=True)
    >>> u
    array(['a', 'b', 'c'],
           dtype='|S1')
    >>> indices
    array([0, 1, 3])
    >>> a[indices]
    array(['a', 'b', 'c'],
           dtype='|S1')

    Reconstruct the input array from the unique values:

    >>> a = np.array([1, 2, 6, 4, 2, 3, 2])
    >>> u, indices = np.unique(a, return_inverse=True)
    >>> u
    array([1, 2, 3, 4, 6])
    >>> indices
    array([0, 1, 4, 3, 1, 2, 1])
    >>> u[indices]
    array([1, 2, 6, 4, 2, 3, 2])
    
    """
    raise NotImplementedError


def unpackbits(myarray, axis=None):
    """unpackbits(myarray, axis=None)

    Unpacks elements of a uint8 array into a binary-valued output array.

    Each element of `myarray` represents a bit-field that should be unpacked
    into a binary-valued output array. The shape of the output array is either
    1-D (if `axis` is None) or the same shape as the input array with unpacking
    done along the axis specified.

    Parameters
    ----------
    myarray : ndarray, uint8 type
       Input array.
    axis : int, optional
       Unpacks along this axis.

    Returns
    -------
    unpacked : ndarray, uint8 type
       The elements are binary-valued (0 or 1).

    See Also
    --------
    packbits : Packs the elements of a binary-valued array into bits in a uint8
               array.

    Examples
    --------
    >>> a = np.array([[2], [7], [23]], dtype=np.uint8)
    >>> a
    array([[ 2],
           [ 7],
           [23]], dtype=uint8)
    >>> b = np.unpackbits(a, axis=1)
    >>> b
    array([[0, 0, 0, 0, 0, 0, 1, 0],
           [0, 0, 0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 0, 1, 1, 1]], dtype=uint8)
    
    """
    # BUILTIN
    raise NotImplementedError


def unravel_index(indices, dims, order='C'):
    """unravel_index(indices, dims, order='C')

    Converts a flat index or array of flat indices into a tuple
    of coordinate arrays.

    Parameters
    ----------
    indices : array_like
        An integer array whose elements are indices into the flattened
        version of an array of dimensions ``dims``. Before version 1.6.0,
        this function accepted just one index value.
    dims : tuple of ints
        The shape of the array to use for unraveling ``indices``.
    order : {'C', 'F'}, optional
        .. versionadded:: 1.6.0

        Determines whether the indices should be viewed as indexing in
        C (row-major) order or FORTRAN (column-major) order.

    Returns
    -------
    unraveled_coords : tuple of ndarray
        Each array in the tuple has the same shape as the ``indices``
        array.

    See Also
    --------
    ravel_multi_index

    Examples
    --------
    >>> np.unravel_index([22, 41, 37], (7,6))
    (array([3, 6, 6]), array([4, 5, 1]))
    >>> np.unravel_index([31, 41, 13], (7,6), order='F')
    (array([3, 6, 6]), array([4, 5, 1]))

    >>> np.unravel_index(1621, (6,7,8,9))
    (3, 1, 4, 1)
    
    """
    # BUILTIN
    raise NotImplementedError


def unwrap(p, discont=3.14159265359, axis=-1):
    """Unwrap by changing deltas between values to 2*pi complement.

    Unwrap radian phase `p` by changing absolute jumps greater than
    `discont` to their 2*pi complement along the given axis.

    Parameters
    ----------
    p : array_like
        Input array.
    discont : float, optional
        Maximum discontinuity between values, default is ``pi``.
    axis : int, optional
        Axis along which unwrap will operate, default is the last axis.

    Returns
    -------
    out : ndarray
        Output array.

    See Also
    --------
    rad2deg, deg2rad

    Notes
    -----
    If the discontinuity in `p` is smaller than ``pi``, but larger than
    `discont`, no unwrapping is done because taking the 2*pi complement
    would only make the discontinuity larger.

    Examples
    --------
    >>> phase = np.linspace(0, np.pi, num=5)
    >>> phase[3:] += np.pi
    >>> phase
    array([ 0.        ,  0.78539816,  1.57079633,  5.49778714,  6.28318531])
    >>> np.unwrap(phase)
    array([ 0.        ,  0.78539816,  1.57079633, -0.78539816,  0.        ])
    
    """
    raise NotImplementedError


def vander(x, N=None):
    """Generate a Van der Monde matrix.

    The columns of the output matrix are decreasing powers of the input
    vector.  Specifically, the `i`-th output column is the input vector
    raised element-wise to the power of ``N - i - 1``.  Such a matrix with
    a geometric progression in each row is named for Alexandre-Theophile
    Vandermonde.

    Parameters
    ----------
    x : array_like
        1-D input array.
    N : int, optional
        Order of (number of columns in) the output.  If `N` is not specified,
        a square array is returned (``N = len(x)``).

    Returns
    -------
    out : ndarray
        Van der Monde matrix of order `N`.  The first column is ``x^(N-1)``,
        the second ``x^(N-2)`` and so forth.

    Examples
    --------
    >>> x = np.array([1, 2, 3, 5])
    >>> N = 3
    >>> np.vander(x, N)
    array([[ 1,  1,  1],
           [ 4,  2,  1],
           [ 9,  3,  1],
           [25,  5,  1]])

    >>> np.column_stack([x**(N-1-i) for i in range(N)])
    array([[ 1,  1,  1],
           [ 4,  2,  1],
           [ 9,  3,  1],
           [25,  5,  1]])

    >>> x = np.array([1, 2, 3, 5])
    >>> np.vander(x)
    array([[  1,   1,   1,   1],
           [  8,   4,   2,   1],
           [ 27,   9,   3,   1],
           [125,  25,   5,   1]])

    The determinant of a square Vandermonde matrix is the product
    of the differences between the values of the input vector:

    >>> np.linalg.det(np.vander(x))
    48.000000000000043
    >>> (5-3)*(5-2)*(5-1)*(3-2)*(3-1)*(2-1)
    48
    
    """
    raise NotImplementedError


def var(a, axis=None, dtype=None, out=None, ddof=0):
    """Compute the variance along the specified axis.

    Returns the variance of the array elements, a measure of the spread of a
    distribution.  The variance is computed for the flattened array by
    default, otherwise over the specified axis.

    Parameters
    ----------
    a : array_like
        Array containing numbers whose variance is desired.  If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the variance is computed.  The default is to compute
        the variance of the flattened array.
    dtype : data-type, optional
        Type to use in computing the variance.  For arrays of integer type
        the default is `float32`; for arrays of float types it is the same as
        the array type.
    out : ndarray, optional
        Alternate output array in which to place the result.  It must have
        the same shape as the expected output, but the type is cast if
        necessary.
    ddof : int, optional
        "Delta Degrees of Freedom": the divisor used in the calculation is
        ``N - ddof``, where ``N`` represents the number of elements. By
        default `ddof` is zero.

    Returns
    -------
    variance : ndarray, see dtype parameter above
        If ``out=None``, returns a new array containing the variance;
        otherwise, a reference to the output array is returned.

    See Also
    --------
    std : Standard deviation
    mean : Average
    numpy.doc.ufuncs : Section "Output arguments"

    Notes
    -----
    The variance is the average of the squared deviations from the mean,
    i.e.,  ``var = mean(abs(x - x.mean())**2)``.

    The mean is normally calculated as ``x.sum() / N``, where ``N = len(x)``.
    If, however, `ddof` is specified, the divisor ``N - ddof`` is used
    instead.  In standard statistical practice, ``ddof=1`` provides an
    unbiased estimator of the variance of a hypothetical infinite population.
    ``ddof=0`` provides a maximum likelihood estimate of the variance for
    normally distributed variables.

    Note that for complex numbers, the absolute value is taken before
    squaring, so that the result is always real and nonnegative.

    For floating-point input, the variance is computed using the same
    precision the input has.  Depending on the input data, this can cause
    the results to be inaccurate, especially for `float32` (see example
    below).  Specifying a higher-accuracy accumulator using the ``dtype``
    keyword can alleviate this issue.

    Examples
    --------
    >>> a = np.array([[1,2],[3,4]])
    >>> np.var(a)
    1.25
    >>> np.var(a,0)
    array([ 1.,  1.])
    >>> np.var(a,1)
    array([ 0.25,  0.25])

    In single precision, var() can be inaccurate:

    >>> a = np.zeros((2,512*512), dtype=np.float32)
    >>> a[0,:] = 1.0
    >>> a[1,:] = 0.1
    >>> np.var(a)
    0.20405951142311096

    Computing the standard deviation in float64 is more accurate:

    >>> np.var(a, dtype=np.float64)
    0.20249999932997387
    >>> ((1-0.55)**2 + (0.1-0.55)**2)/2
    0.20250000000000001
    
    """
    raise NotImplementedError


def vdot():
    """Return the dot product of two vectors.

    The vdot(`a`, `b`) function handles complex numbers differently than
    dot(`a`, `b`).  If the first argument is complex the complex conjugate
    of the first argument is used for the calculation of the dot product.

    Note that `vdot` handles multidimensional arrays differently than `dot`:
    it does *not* perform a matrix product, but flattens input arguments
    to 1-D vectors first. Consequently, it should only be used for vectors.

    Parameters
    ----------
    a : array_like
        If `a` is complex the complex conjugate is taken before calculation
        of the dot product.
    b : array_like
        Second argument to the dot product.

    Returns
    -------
    output : ndarray
        Dot product of `a` and `b`.  Can be an int, float, or
        complex depending on the types of `a` and `b`.

    See Also
    --------
    dot : Return the dot product without using the complex conjugate of the
          first argument.

    Examples
    --------
    >>> a = np.array([1+2j,3+4j])
    >>> b = np.array([5+6j,7+8j])
    >>> np.vdot(a, b)
    (70-8j)
    >>> np.vdot(b, a)
    (70+8j)

    Note that higher-dimensional arrays are flattened!

    >>> a = np.array([[1, 4], [5, 6]])
    >>> b = np.array([[4, 1], [2, 2]])
    >>> np.vdot(a, b)
    30
    >>> np.vdot(b, a)
    30
    >>> 1*4 + 4*1 + 5*2 + 6*2
    30
    
    """
    # BUILTIN
    raise NotImplementedError


def vsplit(ary, indices_or_sections):
    """Split an array into multiple sub-arrays vertically (row-wise).

    Please refer to the ``split`` documentation.  ``vsplit`` is equivalent
    to ``split`` with `axis=0` (default), the array is always split along the
    first axis regardless of the array dimension.

    See Also
    --------
    split : Split an array into multiple sub-arrays of equal size.

    Examples
    --------
    >>> x = np.arange(16.0).reshape(4, 4)
    >>> x
    array([[  0.,   1.,   2.,   3.],
           [  4.,   5.,   6.,   7.],
           [  8.,   9.,  10.,  11.],
           [ 12.,  13.,  14.,  15.]])
    >>> np.vsplit(x, 2)
    [array([[ 0.,  1.,  2.,  3.],
           [ 4.,  5.,  6.,  7.]]),
     array([[  8.,   9.,  10.,  11.],
           [ 12.,  13.,  14.,  15.]])]
    >>> np.vsplit(x, np.array([3, 6]))
    [array([[  0.,   1.,   2.,   3.],
           [  4.,   5.,   6.,   7.],
           [  8.,   9.,  10.,  11.]]),
     array([[ 12.,  13.,  14.,  15.]]),
     array([], dtype=float64)]

    With a higher dimensional array the split is still along the first axis.

    >>> x = np.arange(8.0).reshape(2, 2, 2)
    >>> x
    array([[[ 0.,  1.],
            [ 2.,  3.]],
           [[ 4.,  5.],
            [ 6.,  7.]]])
    >>> np.vsplit(x, 2)
    [array([[[ 0.,  1.],
            [ 2.,  3.]]]),
     array([[[ 4.,  5.],
            [ 6.,  7.]]])]
    
    """
    raise NotImplementedError


def vstack(tup):
    """Stack arrays in sequence vertically (row wise).

    Take a sequence of arrays and stack them vertically to make a single
    array. Rebuild arrays divided by `vsplit`.

    Parameters
    ----------
    tup : sequence of ndarrays
        Tuple containing arrays to be stacked. The arrays must have the same
        shape along all but the first axis.

    Returns
    -------
    stacked : ndarray
        The array formed by stacking the given arrays.

    See Also
    --------
    hstack : Stack arrays in sequence horizontally (column wise).
    dstack : Stack arrays in sequence depth wise (along third dimension).
    concatenate : Join a sequence of arrays together.
    vsplit : Split array into a list of multiple sub-arrays vertically.


    Notes
    -----
    Equivalent to ``np.concatenate(tup, axis=0)``

    Examples
    --------
    >>> a = np.array([1, 2, 3])
    >>> b = np.array([2, 3, 4])
    >>> np.vstack((a,b))
    array([[1, 2, 3],
           [2, 3, 4]])

    >>> a = np.array([[1], [2], [3]])
    >>> b = np.array([[2], [3], [4]])
    >>> np.vstack((a,b))
    array([[1],
           [2],
           [3],
           [2],
           [3],
           [4]])
    
    """
    raise NotImplementedError


def where():
    """where(condition, [x, y])

    Return elements, either from `x` or `y`, depending on `condition`.

    If only `condition` is given, return ``condition.nonzero()``.

    Parameters
    ----------
    condition : array_like, bool
        When True, yield `x`, otherwise yield `y`.
    x, y : array_like, optional
        Values from which to choose. `x` and `y` need to have the same
        shape as `condition`.

    Returns
    -------
    out : ndarray or tuple of ndarrays
        If both `x` and `y` are specified, the output array contains
        elements of `x` where `condition` is True, and elements from
        `y` elsewhere.

        If only `condition` is given, return the tuple
        ``condition.nonzero()``, the indices where `condition` is True.

    See Also
    --------
    nonzero, choose

    Notes
    -----
    If `x` and `y` are given and input arrays are 1-D, `where` is
    equivalent to::

        [xv if c else yv for (c,xv,yv) in zip(condition,x,y)]

    Examples
    --------
    >>> np.where([[True, False], [True, True]],
    ...          [[1, 2], [3, 4]],
    ...          [[9, 8], [7, 6]])
    array([[1, 8],
           [3, 4]])

    >>> np.where([[0, 1], [1, 0]])
    (array([0, 1]), array([1, 0]))

    >>> x = np.arange(9.).reshape(3, 3)
    >>> np.where( x > 5 )
    (array([2, 2, 2]), array([0, 1, 2]))
    >>> x[np.where( x > 3.0 )]               # Note: result is 1D.
    array([ 4.,  5.,  6.,  7.,  8.])
    >>> np.where(x < 5, x, -1)               # Note: broadcasting.
    array([[ 0.,  1.,  2.],
           [ 3.,  4., -1.],
           [-1., -1., -1.]])
    
    """
    # BUILTIN
    raise NotImplementedError


def who(vardict=None):
    """Print the Numpy arrays in the given dictionary.

    If there is no dictionary passed in or `vardict` is None then returns
    Numpy arrays in the globals() dictionary (all Numpy arrays in the
    namespace).

    Parameters
    ----------
    vardict : dict, optional
        A dictionary possibly containing ndarrays.  Default is globals().

    Returns
    -------
    out : None
        Returns 'None'.

    Notes
    -----
    Prints out the name, shape, bytes and type of all of the ndarrays present
    in `vardict`.

    Examples
    --------
    >>> a = np.arange(10)
    >>> b = np.ones(20)
    >>> np.who()
    Name            Shape            Bytes            Type
    ===========================================================
    a               10               40               int32
    b               20               160              float64
    Upper bound on total bytes  =       200

    >>> d = {'x': np.arange(2.0), 'y': np.arange(3.0), 'txt': 'Some str',
    ... 'idx':5}
    >>> np.who(d)
    Name            Shape            Bytes            Type
    ===========================================================
    y               3                24               float64
    x               2                16               float64
    Upper bound on total bytes  =       40
    
    """
    raise NotImplementedError

