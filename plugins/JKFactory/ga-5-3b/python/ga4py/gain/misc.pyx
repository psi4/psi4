# cython: profile=True
import math 

from ga4py import ga
from core import asarray
from core import flatiter
from core import get_dtype
from core import get_shape
from core import is_array
from core import is_distributed
from core import multiply
from core import ndarray
from core import _npin_piece_based_on_out
from core import should_distribute
from core import sync
from core import zeros

import numpy as np
cimport numpy as np

def zeros_like(a, dtype=None, order='K', subok=True):
    """Return an array of zeros with the same shape and type as a given array.

    Equivalent to ``a.copy().fill(0)``.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define the parameters of
        the returned array.

    Returns
    -------
    out : ndarray
        Array of zeros with same shape and type as `a`.

    See Also
    --------
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    zeros : Return a new array setting values to zero.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.

    Examples
    --------
    >>> x = np.arange(6)
    >>> x = x.reshape((2, 3))
    >>> x
    array([[0, 1, 2],
           [3, 4, 5]])
    >>> np.zeros_like(x)
    array([[0, 0, 0],
           [0, 0, 0]])

    >>> y = np.arange(3, dtype=np.float)
    >>> y
    array([ 0.,  1.,  2.])
    >>> np.zeros_like(y)
    array([ 0.,  0.,  0.])

    """
    return zeros(a.shape, dtype=a.dtype)

def ones(shape, dtype=np.float, order='C'):
    """Return a new array of given shape and type, filled with ones.

    Please refer to the documentation for `zeros`.

    See Also
    --------
    zeros

    Examples
    --------
    >>> np.ones(5)
    array([ 1.,  1.,  1.,  1.,  1.])

    >>> np.ones((5,), dtype=np.int)
    array([1, 1, 1, 1, 1])

    >>> np.ones((2, 1))
    array([[ 1.],
           [ 1.]])

    >>> s = (2,2)
    >>> np.ones(s)
    array([[ 1.,  1.],
           [ 1.,  1.]])
    
    """
    if not should_distribute(shape):
        return np.ones(shape, dtype, order)
    a = ndarray(shape, dtype)
    buf = a.access()
    if buf is not None:
        buf[:] = 1
        a.release_update()
    return a

def ones_like(x):
    """ones_like(x[, out])

    Returns an array of ones with the same shape and type as a given array.

    Equivalent to ``a.copy().fill(1)``.

    Please refer to the documentation for `zeros_like`.

    See Also
    --------
    zeros_like

    Examples
    --------
    >>> a = np.array([[1, 2, 3], [4, 5, 6]])
    >>> np.ones_like(a)
    array([[1, 1, 1],
           [1, 1, 1]])

    """
    return ones(x.shape, dtype=x.dtype)

def empty(shape, dtype=float, order='C'):
    """empty(shape, dtype=float, order='C')

    Return a new array of given shape and type, without initializing entries.

    Parameters
    ----------
    shape : int or tuple of int
        Shape of the empty array
    dtype : data-type, optional
        Desired output data-type.
    order : {'C', 'F'}, optional
        Whether to store multi-dimensional data in C (row-major) or
        Fortran (column-major) order in memory.

    See Also
    --------
    empty_like, zeros, ones

    Notes
    -----
    `empty`, unlike `zeros`, does not set the array values to zero,
    and may therefore be marginally faster.  On the other hand, it requires
    the user to manually set all the values in the array, and should be
    used with caution.

    Examples
    --------
    >>> np.empty([2, 2])
    array([[ -9.74499359e+001,   6.69583040e-309],  #random data
           [  2.13182611e-314,   3.06959433e-309]])

    >>> np.empty([2, 2], dtype=int)
    array([[-1073741821, -1067949133],  #random data
           [  496041986,    19249760]])

    """
    return ndarray(shape, dtype)

def empty_like(a, dtype=None, order='K', subok=True):
    """    Return a new array with the same shape and type as a given array.

    Parameters
    ----------
    a : array_like
        The shape and data-type of `a` define the parameters of the
        returned array.

    Returns
    -------
    out : ndarray
        Array of random data with the same shape and type as `a`.

    See Also
    --------
    ones_like : Return an array of ones with shape and type of input.
    zeros_like : Return an array of zeros with shape and type of input.
    empty : Return a new uninitialized array.
    ones : Return a new array setting values to one.
    zeros : Return a new array setting values to zero.

    Notes
    -----
    This function does *not* initialize the returned array; to do that use
    `zeros_like` or `ones_like` instead. It may be marginally faster than the
    functions that do set the array values.

    Examples
    --------
    >>> a = ([1,2,3], [4,5,6])                         # a is array-like
    >>> np.empty_like(a)
    array([[-1073741821, -1073741821,           3],    #random
           [          0,           0, -1073741821]])
    >>> a = np.array([[1., 2., 3.],[4.,5.,6.]])
    >>> np.empty_like(a)
    array([[ -2.00000715e+000,   1.48219694e-323,  -2.00000572e+000], #random
           [  4.38791518e-305,  -2.00000715e+000,   4.17269252e-309]])
    
    """
    return empty(a.shape, dtype or a.dtype)

def eye(N, M=None, k=0, dtype=float):
    """Return a 2-D array with ones on the diagonal and zeros elsewhere.

    Parameters
    ----------
    N : int
      Number of rows in the output.
    M : int, optional
      Number of columns in the output. If None, defaults to `N`.
    k : int, optional
      Index of the diagonal: 0 refers to the main diagonal, a positive value
      refers to an upper diagonal, and a negative value to a lower diagonal.
    dtype : dtype, optional
      Data-type of the returned array.

    Returns
    -------
    I : ndarray (N,M)
      An array where all elements are equal to zero, except for the `k`-th
      diagonal, whose values are equal to one.

    See Also
    --------
    diag : Return a diagonal 2-D array using a 1-D array specified by the user.

    Examples
    --------
    >>> np.eye(2, dtype=int)
    array([[1, 0],
           [0, 1]])
    >>> np.eye(3, k=1)
    array([[ 0.,  1.,  0.],
           [ 0.,  0.,  1.],
           [ 0.,  0.,  0.]])
    
    """
    if M is None:
        M = N
    if not should_distribute((N,M)):
        return np.eye(N,M,k,dtype)
    a = zeros((N,M), dtype=dtype)
    nda = a.access()
    if nda is not None:
        lo,hi = a.distribution()
        indices = np.indices(nda.shape)
        indices[0] += lo[0]
        indices[1] += lo[1]-k
        bindex = (indices[0] == indices[1])
        nda[bindex] = 1
        a.release_update()
    return a

def identity(n, dtype=None):
    """Return the identity array.

    The identity array is a square array with ones on
    the main diagonal.

    Parameters
    ----------
    n : int
        Number of rows (and columns) in `n` x `n` output.
    dtype : data-type, optional
        Data-type of the output.  Defaults to ``float``.

    Returns
    -------
    out : ndarray
        `n` x `n` array with its main diagonal set to one,
        and all other elements 0.

    Examples
    --------
    >>> np.identity(3)
    array([[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]])

    """
    if dtype is None:
        dtype = np.dtype(float)
    return eye(n,n,dtype=dtype)

def fromfunction(func, shape, **kwargs):
    """Construct an array by executing a function over each coordinate.

    The resulting array therefore has a value ``fn(x, y, z)`` at
    coordinate ``(x, y, z)``.

    Parameters
    ----------
    function : callable
        The function is called with N parameters, each of which
        represents the coordinates of the array varying along a
        specific axis.  For example, if `shape` were ``(2, 2)``, then
        the parameters would be two arrays, ``[[0, 0], [1, 1]]`` and
        ``[[0, 1], [0, 1]]``.  `function` must be capable of operating on
        arrays, and should return a scalar value.
    shape : (N,) tuple of ints
        Shape of the output array, which also determines the shape of
        the coordinate arrays passed to `function`.
    dtype : data-type, optional
        Data-type of the coordinate arrays passed to `function`.
        By default, `dtype` is float.

    Returns
    -------
    out : any
        The result of the call to `function` is passed back directly.
        Therefore the type and shape of `out` is completely determined by
        `function`.

    See Also
    --------
    indices, meshgrid

    Notes
    -----
    Keywords other than `shape` and `dtype` are passed to `function`.

    Examples
    --------
    >>> np.fromfunction(lambda i, j: i == j, (3, 3), dtype=int)
    array([[ True, False, False],
           [False,  True, False],
           [False, False,  True]], dtype=bool)

    >>> np.fromfunction(lambda i, j: i + j, (3, 3), dtype=int)
    array([[0, 1, 2],
           [1, 2, 3],
           [2, 3, 4]])
    
    """
    if not should_distribute(shape):
        return np.fromfunction(func, shape, **kwargs)
    dtype = kwargs.pop('dtype', np.float32)
    # create the new GA (collective operation)
    a = ndarray(shape, dtype)
    # determine which part of 'a' we maintain
    local_array = ga.access(a.handle)
    if local_array is not None:
        lo,hi = a.distribution()
        local_shape = hi-lo
        # create a numpy indices array
        args = np.indices(local_shape, dtype=dtype)
        # modify the indices arrays based on our distribution
        for index in xrange(len(lo)):
            args[index] += lo[index]
        # call the passed function
        buf = func(*args, **kwargs)
        # now put the data into the global array
        local_array[:] = buf
    #sync()
    return a

def arange(start, stop=None, step=None, dtype=None, shape=None):
    """Return evenly spaced values within a given interval.

    Values are generated within the half-open interval ``[start, stop)``
    (in other words, the interval including `start` but excluding `stop`).
    For integer arguments the function is equivalent to the Python built-in
    `range <http://docs.python.org/lib/built-in-funcs.html>`_ function,
    but returns a ndarray rather than a list.

    Parameters
    ----------
    start : number, optional
        Start of interval.  The interval includes this value.  The default
        start value is 0.
    stop : number
        End of interval.  The interval does not include this value.
    step : number, optional
        Spacing between values.  For any output `out`, this is the distance
        between two adjacent values, ``out[i+1] - out[i]``.  The default
        step size is 1.  If `step` is specified, `start` must also be given.
    dtype : dtype
        The type of the output array.  If `dtype` is not given, infer the data
        type from the other input arguments.

    Hack Parameters
    ---------------
    shape : tuple of ints
        Useful shortcut when doing something like np.arange(...).reshape(...)

    Returns
    -------
    out : ndarray
        Array of evenly spaced values.

        For floating point arguments, the length of the result is
        ``ceil((stop - start)/step)``.  Because of floating point overflow,
        this rule may result in the last element of `out` being greater
        than `stop`.

    See Also
    --------
    linspace : Evenly spaced numbers with careful handling of endpoints.
    ogrid: Arrays of evenly spaced numbers in N-dimensions
    mgrid: Grid-shaped arrays of evenly spaced numbers in N-dimensions

    Examples
    --------
    >>> np.arange(3)
    array([0, 1, 2])
    >>> np.arange(3.0)
    array([ 0.,  1.,  2.])
    >>> np.arange(3,7)
    array([3, 4, 5, 6])
    >>> np.arange(3,7,2)
    array([3, 5])
    
    """
    if step == 0:
        raise ValueError, "step size of 0 not allowed"
    if not step:
        step = 1
    if not stop:
        start,stop = 0,start
    length = 0
    if ((step < 0 and stop >= start) or (step > 0 and start >= stop)):
        length = 0
    else:
        # true division, otherwise off by one
        length = int(math.ceil((stop-start)/step))
    # bail if threshold not met
    if not should_distribute(length):
        return np.arange(start,stop,step,dtype)
    if dtype is None:
        if (isinstance(start, (int,long))
                and isinstance(stop, (int,long))
                and isinstance(step, (int,long))):
            dtype = np.int64
        else:
            dtype = np.float64
    a = None
    if shape is not None:
        shape = np.asarray(shape,dtype=np.int64)
        if np.prod(shape) != length:
            raise ValueError, "total size of new array must be unchanged"
        a = ndarray(shape, dtype)
        a_local = a.access()
        if a_local is not None:
            lo,hi = a.distribution()
            lshape = hi-lo
            v = np.add.reduce(
                    (np.indices(lshape).reshape(len(lshape),-1).T + lo)
                    * (np.asarray(a.strides)/a.itemsize), axis=1)
            a_local.flat = v*step + start
            a.release_update()
    else:
        a = ndarray(length, dtype)
        a_local = a.access()
        if a_local is not None:
            lo,hi = a.distribution()
            a_local[...] = np.arange(lo[0],hi[0])
            a_local *= step
            a_local += start
            a.release_update()
    return a

def linspace(start, stop, num=50, endpoint=True, retstep=False):
    """Return evenly spaced numbers over a specified interval.

    Returns `num` evenly spaced samples, calculated over the
    interval [`start`, `stop` ].

    The endpoint of the interval can optionally be excluded.

    Parameters
    ----------
    start : scalar
        The starting value of the sequence.
    stop : scalar
        The end value of the sequence, unless `endpoint` is set to False.
        In that case, the sequence consists of all but the last of ``num + 1``
        evenly spaced samples, so that `stop` is excluded.  Note that the step
        size changes when `endpoint` is False.
    num : int, optional
        Number of samples to generate. Default is 50.
    endpoint : bool, optional
        If True, `stop` is the last sample. Otherwise, it is not included.
        Default is True.
    retstep : bool, optional
        If True, return (`samples`, `step`), where `step` is the spacing
        between samples.

    Returns
    -------
    samples : ndarray
        There are `num` equally spaced samples in the closed interval
        ``[start, stop]`` or the half-open interval ``[start, stop)``
        (depending on whether `endpoint` is True or False).
    step : float (only if `retstep` is True)
        Size of spacing between samples.


    See Also
    --------
    arange : Similiar to `linspace`, but uses a step size (instead of the
             number of samples).
    logspace : Samples uniformly distributed in log space.

    Examples
    --------
    >>> np.linspace(2.0, 3.0, num=5)
        array([ 2.  ,  2.25,  2.5 ,  2.75,  3.  ])
    >>> np.linspace(2.0, 3.0, num=5, endpoint=False)
        array([ 2. ,  2.2,  2.4,  2.6,  2.8])
    >>> np.linspace(2.0, 3.0, num=5, retstep=True)
        (array([ 2.  ,  2.25,  2.5 ,  2.75,  3.  ]), 0.25)

    Graphical illustration:

    >>> import matplotlib.pyplot as plt
    >>> N = 8
    >>> y = np.zeros(N)
    >>> x1 = np.linspace(0, 10, N, endpoint=True)
    >>> x2 = np.linspace(0, 10, N, endpoint=False)
    >>> plt.plot(x1, y, 'o')
    >>> plt.plot(x2, y + 0.5, 'o')
    >>> plt.ylim([-0.5, 1])
    >>> plt.show()

    """
    # bail if threshold not met
    if not should_distribute(num):
        return np.linspace(start,stop,num,endpoint,retstep)
    a = ndarray(num)
    step = None
    if endpoint:
        step = (stop-start)/(num-1)
    else:
        step = (stop-start)/num
    buf = a.access()
    if buf is not None:
        lo,hi = a.distribution()
        lo,hi = lo[0],hi[0]
        buf[:] = np.arange(lo,hi)*step+start
        a.release_update()
    #sync()
    if retstep:
        return a,step
    return a

def logspace(start, stop, num=50, endpoint=True, base=10.0):
    """Return numbers spaced evenly on a log scale.

    In linear space, the sequence starts at ``base ** start``
    (`base` to the power of `start`) and ends with ``base ** stop``
    (see `endpoint` below).

    Parameters
    ----------
    start : float
        ``base ** start`` is the starting value of the sequence.
    stop : float
        ``base ** stop`` is the final value of the sequence, unless `endpoint`
        is False.  In that case, ``num + 1`` values are spaced over the
        interval in log-space, of which all but the last (a sequence of
        length ``num``) are returned.
    num : integer, optional
        Number of samples to generate.  Default is 50.
    endpoint : boolean, optional
        If true, `stop` is the last sample. Otherwise, it is not included.
        Default is True.
    base : float, optional
        The base of the log space. The step size between the elements in
        ``ln(samples) / ln(base)`` (or ``log_base(samples)``) is uniform.
        Default is 10.0.

    Returns
    -------
    samples : ndarray
        `num` samples, equally spaced on a log scale.

    See Also
    --------
    arange : Similiar to linspace, with the step size specified instead of the
             number of samples. Note that, when used with a float endpoint, the
             endpoint may or may not be included.
    linspace : Similar to logspace, but with the samples uniformly distributed
               in linear space, instead of log space.

    Notes
    -----
    Logspace is equivalent to the code

    >>> y = linspace(start, stop, num=num, endpoint=endpoint)
    >>> power(base, y)

    Examples
    --------
    >>> np.logspace(2.0, 3.0, num=4)
        array([  100.        ,   215.443469  ,   464.15888336,  1000.        ])
    >>> np.logspace(2.0, 3.0, num=4, endpoint=False)
        array([ 100.        ,  177.827941  ,  316.22776602,  562.34132519])
    >>> np.logspace(2.0, 3.0, num=4, base=2.0)
        array([ 4.        ,  5.0396842 ,  6.34960421,  8.        ])

    Graphical illustration:

    >>> import matplotlib.pyplot as plt
    >>> N = 10
    >>> x1 = np.logspace(0.1, 1, N, endpoint=True)
    >>> x2 = np.logspace(0.1, 1, N, endpoint=False)
    >>> y = np.zeros(N)
    >>> plt.plot(x1, y, 'o')
    >>> plt.plot(x2, y + 0.5, 'o')
    >>> plt.ylim([-0.5, 1])
    >>> plt.show()
    
    """
    # bail if threshold not met
    if not should_distribute(num):
        return np.logspace(start,stop,num,endpoint,base)
    a = ndarray(num)
    step = None
    if endpoint:
        step = (stop-start)/(num-1)
    else:
        step = (stop-start)/num
    buf = a.access()
    if buf is not None:
        lo,hi = a.distribution()
        lo,hi = lo[0],hi[0]
        buf[:] = base**(np.arange(lo,hi)*step+start)
        a.release_update()
    #sync()
    return a

def dot(a, b, out=None):
    """dot(a, b)

    Dot product of two arrays.

    For 2-D arrays it is equivalent to matrix multiplication, and for 1-D
    arrays to inner product of vectors (without complex conjugation). For
    N dimensions it is a sum product over the last axis of `a` and
    the second-to-last of `b`::

        dot(a, b)[i,j,k,m] = sum(a[i,j,:] * b[k,:,m])

    Parameters
    ----------
    a : array_like
        First argument.
    b : array_like
        Second argument.

    Returns
    -------
    output : ndarray
        Returns the dot product of `a` and `b`.  If `a` and `b` are both
        scalars or both 1-D arrays then a scalar is returned; otherwise
        an array is returned.

    Raises
    ------
    ValueError
        If the last dimension of `a` is not the same size as
        the second-to-last dimension of `b`.

    See Also
    --------
    vdot : Complex-conjugating dot product.
    tensordot : Sum products over arbitrary axes.

    Examples
    --------
    >>> np.dot(3, 4)
    12

    Neither argument is complex-conjugated:

    >>> np.dot([2j, 3j], [2j, 3j])
    (-13+0j)

    For 2-D arrays it's the matrix product:

    >>> a = [[1, 0], [0, 1]]
    >>> b = [[4, 1], [2, 2]]
    >>> np.dot(a, b)
    array([[4, 1],
           [2, 2]])

    >>> a = np.arange(3*4*5*6).reshape((3,4,5,6))
    >>> b = np.arange(3*4*5*6)[::-1].reshape((5,4,6,3))
    >>> np.dot(a, b)[2,3,2,1,2,2]
    499128
    >>> sum(a[2,3,2,:] * b[1,2,:,2])
    499128

    """
    a = asarray(a)
    b = asarray(b)
    if not (is_distributed(a) or is_distributed(b)):
        # numpy pass through
        return np.dot(a,b)
    # working with flatiter instances can be expensive, try this opt
    if (isinstance(a,flatiter)
            and isinstance(b,flatiter)
            and a._base is b._base):
        return (a._base * b._base).sum()
    if ((isinstance(a,flatiter) or a.ndim == 1)
            and (isinstance(b,flatiter) or b.ndim == 1)):
        if len(a) != len(b):
            raise ValueError, "objects are not aligned"
        tmp = multiply(a,b)
        ndtmp = tmp.access()
        local_sum = None
        if ndtmp is None:
            local_sum = np.add.reduce(np.asarray([0], dtype=tmp.dtype))
        else:
            local_sum = np.add.reduce(ndtmp)
        return ga.gop_add(local_sum)
    elif a.ndim == 2 and b.ndim == 2:
        if a.shape[1] != b.shape[0]:
            raise ValueError, "objects are not aligned"
        # use GA gemm if certain conditions apply
        valid_types = [np.dtype(np.float32),
                np.dtype(np.float64),
                np.dtype(np.complex64),
                np.dtype(np.complex128)]
        if (a.base is None and b.base is None
                and a.dtype == b.dtype and a.dtype in valid_types):
            out = zeros((a.shape[0],b.shape[1]), a.dtype)
            ga.gemm(False, False, a.shape[0], b.shape[1], b.shape[0],
                    1, a.handle, b.handle, 1, out.handle)
            return out
        else:
            raise NotImplementedError
    elif isinstance(a,(ndarray,flatiter)) and isinstance(b,(ndarray,flatiter)):
        if a.shape[1] != b.shape[0]:
            raise ValueError, "objects are not aligned"
        raise NotImplementedError, "arbitrary dot"
    else:
        # assume we have a scalar somewhere, so just multiply
        return multiply(a,b)

def diag(v, k=0):
    """Extract a diagonal or construct a diagonal array.

    Parameters
    ----------
    v : array_like
        If `v` is a 2-D array, return a copy of its `k`-th diagonal.
        If `v` is a 1-D array, return a 2-D array with `v` on the `k`-th
        diagonal.
    k : int, optional
        Diagonal in question. The default is 0. Use `k>0` for diagonals
        above the main diagonal, and `k<0` for diagonals below the main
        diagonal.

    Returns
    -------
    out : ndarray
        The extracted diagonal or constructed diagonal array.

    See Also
    --------
    diagonal : Return specified diagonals.
    diagflat : Create a 2-D array with the flattened input as a diagonal.
    trace : Sum along diagonals.
    triu : Upper triangle of an array.
    tril : Lower triange of an array.

    Examples
    --------
    >>> x = np.arange(9).reshape((3,3))
    >>> x
    array([[0, 1, 2],
           [3, 4, 5],
           [6, 7, 8]])

    >>> np.diag(x)
    array([0, 4, 8])
    >>> np.diag(x, k=1)
    array([1, 5])
    >>> np.diag(x, k=-1)
    array([3, 7])

    >>> np.diag(np.diag(x))
    array([[0, 0, 0],
           [0, 4, 0],
           [0, 0, 8]])

    """
    v = asarray(v)
    if isinstance(v, ndarray):
        raise NotImplementedError, "TODO"
        # the following isn't right.
        # We want to scatter the values from the given diagonal into a brand
        # new distributed array, but how to compute the indices for the
        # scatter operation?  Or should we "access" the newly created array
        # and "gather" values from the given diagonal?
        #if v.ndim == 1:
        #    k_fabs = math.fabs(k)
        #    N = k_fabs + len(v)
        #    a = zeros((N,N), dtype=v.dtype)
        #    ndv = v.access()
        #    if ndv is not None:
        #        lo,hi = v.distribution()
        #        count = hi[0]-lo[0]
        #        indices = np.ndarray(count*2,dtype=int)
        #        if k >= 0:
        #            indices[0::2] = np.arange(count)+lo[0]
        #            indices[1::2] = np.arange(count)+lo[0]+k
        #        else:
        #            indices[0::2] = np.arange(count)+lo[0]+k_fabs
        #            indices[1::2] = np.arange(count)+lo[0]
        #        a.scatter(
        #    return a
        #elif v.ndim == 2:
        #    pass
        #else:
        #    raise ValueError, "Input must be 1- or 2-d."
    else:
        return np.diag(v,k)

def clip(a, a_min, a_max, out=None):
    """Clip (limit) the values in an array.

    Given an interval, values outside the interval are clipped to
    the interval edges.  For example, if an interval of ``[0, 1]``
    is specified, values smaller than 0 become 0, and values larger
    than 1 become 1.

    Parameters
    ----------
    a : array_like
        Array containing elements to clip.
    a_min : scalar or array_like
        Minimum value.
    a_max : scalar or array_like
        Maximum value.  If `a_min` or `a_max` are array_like, then they will
        be broadcasted to the shape of `a`.
    out : ndarray, optional
        The results will be placed in this array. It may be the input
        array for in-place clipping.  `out` must be of the right shape
        to hold the output.  Its type is preserved.

    Returns
    -------
    clipped_array : ndarray
        An array with the elements of `a`, but where values
        < `a_min` are replaced with `a_min`, and those > `a_max`
        with `a_max`.

    See Also
    --------
    numpy.doc.ufuncs : Section "Output arguments"

    Examples
    --------
    >>> a = np.arange(10)
    >>> np.clip(a, 1, 8)
    array([1, 1, 2, 3, 4, 5, 6, 7, 8, 8])
    >>> a
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> np.clip(a, 3, 6, out=a)
    array([3, 3, 3, 3, 4, 5, 6, 6, 6, 6])
    >>> a = np.arange(10)
    >>> a
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> np.clip(a, [3,4,1,1,1,4,4,4,4,4], 8)
    array([3, 4, 2, 3, 4, 5, 6, 7, 8, 8])

    """
    # just in case
    a = asarray(a)
    a_min = asarray(a_min)
    a_max = asarray(a_max)
    if not (is_distributed(a)
            or is_distributed(a_min)
            or is_distributed(a_max)
            or is_distributed(out)):
        # no ndarray instances used, pass through immediately to numpy
        return np.clip(a, a_min, a_max, out)
    a_shape = get_shape(a)
    a_min_shape = get_shape(a_min)
    a_max_shape = get_shape(a_max)
    if out is None:
        out = ndarray(a_shape, get_dtype(a))
    # sanity checks
    if not is_array(out):
        raise TypeError, "output must be an array"
    if out.shape != a.shape:
        raise ValueError, ("clip: Output array must have thesame shape as "
                "the input.")
    # Now figure out what to do...
    if isinstance(out, ndarray):
        sync()
        # get out as an np.ndarray first
        npout = out.access()
        if npout is not None: # this proc owns data
            # get matching and compatible portions of input arrays
            # broadcasting rules (may) apply
            if a is out:
                npa,release_a = npout,False
            else:
                npa,release_a = _npin_piece_based_on_out(a,out,a_shape)
            if a_min is out:
                npa_min,release_a_min = npout,False
            elif a_min is a:
                npa_min,release_a_min = npa,False
            else:
                npa_min,release_a_min = _npin_piece_based_on_out(
                        a_min,out,a_min_shape)
            if a_max is out:
                npa_max,release_a_max = npout,False
            elif a_max is a:
                npa_max,release_a_max = npa,False
            elif a_max is a_min:
                npa_max,release_a_max = npa_min,False
            else:
                npa_max,release_a_max = _npin_piece_based_on_out(
                        a_max,out,a_max_shape)
            np.clip(npa, npa_min, npa_max, npout)
            if release_a:
                a.release()
            if release_a_min:
                a_min.release()
            if release_a_max:
                a_max.release()
            out.release_update()
        #sync()
    elif isinstance(out, flatiter):
        raise NotImplementedError, "flatiter version of clip"
        #sync()
        ## first op: first and second and out are same object
        #if first is second is out:
        #    self._binary_call(out.base,out.base,out.base,*args,**kwargs)
        #    return out.copy()
        #else:
        #    npout = out.access()
        #    if npout is not None: # this proc 'owns' data
        #        if is_distributed(first):
        #            npfirst = first.get(out._range)
        #        else:
        #            npfirst = first[out._range]
        #        if second is first:
        #            npsecond = npfirst
        #        elif is_distributed(second):
        #            npsecond = second.get(out._range)
        #        else:
        #            npsecond = second[out._range]
        #        self.func(npfirst, npsecond, npout, *args, **kwargs)
        #        out.release_update()
        #sync()
    else:
        sync()
        # out is not distributed
        nda = a
        if is_distributed(a):
            nda = a.allget()
        nda_min = a_min
        if a is a_min:
            nda_min = nda
        elif is_distributed(a_min):
            nda_min = a_min.allget()
        nda_max = a_max
        if a is a_max:
            nda_max = nda
        elif a_max is a_min:
            nda_max = nda_min
        elif is_distributed(a_max):
            nda_max = a_max.allget()
        np.clip(a, nda_min, nda_max, out)
        #sync() # I don't think we need this one
    return out

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
    orig_shape = [dim for dim in dimensions]
    shape = [len(orig_shape)] + orig_shape
    if should_distribute(shape):
        a = zeros(shape, dtype=dtype)
        buf = a.access()
        if buf is not None:
            lo,hi = a.distribution()
            lohi_shape = hi-lo
            for i in range(lo[0],hi[0]):
                vec = np.arange(lohi_shape[i+1])+lo[i+1]
                vec_mod = [None]*len(orig_shape)
                vec_mod[i] = slice(None,None,None)
                values = vec[vec_mod]
                buf[i-lo[0]][:] = values
            a.release_update()
    else:
        return np.indices(orig_shape,dtype)
    return a

def shape(a):
    """Return the shape of an array.

    Parameters
    ----------
    a : array_like
        Input array.

    Returns
    -------
    shape : tuple of ints
        The elements of the shape tuple give the lengths of the
        corresponding array dimensions.

    See Also
    --------
    alen
    ndarray.shape : Equivalent array method.

    Examples
    --------
    >>> np.shape(np.eye(3))
    (3, 3)
    >>> np.shape([[1, 2]])
    (1, 2)
    >>> np.shape([0])
    (1,)
    >>> np.shape(0)
    ()

    >>> a = np.array([(1, 2), (3, 4)], dtype=[('x', 'i4'), ('y', 'i4')])
    >>> np.shape(a)
    (2,)
    >>> a.shape
    (2,)
    
    """
    return asarray(a).shape

def diagonal(a, offset=0, axis1=0, axis2=1):
    """Return specified diagonals.

    If `a` is 2-D, returns the diagonal of `a` with the given offset,
    i.e., the collection of elements of the form ``a[i, i+offset]``.  If
    `a` has more than two dimensions, then the axes specified by `axis1`
    and `axis2` are used to determine the 2-D sub-array whose diagonal is
    returned.  The shape of the resulting array can be determined by
    removing `axis1` and `axis2` and appending an index to the right equal
    to the size of the resulting diagonals.

    Parameters
    ----------
    a : array_like
        Array from which the diagonals are taken.
    offset : int, optional
        Offset of the diagonal from the main diagonal.  Can be positive or
        negative.  Defaults to main diagonal (0).
    axis1 : int, optional
        Axis to be used as the first axis of the 2-D sub-arrays from which
        the diagonals should be taken.  Defaults to first axis (0).
    axis2 : int, optional
        Axis to be used as the second axis of the 2-D sub-arrays from
        which the diagonals should be taken. Defaults to second axis (1).

    Returns
    -------
    array_of_diagonals : ndarray
        If `a` is 2-D, a 1-D array containing the diagonal is returned.
        If the dimension of `a` is larger, then an array of diagonals is
        returned, "packed" from left-most dimension to right-most (e.g.,
        if `a` is 3-D, then the diagonals are "packed" along rows).

    Raises
    ------
    ValueError
        If the dimension of `a` is less than 2.

    See Also
    --------
    diag : MATLAB work-a-like for 1-D and 2-D arrays.
    diagflat : Create diagonal arrays.
    trace : Sum along diagonals.

    Examples
    --------
    >>> a = np.arange(4).reshape(2,2)
    >>> a
    array([[0, 1],
           [2, 3]])
    >>> a.diagonal()
    array([0, 3])
    >>> a.diagonal(1)
    array([1])

    A 3-D example:

    >>> a = np.arange(8).reshape(2,2,2); a
    array([[[0, 1],
            [2, 3]],
           [[4, 5],
            [6, 7]]])
    >>> a.diagonal(0, # Main diagonals of two arrays created by skipping
    ...            0, # across the outer(left)-most axis last and
    ...            1) # the "middle" (row) axis first.
    array([[0, 6],
           [1, 7]])

    The sub-arrays whose main diagonals we just obtained; note that each
    corresponds to fixing the right-most (column) axis, and that the
    diagonals are "packed" in rows.

    >>> a[:,:,0] # main diagonal is [0 6]
    array([[0, 2],
           [4, 6]])
    >>> a[:,:,1] # main diagonal is [1 7]
    array([[1, 3],
           [5, 7]])
    
    """
    return asarray(a).diagonal(offset, axis1, axis2)
