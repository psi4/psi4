# cython: profile=True
cimport mpi4py.MPI as MPI
from mpi4py.mpi_c cimport *
import mpi4py.MPI as MPI

from ga4py import ga
import util

import numpy as np
cimport numpy as np
import os

cpdef int me():
    return ga.pgroup_nodeid(ga.pgroup_get_default())

cpdef int nproc():
    return ga.pgroup_nnodes(ga.pgroup_get_default())

# at what point do we distribute arrays versus leaving as np.ndarray?
cdef int SIZE_THRESHOLD = 1
cpdef int get_size_threshold():
    global SIZE_THRESHOLD
    return SIZE_THRESHOLD
cpdef set_size_threshold(int threshold):
    global SIZE_THRESHOLD
    SIZE_THRESHOLD = threshold
cpdef bint should_distribute(shape):
    the_shape = listify(shape)
    if len(the_shape) == 0:
        return False
    return np.multiply.reduce(the_shape) >= get_size_threshold()
cpdef bint is_distributed(thing):
    return isinstance(thing, (ndarray,flatiter))
cpdef bint is_array(thing):
    return isinstance(thing, (ndarray,flatiter,np.ndarray,np.flatiter))
cpdef get_shape(thing):
    try:
        return thing.shape # an ndarray
    except AttributeError:
        return (len(thing),) # a flatiter
cpdef get_dtype(thing):
    try:
        return thing.dtype # an ndarray
    except AttributeError:
        return thing.base.dtype # a flatiter

cpdef list listify(thing):
    try:
        return list(thing)
    except:
        return [thing]

cpdef tuple tuplify(thing):
    try:
        return tuple(thing)
    except:
        return (thing,)

gatypes = {
np.dtype(np.int8):       ga.C_CHAR,
np.dtype(np.int32):      ga.C_INT,
np.dtype(np.int64):      ga.C_LONG,
np.dtype(np.float32):    ga.C_FLOAT,
np.dtype(np.float64):    ga.C_DBL,
np.dtype(np.complex64):  ga.C_SCPL,
np.dtype(np.complex128): ga.C_DCPL,
}
# numpy doesn't always have these types depending on the system
cdef bint float128_in_np = ('float128' in dir(np))
cdef bint complex256_in_np = ('complex256' in dir(np))
if float128_in_np:
    gatypes[np.dtype(np.float128)] = ga.C_LDBL
if complex256_in_np:
    gatypes[np.dtype(np.complex256)] = ga.C_LDCPL

cdef bint _mask_sync = False
cpdef inline mask_sync():
    global _mask_sync
    _mask_sync = True

cpdef inline sync():
    #print "syncing over group %s" % ga.pgroup_get_default()
    #ga.pgroup_sync(ga.pgroup_get_default())
    if not _mask_sync:
        ga.sync() # internally it checks for the default group

class flagsobj(object):
    def __init__(self):
        self._c = True
        self._f = False
        self._o = True
        self._w = True
        self._a = True
        self._u = False
    def _get_c(self):
        return self._c
    c_contiguous = property(_get_c)
    def _get_f(self):
        return self._f
    f_contiguous = property(_get_f)
    def _get_o(self):
        return self._o
    owndata = property(_get_o)
    def _get_w(self):
        return self._w
    writeable = property(_get_w)
    def _get_a(self):
        return self._a
    aligned = property(_get_a)
    def _get_u(self):
        return self._u
    updateifcopy = property(_get_u)
    def __getitem__(self, item):
        if isinstance(item, str):
            if item == "C" or item == "C_CONTIGUOUS":
                return self._c
            if item == "F" or item == "F_CONTIGUOUS":
                return self._f
            if item == "O" or item == "OWNDATA":
                return self._o
            if item == "W" or item == "WRITEABLE":
                return self._w
            if item == "A" or item == "ALIGNED":
                return self._a
            if item == "U" or item == "UPDATEIFCOPY":
                return self._u
        raise KeyError, "Unknown flag"
    def __repr__(self):
        return """  C_CONTIGUOUS : %s
  F_CONTIGUOUS : %s
  OWNDATA : %s
  WRITEABLE : %s
  ALIGNED : %s
  UPDATEIFCOPY : %s""" % (self._c, self._f, self._o,
                          self._w, self._a, self._u)

class GlobalArrayCache(object):
    """When ndarray instances are removed, their GA handles are preserved.

    This class abstracts away various caching schemes in use and proides a
    consistent interface. The scheme that inspired this cache was to preserve
    the last three arrays with the same shape and type. Using the cache avoids
    many create/destroy cycles for GAs which occur as part of temporary array
    creation during numpy codes.

    """
    def __init__(self, level=3):
        self.cache = {}
        self.level = level

    def __del__(self):
        for value in self.cache.values():
            for g_a in value:
                ga.destroy(g_a)

    def __contains__(self, item):
        return (item in self.cache and self.cache[item])

    def __getitem__(self, item):
        if item in self.cache and self.cache[item]:
            return self.cache[item].pop()
        raise KeyError, item

    def __setitem__(self, item, value):
        if item in self.cache:
            self.cache[item].append(value)
        else:
            self.cache[item] = [value]

    def count(self, item):
        if item in self.cache:
            return len(self.cache[item])
        return 0

    def empty(self, item):
        return self.count(item) == 0

    def full(self, item):
        return self.count(item) == self.level

    def pop(self, item):
        return self[item]

    def size(self):
        """Return the size of the cache in bytes."""
        
ga_cache = GlobalArrayCache()

class ndarray(object):
    """ndarray(shape, dtype=float, buffer=None, offset=0,
            strides=None, order=None)

    An array object represents a multidimensional, homogeneous array
    of fixed-size items.  An associated data-type object describes the
    format of each element in the array (its byte-order, how many bytes it
    occupies in memory, whether it is an integer, a floating point number,
    or something else, etc.)

    Arrays should be constructed using `array`, `zeros` or `empty` (refer
    to the See Also section below).  The parameters given here refer to
    a low-level method (`ndarray(...)`) for instantiating an array.

    For more information, refer to the `numpy` module and examine the
    the methods and attributes of an array.

    Parameters
    ----------
    (for the __new__ method; see Notes below)

    shape : tuple of ints
        Shape of created array.
    dtype : data-type, optional
        Any object that can be interpreted as a numpy data type.
    buffer : object exposing buffer interface, optional
        Used to fill the array with data.
    offset : int, optional
        Offset of array data in buffer.
    strides : tuple of ints, optional
        Strides of data in memory.
    order : {'C', 'F'}, optional
        Row-major or column-major order.

    Parameters added for Global Arrays
    ----------------------------------
    base : ndarray
        Should be a "gainarray".  Used during view creation so that a new
        Global Array is not created and other attributes from base are copied
        into the view ndarray.

    Attributes
    ----------
    T : ndarray
        Transpose of the array.
    data : buffer
        The array's elements, in memory.
    dtype : dtype object
        Describes the format of the elements in the array.
    flags : dict
        Dictionary containing information related to memory use, e.g.,
        'C_CONTIGUOUS', 'OWNDATA', 'WRITEABLE', etc.
    flat : numpy.flatiter object
        Flattened version of the array as an iterator.  The iterator
        allows assignments, e.g., ``x.flat = 3`` (See `ndarray.flat` for
        assignment examples; TODO).
    imag : ndarray
        Imaginary part of the array.
    real : ndarray
        Real part of the array.
    size : int
        Number of elements in the array.
    itemsize : int
        The memory use of each array element in bytes.
    nbytes : int
        The total number of bytes required to store the array data,
        i.e., ``itemsize * size``.
    ndim : int
        The array's number of dimensions.
    shape : tuple of ints
        Shape of the array.
    strides : tuple of ints
        The step-size required to move from one element to the next in
        memory. For example, a contiguous ``(3, 4)`` array of type
        ``int16`` in C-order has strides ``(8, 2)``.  This implies that
        to move from element to element in memory requires jumps of 2 bytes.
        To move from row-to-row, one needs to jump 8 bytes at a time
        (``2 * 4``).
    ctypes : ctypes object
        Class containing properties of the array needed for interaction
        with ctypes.
    base : ndarray
        If the array is a view into another array, that array is its `base`
        (unless that array is also a view).  The `base` array is where the
        array data is actually stored.

    Attributes added for Global Arrays
    ----------------------------------
    handle : int
        The Global Arrays handle.
    global_slice: tuple of integers and/or slice objects
        Represents the slice to take from this ndarray. global_slice is
        calculated first based on the shape of the array, then as slices are
        taken from it, slice arithmetic is performed. When an ndarray is
        accessed or converted to an ndarray, the global_slice is used to turn
        the ndarray into its correct shape/strides before returning to the
        caller.
    _is_real : bool
        Whether this is a 'real' view of a complex ndarray.
    _is_imag : bool
        Whether this is an 'imag' view of a complex ndarray.
    _T : None, tuple of ints
        None, or a tuple of ints If a transpose has been applied. 'i' in the
        'j'-th place in the tuple means 'a's 'i'-th axis becomes
        'a.transpose()'s 'j'-th axis.
    _T_inv : None, tuple of ints
        The inverse of _T, or how to reverse the current transpose

    See Also
    --------
    array : Construct an array.
    zeros : Create an array, each element of which is zero.
    empty : Create an array, but leave its allocated memory unchanged (i.e.,
            it contains "garbage").
    dtype : Create a data-type.

    Notes
    -----
    There are two modes of creating an array using ``__new__``:

    1. If `buffer` is None, then only `shape`, `dtype`, and `order`
       are used.
    2. If `buffer` is an object exposing the buffer interface, then
       all keywords are interpreted.

    No ``__init__`` method is needed because the array is fully initialized
    after the ``__new__`` method.

    Examples
    --------
    These examples illustrate the low-level `ndarray` constructor.  Refer
    to the `See Also` section above for easier ways of constructing an
    ndarray.

    First mode, `buffer` is None:

    >>> np.ndarray(shape=(2,2), dtype=float, order='F')
    array([[ -1.13698227e+002,   4.25087011e-303],
           [  2.88528414e-306,   3.27025015e-309]])

    Second mode:

    >>> np.ndarray((2,), buffer=np.array([1,2,3]),
    ...            offset=np.int_().itemsize,
    ...            dtype=int) # offset = 1*itemsize, i.e. skip first element
    array([2, 3])

    """
    def __new__(cls, shape, dtype=float, buffer=None, offset=0,
            strides=None, order=None, base=None):
        if base is None and not should_distribute(shape):
            return np.ndarray(shape, dtype, buffer, offset, strides, order)
        return super(ndarray, cls).__new__(cls)

    def __init__(self, shape, dtype=float, buffer=None, offset=0,
            strides=None, order=None, base=None):
        shape = tuplify(shape)
        if order not in [None,'C','F']:
            raise TypeError, "order not understood"
        if order is None:
            order = 'C'
        if order is 'F':
            raise NotImplementedError, "Fortran order not supported"
        self._dtype = np.dtype(dtype)
        self._order = order
        self._base = base
        self._is_real = False
        self._is_imag = False
        if base is None:
            self.global_slice = util.MasterKey(shape)
            self._flags = flagsobj()
            dtype_ = self._dtype
            gatype = None
            if dtype_ in gatypes:
                gatype = gatypes[dtype_]
            else:
                gatype = ga.register_dtype(dtype_)
                gatypes[dtype_] = gatype
            if (shape,dtype_.str) in ga_cache:
                self.handle = ga_cache.pop((shape,dtype_.str))
            else:
                self.handle = ga.create(gatype, shape,
                        pgroup=ga.pgroup_get_default())
            if buffer is not None:
                local = ga.access(self.handle)
                if local is not None:
                    a = None
                    if isinstance(buffer, np.ndarray):
                        buffer.shape = shape
                        a = buffer
                    else:
                        a = np.ndarray(shape, dtype_, buffer, offset,
                                strides, order)
                    local[:] = a[ga.zip(*self.distribution())]
                    self.release_update()
            self._strides = [self.itemsize]
            for size in shape[-1:0:-1]:
                self._strides = [size*self._strides[0]] + self._strides
        elif isinstance(base, ndarray):
            self.global_slice = base.global_slice
            self.handle = base.handle
            self._strides = strides
            self._is_real = base._is_real
            self._is_imag = base._is_imag
            self._flags = base._flags
            self._flags._c = False
            self._flags._o = False
        else:
            # assume base is a g_a handle
            self.handle = base
            self.global_slice = util.MasterKey(shape)
            self._flags = flagsobj()
            self._strides = [self.itemsize]
            for size in shape[-1:0:-1]:
                self._strides = [size*self._strides[0]] + self._strides

    def __del__(self):
        if self._base is None:
            if ga.initialized():
                shape = self.shape
                if ga_cache.full((shape,self.dtype.str)):
                    ga.destroy(self.handle)
                else:
                    ga_cache[(shape,self.dtype.str)] = self.handle

    ################################################################
    ### ndarray methods added for Global Arrays
    ################################################################
    def distribution(self):
        """Return the bounds of the distribution.

        Always returns the lo/hi distribution based on the original shape of
        the array.

        This operation is local.

        """
        return ga.distribution(self.handle)

    def owns(self):
        """Return True if this process owns some of the data.

        This operation is local.

        """
        lo,hi = self.distribution()
        return np.all(hi>=0)

    def access(self, global_slice=None):
        """Access the local array. Return None if no data is owned.
        
        This operation is local.
        
        """
        # we access the entire block for this process
        a = ga.access(self.handle)
        if a is None:
            return None # bail; this process doesn't own any data
        # get the original GA distribution
        lo,hi = ga.distribution(self.handle)
        # caller may pass in a different global_slice
        if global_slice is None:
            global_slice = self.global_slice
        # transpose, if needed
        b = a.transpose(global_slice.lohi_T())
        # slice away any fixed dimensions
        try:
            fixed_slice = global_slice.access_key(lo,hi)
        except IndexError:
            self.release()
            return None # bail; this process doesn't own relevant piece
        c = b[fixed_slice]
        # at this point the shape of c should match the global_slice shape
        # except for any None/newaxis added
        d = c[global_slice.None_key()]
        return d

    def get(self, key=None):
        """Similar to the __getitem__ built-in, but one-sided (not collective.)

        We sometimes want the semantics of "slicing" an ndarray and then
        immediately calling ga.get() to fetch the result.  We can't use the
        __getitem__ built-in because it is a collective operation.  For
        example, during a ufunc we need to get() the corresponding pieces of
        the arrays and that is where this function is handy.

        This operation is one-sided.

        """
        # caller may modify the global_slice
        global_slice = self.global_slice
        if key is not None:
            global_slice = self.global_slice[key]
        # determine the lo/hi/skip for the ga.strided_get()
        _lo,_hi,_skip,adjust,need_strided = global_slice.get_lo_hi_skip_adjust()
        if need_strided:
            ret = ga.strided_get(self.handle, _lo, _hi, _skip)
        else:
            ret = ga.get(self.handle, _lo, _hi)
        # pull out real or imag parts, if needed
        if self._is_real:
            ret = ret.real
        elif self._is_imag:
            ret = ret.imag
        # apply the adjustment (which only reverses some dimensions, if needed)
        ret = ret[adjust]
        # transpose the result, if needed
        ret = ret.transpose(global_slice.lohi_T())
        # reshape the array based on global_slice
        # this will add any newaxis and remove any fixed dimensions
        ret.shape = global_slice.shape
        return ret

    def allget(self, key=None):
        """Like get(), but when all processes need the same piece.
        
        This operation is collective.
        
        """
        # TODO it's not clear whether this approach is better than having all
        # P processors ga.get() the same piece.
        if not me():
            result = self.get(key)
            return comm().bcast(result)
        else:
            return comm().bcast()

    def gather(self, i):
        """Use ga.gather() but modify the indicies `i' based on global_slice.

        Parameters
        ----------
        i : 1D or 2D array of indices
            see ga.gather() documentation)

        """
        # we need to modify the indices based on the global_slice
        # the underlying GA instance might not have the same shape, may be a
        # view, etc
        # make sure we're always working with a 2D index array
        i = np.asarray(i, dtype=np.int64)
        if i.ndim == 1:
            i.shape = (-1,self.ndim)
        step = []
        start = []
        for gs in self.global_slice:
            if isinstance(gs, util.RangeKey):
                step.append(gs.step)
                start.append(gs.start)
            elif isinstance(gs, util.NoneKey):
                step.append(1)
                start.append(0)
            elif isinstance(gs, util.FixedKey):
                pass
            else:
                raise TypeError, "unhandled piece of global_slice"
        # modify new index array based on global_slice
        new_i = (i*np.asarray(step,dtype=np.int64)[None,:]
                + np.asarray(start, dtype=np.int64)[None,:])
        # get rid of index columns which don't refer to an actual GA dimension
        # add index columns which were missing, single-valued GA dimensions
        columns = np.hsplit(new_i,self.ndim) # returns a list
        column_iter = iter(columns)
        new_columns = [None]*self.global_slice.get_original_ndim()
        for fixed in self.global_slice.fixed:
            new_column = np.empty((len(new_i),1),dtype=np.int64)
            new_column[:] = fixed.value
            new_columns[fixed.origin] = new_column
        for col,gs in zip(columns,self.global_slice.data):
            if isinstance(gs, util.RangeKey):
                # no change to column
                new_columns[gs.origin] = col
            elif isinstance(gs, util.NoneKey):
                # skip column
                pass
            elif isinstance(gs, util.FixedKey):
                raise TypeError, "FixedKey found in MasterKey"
            else:
                raise TypeError, "unhandled piece of global_slice"
        new_i = np.hstack(new_columns)
        return ga.gather(self.handle, new_i)

    def release(self):
        ga.release(self.handle)

    def release_update(self):
        ga.release_update(self.handle)

    ################################################################
    ### ndarray properties
    ################################################################

    def _get_T(self):
        return self.transpose()
    T = property(_get_T)

    def _get_data(self):
        a = self.access()
        if a is not None:
            return a.data
        return None
    data = property(_get_data)

    def _get_dtype(self):
        return self._dtype
    dtype = property(_get_dtype)

    def _get_flags(self):
        return self._flags
    flags = property(_get_flags)

    def _get_flat(self):
        return flatiter(self)
    def _set_flat(self, value):
        a = flatiter(self)
        a[:] = value
    flat = property(_get_flat,_set_flat)

    def _get_imag(self):
        if self._dtype.kind != 'c':
            return zeros(self.shape, self.dtype)
        else:
            ret = self[:]
            ret._is_imag = True
            ret._dtype = np.dtype("float%s" % (self._dtype.itemsize/2*8))
            return ret
    def _set_imag(self, value):
        if self._dtype.kind != 'c':
            raise TypeError, "array does not have imaginary part to set"
        else:
            self._get_imag()[:] = value
    imag = property(_get_imag,_set_imag)

    def _get_real(self):
        if self._dtype.kind != 'c':
            return self
        else:
            ret = self[:]
            ret._is_real = True
            ret._dtype = np.dtype("float%s" % (self._dtype.itemsize/2*8))
            return ret
    def _set_real(self, value):
        self._get_real()[:] = value
    real = property(_get_real,_set_real)

    def _get_size(self):
        return reduce(lambda x,y: x*y, self.shape, 1)
    size = property(_get_size)

    def _get_itemsize(self):
        return self._dtype.itemsize
    itemsize = property(_get_itemsize)

    def _get_nbytes(self):
        return self.itemsize * self.size
    nbytes = property(_get_nbytes)

    def _get_ndim(self):
        return len(self.shape)
    ndim = property(_get_ndim)

    def _get_shape(self):
        return self.global_slice.shape
    def _set_shape(self, value):
        if self.base is None:
            # we can change the underlying GA
            # we do this by creating a new ndarray and copying properties
            a = self.reshape(value)
            ga.destroy(self.handle)
            self.handle = a.handle
            self.global_slice = a.global_slice
            self._strides = a.strides
            self._flags = a._flags
        else:
            raise AttributeError, "incompatible shape for a non-contiguous array"
    shape = property(_get_shape, _set_shape)

    def _get_strides(self):
        strides = [self.itemsize]
        for size in self.shape[-1:0:-1]:
            strides = [size*strides[0]] + strides
        return strides
    strides = property(_get_strides)

    def _get_ctypes(self):
        raise NotImplementedError, "TODO"
    ctypes = property(_get_ctypes)

    def _get_base(self):
        return self._base
    base = property(_get_base)

    ################################################################
    ### ndarray methods
    ################################################################
    def all(self, axis=None, out=None):
        """Returns True if all elements evaluate to True.

        Refer to `numpy.all` for full documentation.

        See Also
        --------
        numpy.all : equivalent function

        """
        raise NotImplementedError

    def any(self, axis=None, out=None):
        """    Returns True if any of the elements of `a` evaluate to True.

        Refer to `numpy.any` for full documentation.

        See Also
        --------
        numpy.any : equivalent function

        """
        raise NotImplementedError

    def argmax(self, axis=None, out=None):
        """Return indices of the maximum values along the given axis.

        Refer to `numpy.argmax` for full documentation.

        See Also
        --------
        numpy.argmax : equivalent function

        """
        raise NotImplementedError

    def argmin(self, axis=None, out=None):
        """Return indices of the minimum values along the given axis of `a`.

        Refer to `numpy.argmin` for detailed documentation.

        See Also
        --------
        numpy.argmin : equivalent function

        """
        raise NotImplementedError

    def argsort(self, axis=-1, kind='quicksort', order=None):
        """Returns the indices that would sort this array.

        Refer to `numpy.argsort` for full documentation.

        See Also
        --------
        numpy.argsort : equivalent function

        """
        raise NotImplementedError

    def astype(self, t):
        """Copy of the array, cast to a specified type.

        Parameters
        ----------
        t : string or dtype
            Typecode or data-type to which the array is cast.

        Examples
        --------
        >>> x = np.array([1, 2, 2.5])
        >>> x
        array([ 1. ,  2. ,  2.5])

        >>> x.astype(int)
        array([1, 2, 2])

        """
        # TODO we can optimize a copy if new and old ndarray instances align
        the_copy = ndarray(self.shape, dtype=t)
        if should_distribute(the_copy.size):
            local = the_copy.access()
            if local is not None:
                lo,hi = the_copy.distribution()
                local[:] = self.get(ga.zip(lo,hi))
                the_copy.release_update()
        else:
            # case where the copy is not distributed but the original was
            the_copy[:] = self.allget()
        return the_copy

    def byteswap(self, inplace=False):
        """Swap the bytes of the array elements

        Toggle between low-endian and big-endian data representation by
        returning a byteswapped array, optionally swapped in-place.

        Parameters
        ----------
        inplace: bool, optional
            If ``True``, swap bytes in-place, default is ``False``.

        Returns
        -------
        out: ndarray
            The byteswapped array. If `inplace` is ``True``, this is
            a view to self.

        Examples
        --------
        >>> A = np.array([1, 256, 8755], dtype=np.int16)
        >>> map(hex, A)
        ['0x1', '0x100', '0x2233']
        >>> A.byteswap(True)
        array([  256,     1, 13090], dtype=int16)
        >>> map(hex, A)
        ['0x100', '0x1', '0x3322']

        Arrays of strings are not swapped

        >>> A = np.array(['ceg', 'fac'])
        >>> A.byteswap()
        array(['ceg', 'fac'],
              dtype='|S3')
            
        """
        raise NotImplementedError

    def choose(self, choices, out=None, mode='raise'):
        """Use an index array to construct a new array from a set of choices.

        Refer to `numpy.choose` for full documentation.

        See Also
        --------
        numpy.choose : equivalent function

        """
        raise NotImplementedError

    def clip(self, a_min, a_max, out=None):
        """Return an array whose values are limited to ``[a_min, a_max]``.

        Refer to `numpy.clip` for full documentation.

        See Also
        --------
        numpy.clip : equivalent function

        """
        raise NotImplementedError

    def compress(self, condition, axis=None, out=None):
        """Return selected slices of this array along given axis.

        Refer to `numpy.compress` for full documentation.

        See Also
        --------
        numpy.compress : equivalent function

        """
        return NotImplementedError

    def conj(self):
        """Complex-conjugate all elements.

        Refer to `numpy.conjugate` for full documentation.

        See Also
        --------
        numpy.conjugate : equivalent function

        """
        raise NotImplementedError

    def conjugate(self):
        """Return the complex conjugate, element-wise.

        Refer to `numpy.conjugate` for full documentation.

        See Also
        --------
        numpy.conjugate : equivalent function

        """
        raise NotImplementedError

    def copy(self, order='C'):
        """Return a copy of the array.

        Parameters
        ----------
        order : {'C', 'F', 'A'}, optional
            By default, the result is stored in C-contiguous (row-major) order in
            memory.  If `order` is `F`, the result has 'Fortran' (column-major)
            order.  If order is 'A' ('Any'), then the result has the same order
            as the input.

        Examples
        --------
        >>> x = np.array([[1,2,3],[4,5,6]], order='F')

        >>> y = x.copy()

        >>> x.fill(0)

        >>> x
        array([[0, 0, 0],
               [0, 0, 0]])

        >>> y
        array([[1, 2, 3],
               [4, 5, 6]])

        >>> y.flags['C_CONTIGUOUS']
        True

        """
        # TODO we can optimize a copy if new and old ndarray instances align
        the_copy = ndarray(self.shape, dtype=self.dtype)
        if should_distribute(the_copy.size):
            local = the_copy.access()
            if local is not None:
                lo,hi = the_copy.distribution()
                local[:] = self.get(ga.zip(lo,hi))
                the_copy.release_update()
        else:
            # case where the copy is not distributed but the original was
            the_copy[:] = self.allget()
        return the_copy

    def cumprod(self, axis=None, dtype=None, out=None):
        """Return the cumulative product of the elements along the given axis.

        Refer to `numpy.cumprod` for full documentation.

        See Also
        --------
        numpy.cumprod : equivalent function

        """
        raise NotImplementedError

    def cumsum(self, axis=None, dtype=None, out=None):
        """Return the cumulative sum of the elements along the given axis.

        Refer to `numpy.cumsum` for full documentation.

        See Also
        --------
        numpy.cumsum : equivalent function

        """
        raise NotImplementedError

    def diagonal(self, offset=0, axis1=0, axis2=1):
        """Return specified diagonals.

        Refer to `numpy.diagonal` for full documentation.

        See Also
        --------
        numpy.diagonal : equivalent function

        """
        if self.ndim < 2:
            raise ValueError, "array.ndim must be >= 2"
        if self.ndim != 2:
            raise NotImplementedError, "diagonal for a.ndim > 2"
        nrow,ncol = self.shape
        if offset >= 0:
            size = min(ncol-offset,nrow)
        else:
            size = min(nrow+offset,ncol)
        out = ndarray(size, dtype=self.dtype)
        if not is_distributed(out):
            raise NotImplementedError, "resulting diagonal is not distributed"
        npout = out.access()
        if npout is not None:
            indices = np.empty((size,2), dtype=np.int64)
            if offset >= 0:
                indices[:,0] = np.arange(size, dtype=np.int64)
                indices[:,1] = np.arange(offset, size+offset, dtype=np.int64)
            else:
                offset = -offset
                indices[:,0] = np.arange(offset, size+offset, dtype=np.int64)
                indices[:,1] = np.arange(size, dtype=np.int64)
            my_range = out.global_slice.bound_by_lohi(*out.distribution())
            my_range = my_range[0].pyobj()
            my_indices = indices[my_range]
            npout[:] = self.gather(my_indices)
        return out

    def dot(self, b, out=None):
        """a.dot(b, out=None)

        Dot product of two arrays.

        Refer to `numpy.dot` for full documentation.

        See Also
        --------
        numpy.dot : equivalent function

        Examples
        --------
        >>> a = np.eye(2)
        >>> b = np.ones((2, 2)) * 2
        >>> a.dot(b)
        array([[ 2.,  2.],
               [ 2.,  2.]])

        This array method can be conveniently chained:

        >>> a.dot(b).dot(b)
        array([[ 8.,  8.],
               [ 8.,  8.]])

        """
        raise NotImplementedError

    def dump(self, file):
        """Dump a pickle of the array to the specified file.

        The array can be read back with pickle.load or numpy.load.

        Parameters
        ----------
        file : str
            A string naming the dump file.

        """
        raise NotImplementedError

    def dumps(self):
        """Returns the pickle of the array as a string.

        pickle.loads or numpy.loads will convert the string back to an array.

        Parameters
        ----------
        None

        """
        raise NotImplementedError

    def fill(self, value):
        """Fill the array with a scalar value.

        Parameters
        ----------
        value : scalar
            All elements of `a` will be assigned this value.

        Examples
        --------
        >>> a = np.array([1, 2])
        >>> a.fill(0)
        >>> a
        array([0, 0])
        >>> a = np.empty(2)
        >>> a.fill(1)
        >>> a
        array([ 1.,  1.])

        """
        a = self.access()
        if a is not None:
            a.fill(value)
            self.release_update()

    def flatten(self, order='C'):
        """Return a copy of the array collapsed into one dimension.

        Parameters
        ----------
        order : {'C', 'F'}, optional
            Whether to flatten in C (row-major) or Fortran (column-major) order.
            The default is 'C'.

        Returns
        -------
        y : ndarray
            A copy of the input array, flattened to one dimension.

        See Also
        --------
        ravel : Return a flattened array.
        flat : A 1-D flat iterator over the array.

        Examples
        --------
        >>> a = np.array([[1,2], [3,4]])
        >>> a.flatten()
        array([1, 2, 3, 4])
        >>> a.flatten('F')
        array([1, 3, 2, 4])

        """
        return self.reshape(self.size)

    def getfield(self, dtype, offset):
        """Returns a field of the given array as a certain type.

        A field is a view of the array data with each itemsize determined
        by the given type and the offset into the current array, i.e. from
        ``offset * dtype.itemsize`` to ``(offset+1) * dtype.itemsize``.

        Parameters
        ----------
        dtype : str
            String denoting the data type of the field.
        offset : int
            Number of `dtype.itemsize`'s to skip before beginning the element view.

        Examples
        --------
        >>> x = np.diag([1.+1.j]*2)
        >>> x
        array([[ 1.+1.j,  0.+0.j],
               [ 0.+0.j,  1.+1.j]])
        >>> x.dtype
        dtype('complex128')

        >>> x.getfield('complex64', 0) # Note how this != x
        array([[ 0.+1.875j,  0.+0.j   ],
               [ 0.+0.j   ,  0.+1.875j]], dtype=complex64)

        >>> x.getfield('complex64',1) # Note how different this is than x
        array([[ 0. +5.87173204e-39j,  0. +0.00000000e+00j],
               [ 0. +0.00000000e+00j,  0. +5.87173204e-39j]], dtype=complex64)

        >>> x.getfield('complex128', 0) # == x
        array([[ 1.+1.j,  0.+0.j],
               [ 0.+0.j,  1.+1.j]])

        If the argument dtype is the same as x.dtype, then offset != 0 raises
        a ValueError:

        >>> x.getfield('complex128', 1)
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
        ValueError: Need 0 <= offset <= 0 for requested type but received offset = 1

        >>> x.getfield('float64', 0)
        array([[ 1.,  0.],
               [ 0.,  1.]])

        >>> x.getfield('float64', 1)
        array([[  1.77658241e-307,   0.00000000e+000],
               [  0.00000000e+000,   1.77658241e-307]])

        """
        raise NotImplementedError

    def item(self, *args):
        """Copy an element of an array to a standard Python scalar and return it.

        Parameters
        ----------
        \*args : Arguments (variable number and type)

            * none: in this case, the method only works for arrays
              with one element (`a.size == 1`), which element is
              copied into a standard Python scalar object and returned.

            * int_type: this argument is interpreted as a flat index into
              the array, specifying which element to copy and return.

            * tuple of int_types: functions as does a single int_type argument,
              except that the argument is interpreted as an nd-index into the
              array.

        Returns
        -------
        z : Standard Python scalar object
            A copy of the specified element of the array as a suitable
            Python scalar

        Notes
        -----
        When the data type of `a` is longdouble or clongdouble, item() returns
        a scalar array object because there is no available Python scalar that
        would not lose information. Void arrays return a buffer object for item(),
        unless fields are defined, in which case a tuple is returned.

        `item` is very similar to a[args], except, instead of an array scalar,
        a standard Python scalar is returned. This can be useful for speeding up
        access to elements of the array and doing arithmetic on elements of the
        array using Python's optimized math.

        Examples
        --------
        >>> x = np.random.randint(9, size=(3, 3))
        >>> x
        array([[3, 1, 7],
               [2, 8, 3],
               [8, 5, 3]])
        >>> x.item(3)
        2
        >>> x.item(7)
        5
        >>> x.item((0, 1))
        1
        >>> x.item((2, 2))
        3

        """
        raise NotImplementedError

    def itemset(self, *args):
        """Insert scalar into an array (scalar is cast to array's dtype, if possible)

        There must be at least 1 argument, and define the last argument
        as *item*.  Then, ``a.itemset(*args)`` is equivalent to but faster
        than ``a[args] = item``.  The item should be a scalar value and `args`
        must select a single item in the array `a`.

        Parameters
        ----------
        \*args : Arguments
            If one argument: a scalar, only used in case `a` is of size 1.
            If two arguments: the last argument is the value to be set
            and must be a scalar, the first argument specifies a single array
            element location. It is either an int or a tuple.

        Notes
        -----
        Compared to indexing syntax, `itemset` provides some speed increase
        for placing a scalar into a particular location in an `ndarray`,
        if you must do this.  However, generally this is discouraged:
        among other problems, it complicates the appearance of the code.
        Also, when using `itemset` (and `item`) inside a loop, be sure
        to assign the methods to a local variable to avoid the attribute
        look-up at each loop iteration.

        Examples
        --------
        >>> x = np.random.randint(9, size=(3, 3))
        >>> x
        array([[3, 1, 7],
               [2, 8, 3],
               [8, 5, 3]])
        >>> x.itemset(4, 0)
        >>> x.itemset((2, 2), 9)
        >>> x
        array([[3, 1, 7],
               [2, 0, 3],
               [8, 5, 9]])

        """
        raise NotImplementedError

    def max(self, axis=None, out=None):
        """Return the maximum along a given axis.

        Refer to `numpy.amax` for full documentation.

        See Also
        --------
        numpy.amax : equivalent function

        """
        if axis is None:
            a = self.access()
            if a is not None:
                result = np.maximum.reduce(a)
                while result.ndim > 0:
                    result = np.maximum.reduce(result)
                result = comm().allgather(result)
                result = np.maximum.reduce(result)
            else:
                # we don't own anything, so get an actual value from the array
                result = self.get([0]*self.ndim)
                result = comm().allgather(result)
                result = np.maximum.reduce(result)
            if out is None:
                return result
            else:
                out[:] = result
                return out
        else:
            return maximum.reduce(self, axis, out)

    def mean(self, axis=None, dtype=None, out=None):
        """Returns the average of the array elements along given axis.

        Refer to `numpy.mean` for full documentation.

        See Also
        --------
        numpy.mean : equivalent function

        """
        raise NotImplementedError

    def min(self, axis=None, out=None):
        """Return the minimum along a given axis.

        Refer to `numpy.amin` for full documentation.

        See Also
        --------
        numpy.amin : equivalent function

        """
        if axis is None:
            a = self.access()
            if a is not None:
                result = np.minimum.reduce(a)
                while result.ndim > 0:
                    result = np.minimum.reduce(result)
                result = comm().allgather(result)
                result = np.minimum.reduce(result)
            else:
                # we don't own anything, so get an actual value from the array
                result = self.get([0]*self.ndim)
                result = comm().allgather(result)
                result = np.minimum.reduce(result)
            if out is None:
                return result
            else:
                out[:] = result
                return out
        else:
            return minimum.reduce(self, axis, out)

    def newbyteorder(self, new_order='S'):
        """Return the array with the same data viewed with a different byte order.

        Equivalent to::

            arr.view(arr.dtype.newbytorder(new_order))

        Changes are also made in all fields and sub-arrays of the array data
        type.



        Parameters
        ----------
        new_order : string, optional
            Byte order to force; a value from the byte order specifications
            above. `new_order` codes can be any of::

             * 'S' - swap dtype from current to opposite endian
             * {'<', 'L'} - little endian
             * {'>', 'B'} - big endian
             * {'=', 'N'} - native order
             * {'|', 'I'} - ignore (no change to byte order)

            The default value ('S') results in swapping the current
            byte order. The code does a case-insensitive check on the first
            letter of `new_order` for the alternatives above.  For example,
            any of 'B' or 'b' or 'biggish' are valid to specify big-endian.


        Returns
        -------
        new_arr : array
            New array object with the dtype reflecting given change to the
            byte order.

        """
        raise NotImplementedError

    def nonzero(self):
        """Return the indices of the elements that are non-zero.

        Refer to `numpy.nonzero` for full documentation.

        See Also
        --------
        numpy.nonzero : equivalent function

        """
        raise NotImplementedError

    def prod(self, axis=None, dtype=None, out=None):
        """Return the product of the array elements over the given axis

        Refer to `numpy.prod` for full documentation.

        See Also
        --------
        numpy.prod : equivalent function

        """
        raise NotImplementedError

    def ptp(self, axis=None, out=None):
        """Peak to peak (maximum - minimum) value along a given axis.

        Refer to `numpy.ptp` for full documentation.

        See Also
        --------
        numpy.ptp : equivalent function

        """
        raise NotImplementedError

    def put(self, indices, values, mode='raise'):
        """Set ``a.flat[n] = values[n]`` for all `n` in indices.

        Refer to `numpy.put` for full documentation.

        See Also
        --------
        numpy.put : equivalent function
        
        """
        raise NotImplementedError

    def ravel(self, order=None):
        """Return a flattened array.

        Refer to `numpy.ravel` for full documentation.

        See Also
        --------
        numpy.ravel : equivalent function

        ndarray.flat : a flat iterator on the array.

        """
        raise NotImplementedError

    def repeat(self, repeats, axis=None):
        """Repeat elements of an array.

        Refer to `numpy.repeat` for full documentation.

        See Also
        --------
        numpy.repeat : equivalent function

        """
        raise NotImplementedError

    def reshape(self, shape, order='C', *args):
        """Returns an array containing the same data with a new shape.

        Refer to `numpy.reshape` for full documentation.

        See Also
        --------
        numpy.reshape : equivalent function

        """
        if type(order) == type(""):
            # assume user specified a shape tuple
            shape = listify(shape)
        else:
            # assume user specified an unpacked shape e.g. reshape(3,4,5)
            shape = listify([shape]+[order]+list(args))
            order = 'C'
        shape = np.asarray(shape, dtype=np.int64)
        count_neg_ones = shape[shape==-1].size
        if count_neg_ones > 1:
            raise ValueError, "can only specify one unknown dimension"
        if count_neg_ones == 1:
            shape_product = np.prod(shape)*(-1)
            if np.size%shape_product != 0:
                raise ValueError, "total size of new array must be unchanged"
            shape[shape==-1] = np.size//shape_product
        # now that shape is established, create new array
        a = ndarray(shape, dtype=self.dtype)
        # based on distribution, gather same indices from self
        nda = a.access()
        if nda is not None:
            lo,hi = a.distribution()
            lohi_shape = hi-lo
            i = np.indices(lohi_shape).reshape(len(lohi_shape),-1).T + lo
            # get the flattened indices
            strides = [1]
            for size in shape[-1:0:-1]:
                strides = [size*strides[0]] + strides
            i_flat = np.sum(i*strides, axis=1)
            # unravel the flattened indices based on the original array
            i_unravel = util.unravel_index(i_flat, self.shape)
            nda.flat = self.gather(i_unravel)
            a.release_update()
        return a

    def resize(self, new_shape, refcheck=True):
        """Change shape and size of array in-place.

        Parameters
        ----------
        new_shape : tuple of ints, or `n` ints
            Shape of resized array.
        refcheck : bool, optional
            If False, reference count will not be checked. Default is True.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If `a` does not own its own data or references or views to it exist,
            and the data memory must be changed.

        SystemError
            If the `order` keyword argument is specified. This behaviour is a
            bug in NumPy.

        See Also
        --------
        resize : Return a new array with the specified shape.

        Notes
        -----
        This reallocates space for the data area if necessary.

        Only contiguous arrays (data elements consecutive in memory) can be
        resized.

        The purpose of the reference count check is to make sure you
        do not use this array as a buffer for another Python object and then
        reallocate the memory. However, reference counts can increase in
        other ways so if you are sure that you have not shared the memory
        for this array with another Python object, then you may safely set
        `refcheck` to False.

        Examples
        --------
        Shrinking an array: array is flattened (in the order that the data are
        stored in memory), resized, and reshaped:

        >>> a = np.array([[0, 1], [2, 3]], order='C')
        >>> a.resize((2, 1))
        >>> a
        array([[0],
               [1]])

        >>> a = np.array([[0, 1], [2, 3]], order='F')
        >>> a.resize((2, 1))
        >>> a
        array([[0],
               [2]])

        Enlarging an array: as above, but missing entries are filled with zeros:

        >>> b = np.array([[0, 1], [2, 3]])
        >>> b.resize(2, 3) # new_shape parameter doesn't have to be a tuple
        >>> b
        array([[0, 1, 2],
               [3, 0, 0]])

        Referencing an array prevents resizing...

        >>> c = a
        >>> a.resize((1, 1))
        Traceback (most recent call last):
        ...
        ValueError: cannot resize an array that has been referenced ...

        Unless `refcheck` is False:

        >>> a.resize((1, 1), refcheck=False)
        >>> a
        array([[0]])
        >>> c
        array([[0]])

        """
        raise NotImplementedError

    def round(self, decimals=0, out=None):
        """Return `a` with each element rounded to the given number of decimals.

        Refer to `numpy.around` for full documentation.

        See Also
        --------
        numpy.around : equivalent function

        """
        raise NotImplementedError

    def searchsorted(self, v, side='left'):
        """Find indices where elements of v should be inserted in a to maintain order.

        For full documentation, see `numpy.searchsorted`

        See Also
        --------
        numpy.searchsorted : equivalent function

        """
        raise NotImplementedError

    def setfield(self, val, dtype, offset=0):
        """Put a value into a specified place in a field defined by a data-type.

        Place `val` into `a`'s field defined by `dtype` and beginning `offset`
        bytes into the field.

        Parameters
        ----------
        val : object
            Value to be placed in field.
        dtype : dtype object
            Data-type of the field in which to place `val`.
        offset : int, optional
            The number of bytes into the field at which to place `val`.

        Returns
        -------
        None

        See Also
        --------
        getfield

        Examples
        --------
        >>> x = np.eye(3)
        >>> x.getfield(np.float64)
        array([[ 1.,  0.,  0.],
               [ 0.,  1.,  0.],
               [ 0.,  0.,  1.]])
        >>> x.setfield(3, np.int32)
        >>> x.getfield(np.int32)
        array([[3, 3, 3],
               [3, 3, 3],
               [3, 3, 3]])
        >>> x
        array([[  1.00000000e+000,   1.48219694e-323,   1.48219694e-323],
               [  1.48219694e-323,   1.00000000e+000,   1.48219694e-323],
               [  1.48219694e-323,   1.48219694e-323,   1.00000000e+000]])
        >>> x.setfield(np.eye(3), np.int32)
        >>> x
        array([[ 1.,  0.,  0.],
               [ 0.,  1.,  0.],
               [ 0.,  0.,  1.]])

        """
        raise NotImplementedError

    def setflags(self, write=None, align=None, uic=None):
        """Set array flags WRITEABLE, ALIGNED, and UPDATEIFCOPY, respectively.

        These Boolean-valued flags affect how numpy interprets the memory
        area used by `a` (see Notes below). The ALIGNED flag can only
        be set to True if the data is actually aligned according to the type.
        The UPDATEIFCOPY flag can never be set to True. The flag WRITEABLE
        can only be set to True if the array owns its own memory, or the
        ultimate owner of the memory exposes a writeable buffer interface,
        or is a string. (The exception for string is made so that unpickling
        can be done without copying memory.)

        Parameters
        ----------
        write : bool, optional
            Describes whether or not `a` can be written to.
        align : bool, optional
            Describes whether or not `a` is aligned properly for its type.
        uic : bool, optional
            Describes whether or not `a` is a copy of another "base" array.

        Notes
        -----
        Array flags provide information about how the memory area used
        for the array is to be interpreted. There are 6 Boolean flags
        in use, only three of which can be changed by the user:
        UPDATEIFCOPY, WRITEABLE, and ALIGNED.

        WRITEABLE (W) the data area can be written to;

        ALIGNED (A) the data and strides are aligned appropriately for the hardware
        (as determined by the compiler);

        UPDATEIFCOPY (U) this array is a copy of some other array (referenced
        by .base). When this array is deallocated, the base array will be
        updated with the contents of this array.

        All flags can be accessed using their first (upper case) letter as well
        as the full name.

        Examples
        --------
        >>> y
        array([[3, 1, 7],
               [2, 0, 0],
               [8, 5, 9]])
        >>> y.flags
          C_CONTIGUOUS : True
          F_CONTIGUOUS : False
          OWNDATA : True
          WRITEABLE : True
          ALIGNED : True
          UPDATEIFCOPY : False
        >>> y.setflags(write=0, align=0)
        >>> y.flags
          C_CONTIGUOUS : True
          F_CONTIGUOUS : False
          OWNDATA : True
          WRITEABLE : False
          ALIGNED : False
          UPDATEIFCOPY : False
        >>> y.setflags(uic=1)
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
        ValueError: cannot set UPDATEIFCOPY flag to True

        """
        raise NotImplementedError

    def sort(self, axis=-1, kind='quicksort', order=None):
        """Sort an array, in-place.

        Parameters
        ----------
        axis : int, optional
            Axis along which to sort. Default is -1, which means sort along the
            last axis.
        kind : {'quicksort', 'mergesort', 'heapsort'}, optional
            Sorting algorithm. Default is 'quicksort'.
        order : list, optional
            When `a` is an array with fields defined, this argument specifies
            which fields to compare first, second, etc.  Not all fields need be
            specified.

        See Also
        --------
        numpy.sort : Return a sorted copy of an array.
        argsort : Indirect sort.
        lexsort : Indirect stable sort on multiple keys.
        searchsorted : Find elements in sorted array.

        Notes
        -----
        See ``sort`` for notes on the different sorting algorithms.

        Examples
        --------
        >>> a = np.array([[1,4], [3,1]])
        >>> a.sort(axis=1)
        >>> a
        array([[1, 4],
               [1, 3]])
        >>> a.sort(axis=0)
        >>> a
        array([[1, 3],
               [1, 4]])

        Use the `order` keyword to specify a field to use when sorting a
        structured array:

        >>> a = np.array([('a', 2), ('c', 1)], dtype=[('x', 'S1'), ('y', int)])
        >>> a.sort(order='y')
        >>> a
        array([('c', 1), ('a', 2)],
              dtype=[('x', '|S1'), ('y', '<i4')])

        """
        raise NotImplementedError

    def squeeze(self):
        """Remove single-dimensional entries from the shape of `a`.

        Refer to `numpy.squeeze` for full documentation.

        See Also
        --------
        numpy.squeeze : equivalent function

        """
        raise NotImplementedError

    def std(self, axis=None, dtype=None, out=None, ddof=0):
        """Remove single-dimensional entries from the shape of `a`.

        Refer to `numpy.squeeze` for full documentation.

        See Also
        --------
        numpy.squeeze : equivalent function

        """
        raise NotImplementedError

    def sum(self, axis=None, dtype=None, out=None):
        """Return the sum of the array elements over the given axis.

        Refer to `numpy.sum` for full documentation.

        See Also
        --------
        numpy.sum : equivalent function

        """
        if axis is None:
            local = self.access()
            if local is not None:
                value = np.sum(local)
                value = comm().allreduce(value, MPI.SUM)
            else:
                value = comm().allreduce(0, MPI.SUM)
            return value
        else:
            raise NotImplementedError

    def swapaxes(self, axis1, axis2):
        """Return a view of the array with `axis1` and `axis2` interchanged.

        Refer to `numpy.swapaxes` for full documentation.

        See Also
        --------
        numpy.swapaxes : equivalent function

        """
        raise NotImplementedError

    def take(self, indices, axis=None, out=None, mode='raise'):
        """Return an array formed from the elements of `a` at the given indices.

        Refer to `numpy.take` for full documentation.

        See Also
        --------
        numpy.take : equivalent function

        """
        raise NotImplementedError

    def tofile(self, fid, sep="", format="%s"):
        """Write array to a file as text or binary (default).

        Data is always written in 'C' order, independent of the order of `a`.
        The data produced by this method can be recovered using the function
        fromfile().

        Parameters
        ----------
        fid : file or str
            An open file object, or a string containing a filename.
        sep : str
            Separator between array items for text output.
            If "" (empty), a binary file is written, equivalent to
            ``file.write(a.tostring())``.
        format : str
            Format string for text file output.
            Each entry in the array is formatted to text by first converting
            it to the closest Python type, and then using "format" % item.

        Notes
        -----
        This is a convenience function for quick storage of array data.
        Information on endianness and precision is lost, so this method is not a
        good choice for files intended to archive data or transport data between
        machines with different endianness. Some of these problems can be overcome
        by outputting the data as text files, at the expense of speed and file
        size.

        """
        raise NotImplementedError

    def tolist(self):
        """Return the array as a (possibly nested) list.

        Return a copy of the array data as a (nested) Python list.
        Data items are converted to the nearest compatible Python type.

        Parameters
        ----------
        none

        Returns
        -------
        y : list
            The possibly nested list of array elements.

        Notes
        -----
        The array may be recreated, ``a = np.array(a.tolist())``.

        Examples
        --------
        >>> a = np.array([1, 2])
        >>> a.tolist()
        [1, 2]
        >>> a = np.array([[1, 2], [3, 4]])
        >>> list(a)
        [array([1, 2]), array([3, 4])]
        >>> a.tolist()
        [[1, 2], [3, 4]]

        """
        raise NotImplementedError

    def tostring(self, order='C'):
        """Construct a Python string containing the raw data bytes in the array.

        Constructs a Python string showing a copy of the raw contents of
        data memory. The string can be produced in either 'C' or 'Fortran',
        or 'Any' order (the default is 'C'-order). 'Any' order means C-order
        unless the F_CONTIGUOUS flag in the array is set, in which case it
        means 'Fortran' order.

        Parameters
        ----------
        order : {'C', 'F', None}, optional
            Order of the data for multidimensional arrays:
            C, Fortran, or the same as for the original array.

        Returns
        -------
        s : str
            A Python string exhibiting a copy of `a`'s raw data.

        Examples
        --------
        >>> x = np.array([[0, 1], [2, 3]])
        >>> x.tostring()
        '\x00\x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x00\x03\x00\x00\x00'
        >>> x.tostring('C') == x.tostring()
        True
        >>> x.tostring('F')
        '\x00\x00\x00\x00\x02\x00\x00\x00\x01\x00\x00\x00\x03\x00\x00\x00'

        """
        raise NotImplementedError

    def trace(self, offset=0, axis1=0, axis2=1, dtype=None, out=None):
        """Return the sum along diagonals of the array.

        Refer to `numpy.trace` for full documentation.

        See Also
        --------
        numpy.trace : equivalent function

        """
        raise NotImplementedError

    def transpose(self, *axes):
        """Returns a view of the array with axes transposed.

        For a 1-D array, this has no effect. (To change between column and
        row vectors, first cast the 1-D array into a matrix object.)
        For a 2-D array, this is the usual matrix transpose.
        For an n-D array, if axes are given, their order indicates how the
        axes are permuted (see Examples). If axes are not provided and
        ``a.shape = (i[0], i[1], ... i[n-2], i[n-1])``, then
        ``a.transpose().shape = (i[n-1], i[n-2], ... i[1], i[0])``.

        Parameters
        ----------
        axes : None, tuple of ints, or `n` ints

         * None or no argument: reverses the order of the axes.

         * tuple of ints: `i` in the `j`-th place in the tuple means `a`'s
           `i`-th axis becomes `a.transpose()`'s `j`-th axis.

         * `n` ints: same as an n-tuple of the same ints (this form is
           intended simply as a "convenience" alternative to the tuple form)

        Returns
        -------
        out : ndarray
            View of `a`, with axes suitably permuted.

        See Also
        --------
        ndarray.T : Array property returning the array transposed.

        Examples
        --------
        >>> a = np.array([[1, 2], [3, 4]])
        >>> a
        array([[1, 2],
               [3, 4]])
        >>> a.transpose()
        array([[1, 3],
               [2, 4]])
        >>> a.transpose((1, 0))
        array([[1, 3],
               [2, 4]])
        >>> a.transpose(1, 0)
        array([[1, 3],
               [2, 4]])

        """
        ret = self[:]
        if self.ndim < 2:
            return ret
        if not axes:
            # empty axes tuple i.e. no axes were passed
            axes = np.arange(self.ndim)[::-1]
        elif len(axes) == 1:
            # we have either None or a tuple of ints
            if axes[0] is None:
                axes = np.arange(self.ndim)[::-1]
            elif isinstance(axes[0], tuple):
                axes = axes[0]
            else:
                raise ValueError, "invalid axis for this array"
        else:
            # assume axes is a tuple of ints
            axes = np.asarray(axes, dtype=np.int64)
        if len(axes) != self.ndim:
            raise ValueError, "axes don't match array"
        ret.global_slice = self.global_slice.transpose(axes)
        return ret

    def var(self, axis=None, dtype=None, out=None, ddof=0):
        """Returns the variance of the array elements, along given axis.

        Refer to `numpy.var` for full documentation.

        See Also
        --------
        numpy.var : equivalent function

        """
        raise NotImplementedError

    def view(self, dtype=None, type=None):
        """New view of array with the same data.

        Parameters
        ----------
        dtype : data-type, optional
            Data-type descriptor of the returned view, e.g., float32 or int16.
            The default, None, results in the view having the same data-type
            as `a`.
        type : Python type, optional
            Type of the returned view, e.g., ndarray or matrix.  Again, the
            default None results in type preservation.

        Notes
        -----
        ``a.view()`` is used two different ways:

        ``a.view(some_dtype)`` or ``a.view(dtype=some_dtype)`` constructs a view
        of the array's memory with a different data-type.  This can cause a
        reinterpretation of the bytes of memory.

        ``a.view(ndarray_subclass)`` or ``a.view(type=ndarray_subclass)`` just
        returns an instance of `ndarray_subclass` that looks at the same array
        (same shape, dtype, etc.)  This does not cause a reinterpretation of the
        memory.


        Examples
        --------
        >>> x = np.array([(1, 2)], dtype=[('a', np.int8), ('b', np.int8)])

        Viewing array data using a different type and dtype:

        >>> y = x.view(dtype=np.int16, type=np.matrix)
        >>> y
        matrix([[513]], dtype=int16)
        >>> print type(y)
        <class 'numpy.matrixlib.defmatrix.matrix'>

        Creating a view on a structured array so it can be used in calculations

        >>> x = np.array([(1, 2),(3,4)], dtype=[('a', np.int8), ('b', np.int8)])
        >>> xv = x.view(dtype=np.int8).reshape(-1,2)
        >>> xv
        array([[1, 2],
               [3, 4]], dtype=int8)
        >>> xv.mean(0)
        array([ 2.,  3.])

        Making changes to the view changes the underlying array

        >>> xv[0,1] = 20
        >>> print x
        [(1, 20) (3, 4)]

        Using a view to convert an array to a record array:

        >>> z = x.view(np.recarray)
        >>> z.a
        array([1], dtype=int8)

        Views share data:

        >>> x[0] = (9, 10)
        >>> z[0]
        (9, 10)

        """
        raise NotImplementedError

    ################################################################
    ### ndarray operator overloading
    ################################################################
    def __abs__(self):
        return abs(self)

    def __add__(self, y):
        return add(self,y)

    def __and__(self, y):
        return logical_and(self,y)

    def __array__(self, dtype):
        raise NotImplementedError

    def __array_finalize__(self, *args, **kwargs):
        raise NotImplementedError

    def __array_interface__(self, *args, **kwargs):
        raise NotImplementedError

    def __array_prepare__(self, *args, **kwargs):
        raise NotImplementedError

    def __array_priority__(self, *args, **kwargs):
        raise NotImplementedError

    def __array_struct__(self, *args, **kwargs):
        raise NotImplementedError

    def __array_wrap__(self, *args, **kwargs):
        raise NotImplementedError

    def __contains__(self, y):
        raise NotImplementedError

    def __copy__(self, order=None):
        """Return a copy of the array.

    Parameters
    ----------
    order : {'C', 'F', 'A'}, optional
        If order is 'C' (False) then the result is contiguous (default).
        If order is 'Fortran' (True) then the result has fortran order.
        If order is 'Any' (None) then the result has fortran order
        only if the array already is in fortran order.

        """
        raise NotImplementedError

    def __deepcopy__(self):
        raise NotImplementedError

    #def __delattr__

    def __delitem__(self, *args, **kwargs):
        raise ValueError, "cannot delete array elements"

    #def __delslice__

    def __div__(self, y):
        return divide(self,y)

    def __divmod__(self, y):
        t = mod(self,y)
        s = subtract(self,t)
        s = divide(s,y,s)
        return s,t

    def __eq__(self, y):
        return equal(self,y)

    def __float__(self):
        raise NotImplementedError

    def __floordiv__(self,y):
        return floor_divide(self,y)

    #def __format__

    def __ge__(self,y):
        return greater_equal(self,y)

    #def __getattribute__

    def __getitem__(self, key):
        # THIS IS A COLLECTIVE OPERATION
        if isinstance(key, (str,unicode)):
            raise NotImplementedError, "str or unicode key"
        if self.ndim == 0:
            raise IndexError, "0-d arrays can't be indexed"
        key = listify(key)
        fancy = False
        for arg in key:
            if isinstance(arg, (ndarray,np.ndarray,list,tuple)):
                fancy = True
                break
        if fancy:
            raise NotImplementedError, "TODO: fancy indexing"
        new_global_slice = self.global_slice[key]
        new_shape = new_global_slice.shape
        a = ndarray(new_shape, self.dtype, base=self)
        a.global_slice = new_global_slice
        if a.size == 1:
            a = a.allget()
            # convert single item to np.generic (scalar)
            a = a.dtype.type(a)
        return a

    #def __getslice__

    def __gt__(self,y):
        return greater(self,y)

    #def __hash__
    #def __hex__

    def __iadd__(self,y):
        return add(self,y,self)

    def __iand__(self,y):
        return logical_and(self,y,self)

    def __idiv__(self,y):
        return divide(self,y,self)

    def __ifloordiv__(self,y):
        return floor_divide(self,y,self)

    def __ilshift__(self,y):
        return left_shift(self,y,self)

    def __imod__(self,y):
        return mod(self,y,self)

    def __imul__(self,y):
        return multiply(self,y,self)

    #def __index__

    def __int__(self, *args, **kwargs):
        if self.size != 1:
            raise TypeError, "only length-1 arrays can be converted to Python scalars"
        return np.dtype(np.int32).type(self.allget())

    def __invert__(self):
        return invert(self)

    def __ior__(self,y):
        return logical_or(self,y,self)

    def __ipow__(self,y):
        return power(self,y,self)

    def __irshift__(self,y):
        return right_shift(self,y,self)

    def __isub__(self,y):
        return subtract(self,y,self)

    def __iter__(self, *args, **kwargs):
        raise NotImplementedError

    def __itruediv__(self,y):
        return true_divide(self,y,self)

    def __ixor__(self,y):
        return logical_xor(self,y,self)

    def __le__(self,y):
        return less_equal(self,y)

    def __len__(self):
        return self.shape[0]

    def __long__(self, *args, **kwargs):
        raise NotImplementedError

    def __lshift__(self,y):
        return left_shift(self,y)

    def __lt__(self,y):
        return less(self,y)

    def __mod__(self,y):
        return mod(self,y)

    def __mul__(self,y):
        return multiply(self,y)

    def __ne__(self,y):
        return not_equal(self,y)

    def __neg__(self):
        return negative(self)

    def __nonzero__(self):
        raise ValueError, ("The truth value of an array with more than one"
            " element is ambiguous. Use a.any() or a.all()")

    def __oct__(self, *args, **kwargs):
        raise NotImplementedError

    def __or__(self,y):
        return logical_or(self,y)

    def __pos__(self):
        return self.copy()

    def __pow__(self,y):
        return power(self,y)

    def __radd__(self,y):
        return add(y,self)

    def __rand__(self,y):
        return logical_and(y,self)

    def __rdiv__(self,y):
        return divide(y,self)

    def __rdivmod__(self,y):
        t = mod(y,self)
        s = subtract(y,t)
        s = divide(s,self,s)
        return s,t

    def __reduce__(self, *args, **kwargs):
        raise NotImplementedError

    def __reduce_ex__(self, *args, **kwargs):
        raise NotImplementedError

    def __repr__(self):
        result = ""
        if 0 == me():
            result = repr(self.get())
        return result

    def __rfloordiv__(self,y):
        return floor_divide(y,self)

    def __rlshift__(self,y):
        return left_shift(y,self)

    def __rmod__(self,y):
        return mod(y,self)

    def __rmul__(self,y):
        return multiply(y,self)

    def __ror__(self,y):
        return logical_or(y,self)

    def __rpow__(self,y):
        return power(y,self)

    def __rrshift__(self,y):
        return right_shift(y,self)

    def __rshift__(self,y):
        return right_shift(self,y)

    def __rsub__(self,y):
        return subtract(y,self)

    def __rtruediv__(self,y):
        return true_divide(y,self)

    def __rxor__(self,y):
        return logical_xor(y,self)

    #def __setattr__

    def __setitem__(self, key, value):
        # THIS IS A COLLECTIVE OPERATION
        sync()
        if isinstance(key, (str,unicode)):
            raise NotImplementedError, "str or unicode key"
        if self.ndim == 0:
            raise IndexError, "0-d arrays can't be indexed"
        key = listify(key)
        fancy = False
        for arg in key:
            if isinstance(arg, (ndarray,np.ndarray,list,tuple)):
                fancy = True
                break
        if fancy:
            raise NotImplementedError, "TODO: fancy indexing"
        global_slice = self.global_slice[key]
        new_shape = global_slice.shape
        value = asarray(value)
        npvalue = None
        release_value = False
        # access based on new global_slice as an ndarray first
        npself = self.access(global_slice)
        if npself is not None:
            if isinstance(value, ndarray):
                if (ga.compare_distr(value.handle, self.handle)
                        and value.global_slice == global_slice
                        and value.global_slice.T == self.global_slice.T):
                    # opt: same distributions and same slicing
                    # in practice this might not happen all that often
                    npvalue = value.access()
                    release_value = True
                else:
                    lo,hi = self.distribution()
                    result = global_slice.get_key(lo,hi)
                    result = util.broadcast_chomp(value.shape, result)
                    npvalue = value.get(result)
            elif isinstance(value, flatiter):
                raise NotImplementedError
            elif value.ndim > 0:
                    lo,hi = self.distribution()
                    result = global_slice.get_key(lo,hi)
                    result = util.broadcast_chomp(value.shape, result)
                    npvalue = value[result]
            else:
                npvalue = value
            npself[:] = npvalue
            if release_value:
                value.release()
            self.release_update()
        sync()

    #def __setslice__
    #def __setstate__
    #def __sizeof__
    
    def __str__(self):
        result = ""
        if 0 == me():
            result = str(self.get())
        return result

    def __sub__(self,y):
        return subtract(self,y)

    #def __subclasshook__

    def __truediv__(self,y):
        return true_divide(self,y)

    def __xor__(self,y):
        return logical_xor(self,y)

def _npin_piece_based_on_out(input, out, shape=None):
    # opt: same distributions and same slicing
    #   we can use local data exclusively
    #   in practice this might not happen all that often
    if (isinstance(input, ndarray)
            and ga.compare_distr(input.handle, out.handle)
            and input.global_slice == out.global_slice
            and input.global_slice.T == out.global_slice.T):
        return input.access(),True
    # no opt: requires copy of remote data
    elif shape is None or len(shape) > 0:
        lo,hi = out.distribution()
        result = out.global_slice.get_key(lo,hi)
        if shape is not None:
            result = util.broadcast_chomp(shape, result)
        if is_distributed(input):
            return input.get(result),False
        else:
            return input[result],False
    else:
        return input,False

class ufunc(object):
    """Functions that operate element by element on whole arrays.

    A detailed explanation of ufuncs can be found in the "ufuncs.rst"
    file in the NumPy reference guide.

    Unary ufuncs:
    =============

    op(X, out=None)
    Apply op to X elementwise

    Parameters
    ----------
    X : array_like
        Input array.
    out : array_like
        An array to store the output. Must be the same shape as `X`.

    Returns
    -------
    r : array_like
        `r` will have the same shape as `X`; if out is provided, `r`
        will be equal to out.

    Binary ufuncs:
    ==============

    op(X, Y, out=None)
    Apply `op` to `X` and `Y` elementwise. May "broadcast" to make
    the shapes of `X` and `Y` congruent.

    The broadcasting rules are:

    * Dimensions of length 1 may be prepended to either array.
    * Arrays may be repeated along dimensions of length 1.

    Parameters
    ----------
    X : array_like
        First input array.
    Y : array_like
        Second input array.
    out : array_like
        An array to store the output. Must be the same shape as the
        output would have.

    Returns
    -------
    r : array_like
        The return value; if out is provided, `r` will be equal to out.

    """
    def __init__(self, func):
        self.func = func
        self.__doc__ = func.__doc__

    def _get_identity(self):
        return self.func.identity
    identity = property(_get_identity)

    def _get_nargs(self):
        return self.func.nargs
    nargs = property(_get_nargs)

    def _get_nin(self):
        return self.func.nin
    nin = property(_get_nin)

    def _get_nout(self):
        return self.func.nout
    nout = property(_get_nout)

    def _get_ntypes(self):
        return self.func.ntypes
    ntypes = property(_get_ntypes)

    def _get_signature(self):
        return self.func.signature
    signature = property(_get_signature)

    def _get_types(self):
        return self.func.types
    types = property(_get_types)

    def __call__(self, *args, **kwargs):
        if self.func.nin == 1:
            return self._unary_call(*args, **kwargs)
        elif self.func.nin == 2:
            return self._binary_call(*args, **kwargs)
        else:
            raise ValueError, "only unary and binary ufuncs supported"

    def _unary_call(self, input, out=None, *args, **kwargs):
        input = asarray(input)
        input_shape = get_shape(input)
        input_dtype = get_dtype(input)
        if not (is_distributed(input) or is_distributed(out)):
            # no ndarray instances used, pass through immediately to numpy
            return self.func(input, out, *args, **kwargs)
        if out is None:
            # input must be an ndarray given previous conditionals
            # TODO okay, is there something better than this?
            ignore = np.ones(1, dtype=input_dtype)
            out_type = self.func(ignore).dtype
            out = ndarray(input_shape, out_type)
        # sanity checks
        if not is_array(out):
            raise TypeError, "return arrays must be of ArrayType"
        out_shape = get_shape(out)
        if input_shape != out_shape:
            # broadcasting doesn't apply to unary operations
            raise ValueError, 'invalid return array shape'
        # Now figure out what to do...
        if isinstance(out, ndarray):
            sync()
            # get out as an np.ndarray first
            npout = out.access()
            if npout is not None: # this proc owns data
                if input is out:
                    npin,release_in = npout,False
                else:
                    npin,release_in = _npin_piece_based_on_out(input,out)
                self.func(npin, npout, *args, **kwargs)
                if release_in:
                    input.release()
                out.release_update()
            #sync()
        elif isinstance(out, flatiter):
            sync()
            # first opt: input and out are same object
            #   we call _unary_call over again with the bases
            #   NOT SURE THAT THIS IS ACTUALLY OPTIMAL -- NEED TO TEST
            if input is out:
                self._unary_call(out.base, out.base, *args, **kwargs)
                return out.copy() # differs from NumPy (should be view)
            else:
                npout = out.access()
                if npout is not None: # this proc 'owns' data
                    lo,hi = out.distribution()
                    result = out.global_slice.get_key(lo,hi)
                    if is_distributed(input):
                        npin = input.get(result)
                    else:
                        npin = input[result]
                    self.func(npin, npout, *args, **kwargs)
                    out.release_update()
            #sync()
        else:
            sync()
            # out is not distributed
            npin = input
            if is_distributed(input):
                npin = input.allget()
            self.func(npin, out, *args, **kwargs)
            #sync() # I don't think we need this one
        return out

    def _binary_call(self, first, second, out=None, *args, **kwargs):
        first_isscalar = np.isscalar(first)
        second_isscalar = np.isscalar(second)
        # just in case
        first = asarray(first)
        second = asarray(second)
        if not (is_distributed(first)
                or is_distributed(second)
                or is_distributed(out)):
            # no ndarray instances used, pass through immediately to numpy
            return self.func(first, second, out, *args, **kwargs)
        first_dtype = get_dtype(first)
        second_dtype = get_dtype(second)
        first_shape = get_shape(first)
        second_shape = get_shape(second)
        if out is None:
            # first and/or second must be ndarrays given previous conditionals
            # TODO okay, is there something better than this?
            dtype = None
            if first_isscalar:
                if second_isscalar:
                    dtype = np.find_common_type([],[first_dtype,second_dtype])
                else:
                    dtype = np.find_common_type([second_dtype],[first_dtype])
            else:
                if second_isscalar:
                    dtype = np.find_common_type([first_dtype],[second_dtype])
                else:
                    dtype = np.find_common_type([first_dtype,second_dtype],[])
            shape = util.broadcast_shape(first_shape, second_shape)
            out = ndarray(shape, dtype)
        # sanity checks
        if not is_array(out):
            raise TypeError, "return arrays must be of ArrayType"
        # Now figure out what to do...
        if isinstance(out, ndarray):
            sync()
            # get out as an np.ndarray first
            npout = out.access()
            if npout is not None: # this proc owns data
                # get matching and compatible portions of input arrays
                # broadcasting rules (may) apply
                if first is out:
                    npfirst,release_first = npout,False
                else:
                    npfirst,release_first = _npin_piece_based_on_out(
                            first,out,first_shape)
                if second is first:
                    # zeroth opt: first and second are same object, so do the
                    # same thing for second that we did for first
                    npsecond,release_second = npfirst,False
                elif second is out:
                    npsecond,release_second = npout,False
                else:
                    npsecond,release_second = _npin_piece_based_on_out(
                            second,out,second_shape)
                self.func(npfirst, npsecond, npout, *args, **kwargs)
                if release_first:
                    first.release()
                if release_second:
                    second.release()
                out.release_update()
            #sync()
        elif isinstance(out, flatiter):
            sync()
            # first op: first and second and out are same object
            if first is second is out:
                self._binary_call(out.base,out.base,out.base,*args,**kwargs)
                return out.copy()
            else:
                npout = out.access()
                if npout is not None: # this proc 'owns' data
                    lo,hi = out.distribution()
                    result = out.global_slice.get_key(lo,hi)
                    if is_distributed(first):
                        npfirst = first.get(result)
                    else:
                        npfirst = first[result]
                    if second is first:
                        npsecond = npfirst
                    elif is_distributed(second):
                        npsecond = second.get(result)
                    else:
                        npsecond = second[result]
                    self.func(npfirst, npsecond, npout, *args, **kwargs)
                    out.release_update()
            #sync()
        else:
            sync()
            # out is not distributed
            ndfirst = first
            if is_distributed(first):
                ndfirst = first.allget()
            ndsecond = second
            if second is first:
                ndsecond = ndfirst
            elif is_distributed(second):
                ndsecond = second.allget()
            self.func(ndfirst, ndsecond, out, *args, **kwargs)
            #sync() # I don't think we need this one
        return out

    def reduce(self, a, axis=0, dtype=None, out=None, *args, **kwargs):
        """reduce(a, axis=0, dtype=None, out=None)

    Reduces `a`'s dimension by one, by applying ufunc along one axis.

    Let :math:`a.shape = (N_0, ..., N_i, ..., N_{M-1})`.  Then
    :math:`ufunc.reduce(a, axis=i)[k_0, ..,k_{i-1}, k_{i+1}, .., k_{M-1}]` =
    the result of iterating `j` over :math:`range(N_i)`, cumulatively applying
    ufunc to each :math:`a[k_0, ..,k_{i-1}, j, k_{i+1}, .., k_{M-1}]`.
    For a one-dimensional array, reduce produces results equivalent to:
    ::

     r = op.identity # op = ufunc
     for i in xrange(len(A)):
       r = op(r, A[i])
     return r

    For example, add.reduce() is equivalent to sum().

    Parameters
    ----------
    a : array_like
        The array to act on.
    axis : int, optional
        The axis along which to apply the reduction.
    dtype : data-type code, optional
        The type used to represent the intermediate results. Defaults
        to the data-type of the output array if this is provided, or
        the data-type of the input array if no output array is provided.
    out : ndarray, optional
        A location into which the result is stored. If not provided, a
        freshly-allocated array is returned.

    Returns
    -------
    r : ndarray
        The reduced array. If `out` was supplied, `r` is a reference to it.

    Examples
    --------
    >>> np.multiply.reduce([2,3,5])
    30

    A multi-dimensional array example:

    >>> X = np.arange(8).reshape((2,2,2))
    >>> X
    array([[[0, 1],
            [2, 3]],
           [[4, 5],
            [6, 7]]])
    >>> np.add.reduce(X, 0)
    array([[ 4,  6],
           [ 8, 10]])
    >>> np.add.reduce(X) # confirm: default axis value is 0
    array([[ 4,  6],
           [ 8, 10]])
    >>> np.add.reduce(X, 1)
    array([[ 2,  4],
           [10, 12]])
    >>> np.add.reduce(X, 2)
    array([[ 1,  5],
           [ 9, 13]])

        """
        if self.func.nin != 2:
            raise ValueError, "reduce only supported for binary functions"
        a = asarray(a)
        if not (isinstance(a, ndarray) or isinstance(out, ndarray)):
            # no ndarray instances used, pass through immediately to numpy
            return self.func.reduce(a, axis, dtype, out, *args, **kwargs)
        if axis < 0:
            axis += a.ndim
        if a.ndim < axis < 0:
            raise ValueError, "axis not in array"
        if out is None:
            shape = list(a.shape)
            del shape[axis]
            if dtype is None:
                dtype = a.dtype
            out = ndarray(shape, dtype=dtype)
        if out.ndim == 0:
            # optimize the 1d reduction
            nda = a.access()
            value = self.func.identity
            if nda is not None:
                value = self.func.reduce(nda)
            everything = comm().allgather(value)
            self.func.reduce(everything, out=out)
        elif out.ndim == 1 and 'OLD_REDUCE' in os.environ:
            # optimize the 2d reduction
            ndout = out.access()
            if ndout is not None:
                my_range = out.global_slice.bound_by_lohi(*out.distribution())
                if axis == 0: # rows
                    my_range = [slice(None,None,None)]+my_range
                elif axis == 1: # cols
                    my_range = my_range+[slice(None,None,None)]
                piece = a.get(my_range)
                self.func.reduce(piece, axis, dtype, ndout, *args, **kwargs)
        elif out.ndim == 1:
            # new algorithm using Isend/Irecv
            send_requests = []
            recv_requests = []
            nda = a.access()
            if nda is not None:
                TAG = 733823
                # we own part of the input
                # reduce it, then find where to send it
                local_out = self.func.reduce(nda, axis, dtype, *args, **kwargs)
                lo,hi = a.distribution()
                if axis == 0:
                    lo,hi = lo[1],hi[1]
                else:
                    lo,hi = lo[0],hi[0]
                m,p = ga.locate_region(out.handle, lo, hi)
                #print_sync("nda m=%s" % str(m))
                #print_sync("nda p=%s" % str(p))
                for i in range(len(p)):
                    rlo,rhi = m[i]
                    #print_sync("rlo,rhi=%s,%s" % (rlo,rhi))
                    rlo = rlo - lo
                    rhi = rhi - lo
                    #print_sync("again rlo,rhi=%s,%s" % (rlo,rhi))
                    send_requests.append(
                            comm().Isend(local_out[rlo:rhi], p[i], TAG))
            ndout = out.access()
            buf = None
            if ndout is not None:
                lo,hi = out.distribution()
                lo,hi = lo[0],hi[0]
                #print_sync("lo,hi=%s,%s" % (lo,hi))
                shape = ga.inquire_dims(a.handle)
                #print_sync("shape=%s" % shape)
                olo = [lo,lo]
                ohi = [hi,hi]
                #print_sync("olo,ohi=%s,%s" % (olo,ohi))
                olo[axis],ohi[axis] = 0,shape[axis]
                #print_sync("again olo,ohi=%s,%s" % (olo,ohi))
                m,p = ga.locate_region(a.handle, olo, ohi)
                #print_sync("nda m=%s" % str(m))
                #print_sync("nda p=%s" % str(p))
                buf = np.ndarray((len(p),hi-lo), dtype=out.dtype)
                buf[:] = self.func.identity
                for i in range(len(p)):
                    rlo,rhi = m[i]
                    #print_sync("rlo,rhi=%s,%s" % (rlo,rhi))
                    if axis == 0:
                        rlo,rhi = rlo[1],rhi[1]
                    else:
                        rlo,rhi = rlo[0],rhi[0]
                    #print_sync("again rlo,rhi=%s,%s" % (rlo,rhi))
                    rlo = rlo - lo
                    rhi = rhi - lo
                    #print_sync("yet again rlo,rhi=%s,%s" % (rlo,rhi))
                    recv_requests.append(
                            comm().Irecv(buf[i][rlo:rhi], p[i], TAG))
            while send_requests or recv_requests:
                if MPI.Request.Testall(send_requests):
                    send_requests = []
                if MPI.Request.Testall(recv_requests):
                    recv_requests = []
            if ndout is not None:
                self.func.reduce(buf, axis=0, out=ndout)
                out.release_update()
            if nda is not None:
                a.release()
        else:
            slicer = [slice(0,None,None)]*a.ndim
            axis_iterator = iter(xrange(a.shape[axis]))
            # copy first loop iteration to 'out'
            slicer[axis] = axis_iterator.next()
            out[:] = a[slicer]
            # remaining loop iterations are appropriately reduced
            for i in axis_iterator:
                slicer[axis] = i
                ai = a[slicer]
                self.__call__(out,ai,out)
        return out

    def accumulate(self, a, axis=0, dtype=None, out=None, *args, **kwargs):
        if self.func.nin != 2:
            raise ValueError, "accumulate only supported for binary functions"
        a = asarray(a)
        if not (isinstance(a, ndarray) or isinstance(out, ndarray)):
            # no ndarray instances used, pass through immediately to numpy
            return self.func.accumulate(a, axis, dtype, out, *args, **kwargs)
        if axis < 0:
            axis += a.ndim
        if a.ndim < axis < 0:
            raise ValueError, "axis not in array"
        if out is None:
            if dtype is None:
                dtype = a.dtype
            out = ndarray(a.shape, dtype=dtype)
        if out.ndim == 1:
            # optimize the 1d accumulate
            if isinstance(out, ndarray):
                ndout = out.access()
                if ndout is not None:
                    result = out.global_slice.get_key(*out.distribution())
                    lo,hi = result[0].start,result[0].stop
                    if is_distributed(a):
                        piece = a.get(result)
                    else:
                        piece = a[lo:hi]
                    self.func.accumulate(piece, out=ndout)
                    # probably more efficient to use allgather and exchange last
                    # values among all procs. We also need ordering information,
                    # so we exchange the 'lo' value.
                    everything = comm().allgather((ndout[-1],lo))
                    reduction = self.func.identity
                    for lvalue,llo in everything:
                        if lvalue is not None and llo < lo:
                            reduction = self.func(reduction,lvalue)
                    self.func(ndout,reduction,ndout)
                else:
                    everything = comm().allgather((None,None))
            else:
                raise NotImplementedError
        else:
            slicer_i = [slice(0,None,None)]*a.ndim
            slicer_i_1 = [slice(0,None,None)]*a.ndim
            axis_iterator = iter(xrange(a.shape[axis]))
            # copy first loop iteration to 'out'
            slicer_i[axis] = axis_iterator.next()
            out[slicer_i] = a[slicer_i]
            # remaining loop iterations are appropriately accumulated
            for i in axis_iterator:
                slicer_i[axis] = i
                slicer_i_1[axis] = i-1
                x = out[slicer_i_1]
                y = a[slicer_i]
                z = out[slicer_i]
                self.__call__(x,y,z)
                #self.__call__(out[slicer_i_1],a[slicer_i],out[slicer_i])
        return out

    def outer(self, *args, **kwargs):
        if self.func.nin != 2:
            raise ValueError, "outer product only supported for binary functions"
        raise NotImplementedError

    def reduceat(self, *args, **kwargs):
        if self.func.nin != 2:
            raise ValueError, "reduceat only supported for binary functions"
        raise NotImplementedError

# unary ufuncs
abs = ufunc(np.abs)
absolute = abs
arccos = ufunc(np.arccos)
arccosh = ufunc(np.arccosh)
arcsin = ufunc(np.arcsin)
arcsinh = ufunc(np.arcsinh)
arctan = ufunc(np.arctan)
arctanh = ufunc(np.arctanh)
bitwise_not = ufunc(np.bitwise_not)
ceil = ufunc(np.ceil)
conj = ufunc(np.conj)
conjugate = ufunc(np.conjugate)
cos = ufunc(np.cos)
cosh = ufunc(np.cosh)
deg2rad = ufunc(np.deg2rad)
degrees = ufunc(np.degrees)
exp = ufunc(np.exp)
exp2 = ufunc(np.exp2)
expm1 = ufunc(np.expm1)
fabs = ufunc(np.fabs)
floor = ufunc(np.floor)
frexp = ufunc(np.frexp)
invert = ufunc(np.invert)
isfinite = ufunc(np.isfinite)
isinf = ufunc(np.isinf)
isnan = ufunc(np.isnan)
log = ufunc(np.log)
log2 = ufunc(np.log2)
log10 = ufunc(np.log10)
log1p = ufunc(np.log1p)
logical_not = ufunc(np.logical_not)
modf = ufunc(np.modf)
negative = ufunc(np.negative)
rad2deg = ufunc(np.rad2deg)
radians = ufunc(np.radians)
reciprocal = ufunc(np.reciprocal)
rint = ufunc(np.rint)
sign = ufunc(np.sign)
signbit = ufunc(np.signbit)
sin = ufunc(np.sin)
sinh = ufunc(np.sinh)
spacing = ufunc(np.spacing)
sqrt = ufunc(np.sqrt)
square = ufunc(np.square)
tan = ufunc(np.tan)
tanh = ufunc(np.tanh)
trunc = ufunc(np.trunc)
# binary ufuncs
add = ufunc(np.add)
arctan2 = ufunc(np.arctan2)
bitwise_and = ufunc(np.bitwise_and)
bitwise_or = ufunc(np.bitwise_or)
bitwise_xor = ufunc(np.bitwise_xor)
copysign = ufunc(np.copysign)
divide = ufunc(np.divide)
equal = ufunc(np.equal)
floor_divide = ufunc(np.floor_divide)
fmax = ufunc(np.fmax)
fmin = ufunc(np.fmin)
fmod = ufunc(np.fmod)
greater = ufunc(np.greater)
greater_equal = ufunc(np.greater_equal)
hypot = ufunc(np.hypot)
ldexp = ufunc(np.ldexp)
left_shift = ufunc(np.left_shift)
less = ufunc(np.less)
less_equal = ufunc(np.less_equal)
logaddexp = ufunc(np.logaddexp)
logaddexp2 = ufunc(np.logaddexp2)
logical_and = ufunc(np.logical_and)
logical_or = ufunc(np.logical_or)
logical_xor = ufunc(np.logical_xor)
maximum = ufunc(np.maximum)
minimum = ufunc(np.minimum)
mod = ufunc(np.mod)
multiply = ufunc(np.multiply)
nextafter = ufunc(np.nextafter)
not_equal = ufunc(np.not_equal)
power = ufunc(np.power)
remainder = ufunc(np.remainder)
right_shift = ufunc(np.right_shift)
subtract = ufunc(np.subtract)
true_divide = ufunc(np.true_divide)

class flatiter(object):
    """Flat iterator object to iterate over arrays.

    A `flatiter` iterator is returned by ``x.flat`` for any array `x`.
    It allows iterating over the array as if it were a 1-D array,
    either in a for-loop or by calling its `next` method.

    Iteration is done in C-contiguous style, with the last index varying the
    fastest. The iterator can also be indexed using basic slicing or
    advanced indexing.

    See Also
    --------
    ndarray.flat : Return a flat iterator over an array.
    ndarray.flatten : Returns a flattened copy of an array.

    Notes
    -----
    A `flatiter` iterator can not be constructed directly from Python code
    by calling the `flatiter` constructor.

    Examples
    --------
    >>> x = np.arange(6).reshape(2, 3)
    >>> fl = x.flat
    >>> type(fl)
    <type 'numpy.flatiter'>
    >>> for item in fl:
    ...     print item
    ...
    0
    1
    2
    3
    4
    5

    >>> fl[2:4]
    array([2, 3])

    """
    def __init__(self, base):
        self._base = base
        self._index = 0
        self._values = None
        size = base.global_slice.size
        self.global_slice = util.MasterKey([size])
        size_per_proc = size//nproc()
        remainder = size%nproc() # remainder gets distributed
        if me() < remainder:
            self._lo = [size_per_proc*me() + me()]
            self._hi = [size_per_proc*me() + me() + size_per_proc + 1]
        else:
            self._lo = [size_per_proc*me() + remainder]
            self._hi = [size_per_proc*me() + remainder + size_per_proc]
        # if lo,hi are the same, this indicates no data on this proc
        if self._lo == self._hi:
            self._lo = [-1]
            self._hi = [-1]
    
    def _get_base(self):
        return self._base
    base = property(_get_base)

    def _get_coords(self):
        return np.unravel_index(self._index, self._base.shape)
    coords = property(_get_coords)

    def copy(self):
        """Get a copy of the iterator as a 1-D array.

        Examples
        --------
        >>> x = np.arange(6).reshape(2, 3)
        >>> x
        array([[0, 1, 2],
               [3, 4, 5]])
        >>> fl = x.flat
        >>> fl.copy()
        array([0, 1, 2, 3, 4, 5])

        """
        return self._base.flatten()
    
    def _get_index(self):
        return self._index
    index = property(_get_index)

    def distribution(self):
        return self._lo,self._hi

    def owns(self):
        return self._hi[0]>=0

    def get(self, key=None):
        # THIS OPERATION IS ONE-SIDED
        # we expect key to be a single value or a slice, but might also be an
        # iterable of length 1
        if key is None:
            # TODO this doesn't feel right -- communicate and then a local copy
            # operation? The answer is correct, but seems like too much extra
            # work.
            return self._base.get().flatten()
        if isinstance(key, (list,tuple)):
            key = key[0]
        if isinstance(key, (slice,util.RangeKey)):
            # get shape of base's global_slice
            shape = self._base.global_slice.shape
            # create index coordinates
            if isinstance(key, util.RangeKey):
                key = key.pyobj()
            i = (np.indices(shape).reshape(len(shape),-1).T)[key]
            return self._base.gather(i)
        else:
            key = long(key)
            index = np.unravel_index(key, self._base.shape)
            return self._base.gather(index)

    def allget(self, key=None):
        """Like get(), but when all processes need the same piece.
        
        This operation is collective.
        
        """
        # TODO it's not clear whether this approach is better than having all
        # P processors ga.get() the same piece.
        if not me():
            result = self.get(key)
            return comm().bcast(result)
        else:
            return comm().bcast()

    def __getitem__(self, key):
        # THIS OPERATION IS COLLECTIVE
        sync()
        if isinstance(key, (list,tuple)):
            if len(key) > 1:
                raise IndexError, "unsupported iterator index"
        # TODO optimize the gather since this is a collective
        return self.allget(key)

    def __len__(self):
        return self.global_slice.size

    def put(self, key, value):
        # THIS OPERATION IS ONE-SIDED
        # we expect key to be a single value or a slice, but might also be an
        # iterable of length 1
        if isinstance(key, (list,tuple)):
            key = key[0]
        if isinstance(key, (slice,util.RangeKey)):
            if isinstance(key, util.RangeKey):
                key = key.pyobj()
            # get shape of global_slice
            shape = self._base.global_slice.shape
            # create index coordinates
            i = (np.indices(shape).reshape(len(shape),-1).T)[key]
            value = asarray(value)
            values = None
            if value.size == 1:
                values = np.zeros(len(i), dtype=self._base.dtype)
                values[:] = value
            else:
                assert value.ndim == 1 and len(value) == len(i)
                if is_distributed(value):
                    values = value.get()
                else:
                    values = value
            ga.scatter(self._base.handle, values, i)
        else:
            key = long(key)
            i = np.unravel_index(key, self._base.shape)
            ga.scatter(self._base.handle, value, i)

    def __setitem__(self, key, value):
        # THIS OPERATION IS COLLECTIVE
        sync()
        if self.owns():
            global_slice = self.global_slice[key]
            if isinstance(key, (slice,util.RangeKey)):
                try:
                    my_range = global_slice.bound_by_lohi(self._lo,self._hi)
                except IndexError:
                    sync()
                    return
                value = asarray(value)
                if value.ndim == 1:
                    self.put(my_range, value[my_range])
                else:
                    self.put(my_range, value)
            else:
                self.put(key, value)
        sync()

    def next(self):
        if self._index < self._len:
            tmp = self.base[self.coords]
            self._index += 1
            return tmp
        else:
            raise StopIteration

    def access(self):
        """Return a copy of a 'local' portion."""
        if self._values is not None:
            raise ValueError, "call release or release_update before access"
        if not self.owns():
            self._values = None
        else:
            lo,hi = self.distribution()
            result = self.global_slice.get_key(lo,hi)
            self._values = self.get(result)
        return self._values

    def release(self):
        self._values = None

    def release_update(self):
        if self._values is not None:
            lo,hi = self.distribution()
            result = self.global_slice.get_key(lo,hi)
            self.put(result, self._values)
        self._values = None

    def __repr__(self):
        result = ""
        if 0 == me():
            result = repr(self.get())
        return result

    def __str__(self):
        result = ""
        if 0 == me():
            result = str(self.get())
        return result

def zeros(shape, dtype=np.float, order='C'):
    """zeros(shape, dtype=float, order='C')

    Return a new array of given shape and type, filled with zeros.

    Parameters
    ----------
    shape : int or sequence of ints
        Shape of the new array, e.g., ``(2, 3)`` or ``2``.
    dtype : data-type, optional
        The desired data-type for the array, e.g., `numpy.int8`.  Default is
        `numpy.float64`.
    order : {'C', 'F'}, optional
        Whether to store multidimensional data in C- or Fortran-contiguous
        (row- or column-wise) order in memory.

    Returns
    -------
    out : ndarray
        Array of zeros with the given shape, dtype, and order.

    See Also
    --------
    zeros_like : Return an array of zeros with shape and type of input.
    ones_like : Return an array of ones with shape and type of input.
    empty_like : Return an empty array with shape and type of input.
    ones : Return a new array setting values to one.
    empty : Return a new uninitialized array.

    Examples
    --------
    >>> np.zeros(5)
    array([ 0.,  0.,  0.,  0.,  0.])

    >>> np.zeros((5,), dtype=numpy.int)
    array([0, 0, 0, 0, 0])

    >>> np.zeros((2, 1))
    array([[ 0.],
           [ 0.]])

    >>> s = (2,2)
    >>> np.zeros(s)
    array([[ 0.,  0.],
           [ 0.,  0.]])

    >>> np.zeros((2,), dtype=[('x', 'i4'), ('y', 'i4')]) # custom dtype
    array([(0, 0), (0, 0)],
          dtype=[('x', '<i4'), ('y', '<i4')])

    """
    if not should_distribute(shape):
        return np.zeros(shape, dtype, order)
    a = ndarray(shape, dtype)
    buf = a.access()
    if buf is not None:
        buf[:] = 0
        a.release_update()
    return a

def asarray(a, dtype=None, order=None):
    if isinstance(a, (ndarray,flatiter,np.ndarray,np.generic)):
        # we return ga.gain.ndarray instances for obvious reasons, but we
        # also return numpy.ndarray instances because they already exist in
        # whole on all procs -- no need to distribute pieces
        return a
    else:
        npa = np.asarray(a, dtype=dtype)
        if should_distribute(npa.size):
            g_a = ndarray(npa.shape, npa.dtype, npa)
            return g_a # distributed using Global Arrays ndarray
        else:
            return npa # possibly a scalar or zero rank array

def from_ga(g_a):
    """Create ndarray from a global array."""
    dtype = ga.inquire_dtype(g_a)
    shape = ga.inquire_dims(g_a)
    return ndarray(shape, dtype, base=g_a)

cdef extern MPI_Comm GA_MPI_Comm_pgroup_default()
def comm():
    """Returns the MPI_Comm instance associated with the process group."""
    cdef MPI.Comm communicator = MPI.Comm()
    communicator.ob_mpi = GA_MPI_Comm_pgroup_default()
    return communicator

#cdef bint DEBUG = False
#cdef bint DEBUG_SYNC = False
#
#def set_debug(val):
#    global DEBUG
#    DEBUG = val
#
#def set_debug_sync(val):
#    global DEBUG_SYNC
#    DEBUG_SYNC = val
#
#def print_debug(s):
#    if DEBUG:
#        print s
#
def print_sync(what):
    #if DEBUG:
    if True:
        sync()
        if 0 == me():
            print "[0] %s" % str(what)
            for proc in xrange(1,nproc()):
                data = comm().recv(source=proc, tag=11)
                print "[%d] %s" % (proc, str(data))
        else:
            comm().send(what, dest=0, tag=11)
        sync()
