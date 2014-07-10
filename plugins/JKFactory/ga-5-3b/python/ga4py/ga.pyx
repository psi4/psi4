#cython: embedsignature=True
"""The Global Arrays (GA) Python interface.

This module exports the GA C API, with a few enhancements.  The notable
exceptions include supporting Pythonic ranges.  The ranges here are half-open
e.g. [0,n) instead of in the C API where they are closed e.g. [0,n].  This
follows the Python convention.

"""
# keep the ga functions alphabetical since this is going to be a huge file!

from libc.stdlib cimport malloc,free
from gah cimport *
import numpy as np
cimport numpy as np
from cpython.ref cimport Py_INCREF
from cpython.ref cimport PyTypeObject
import __builtin__

DEF EXCLUSIVE = 0

np.import_array()

cdef extern from "numpy/arrayobject.h":
    object PyArray_NewFromDescr(PyTypeObject* subtype,
                                np.dtype descr,
                                int nd,
                                np.npy_intp * dims,
                                np.npy_intp * strides,
                                void * data,
                                int flags,
                                object obj)

cdef bint _initialized = False

TYPE_BASE  = 1000
C_CHAR     = (TYPE_BASE + 0)
C_INT      = (TYPE_BASE + 1)
C_LONG     = (TYPE_BASE + 2)
C_FLOAT    = (TYPE_BASE + 3)
C_DBL      = (TYPE_BASE + 4)
C_LDBL     = (TYPE_BASE + 5)
C_SCPL     = (TYPE_BASE + 6)
C_DCPL     = (TYPE_BASE + 7)
C_LDCPL    = (TYPE_BASE + 8)
F_BYTE     = (TYPE_BASE + 9)
F_INT      = (TYPE_BASE + 10)
F_LOG      = (TYPE_BASE + 11)
F_REAL     = (TYPE_BASE + 12)
F_DBL      = (TYPE_BASE + 13)
F_SCPL     = (TYPE_BASE + 14)
F_DCPL     = (TYPE_BASE + 15)
C_LONGLONG = (TYPE_BASE + 16)

WORLD_PROC_GROUP = -1

_to_dtype = {
        C_CHAR:     np.dtype(np.int8),
        C_INT:      np.dtype(np.int32),
        C_LONG:     np.dtype(np.int64),
        C_LONGLONG: np.dtype(np.int64),
        C_FLOAT:    np.dtype(np.float32),
        C_DBL:      np.dtype(np.float64),
        C_SCPL:     np.dtype(np.complex64),
        C_DCPL:     np.dtype(np.complex128),
        }
# numpy doesn't always have these types depending on the system
cdef bint float128_in_np = ('float128' in dir(np))
cdef bint complex256_in_np = ('complex256' in dir(np))
if float128_in_np:
    _to_dtype[C_LDBL] = np.dtype(np.float128)
if complex256_in_np:
    _to_dtype[C_LDCPL] = np.dtype(np.complex256)

#############################################################################
# utility functions
#############################################################################

def dtype(int gatype):
    """Converts the given GA type to a numpy dtype."""
    if gatype in _to_dtype:
        return _to_dtype[gatype]
    raise ValueError, "%d was not a recognized GA type" % gatype

def inquire_dtype(int g_a):
    """Returns the numpy dtype of the given GA."""
    gatype = inquire_type(g_a)
    return dtype(gatype)

cdef inline void* _gapy_malloc(size_t bytes, int align, char *name):
    """Wrapper around C stdlib malloc()."""
    return malloc(bytes)

cdef inline void _gapy_free(void *ptr):
    """Wrapper around C stdlib free()."""
    free(ptr)

cdef inline np.ndarray[np.int32_t, ndim=1] _inta32(array_like):
    """Converts an integer array-like to an ndarray of 32bit integers.

    Functions which take a dimension shape or subscript can use this to
    convert what the user passes to a numpy.ndarray using numpy.asarray.

    As a convenience, single values can be passed as well.

    :Parameters:
        array_like : integer array-like

    :returns: The converted array_like to an ndarray.

    """
    cdef np.ndarray[np.int32_t, ndim=1] array_like_nd
    try:
        array_like_nd = np.asarray(array_like, dtype=np.int32)
    except ValueError: # try again in case array_like is a single value
        array_like_nd = np.asarray([array_like], dtype=np.int32)
    return array_like_nd

cdef inline np.ndarray[np.int64_t, ndim=1] _inta64(array_like):
    """Converts an integer array-like to an ndarray of 64bit integers.

    Functions which take a dimension shape or subscript can use this to
    convert what the user passes to a numpy.ndarray using numpy.asarray.

    As a convenience, single values can be passed as well.

    :Parameters:
        array_like : integer array-like

    :returns: The converted array_like to an ndarray.

    """
    cdef np.ndarray[np.int64_t, ndim=1] array_like_nd
    try:
        array_like_nd = np.asarray(array_like, dtype=np.int64)
    except ValueError: # try again in case array_like is a single value
        array_like_nd = np.asarray([array_like], dtype=np.int64)
    return array_like_nd

cdef inline _lohi(int g_a, lo, hi):
    """Converts and/or prepares a lo/hi combination.

    Functions which take a patch specification can use this to convert the
    given lo and/or hi into ndarrays using numpy.asarray.

    * If neither lo nor hi is given, lo is replaced with an array of zeros and
      hi is replaced with the last index in each dimension (i.e. the shape).

    * If only lo is given, hi is replaced with lo. In other words, this is
      a single value.

    * It is an error to specify hi without lo.

    :Parameters:
        g_a : int
            the array handle
        lo : 1D array-like
            lower bounds of a slice
        hi : 1D array-like
            upper bounds of a slice

    :returns: The converted lo and hi ndarrays.

    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd
    cdef int ndim = GA_Ndim(g_a)
    if lo is None and hi is None:
        lo_nd = np.zeros((ndim), dtype=np.int64)
        hi_nd = inquire_dims(g_a)
    elif lo is not None and hi is None:
        lo_nd = _inta64(lo)
        hi_nd = lo_nd+1
    elif lo is None and hi is not None:
        raise ValueError, 'lo cannot be None if hi is None'
    else: # lo and hi are not None
        lo_nd = _inta64(lo)
        hi_nd = _inta64(hi)
    if len(lo_nd) != ndim:
        raise ValueError, 'len(lo_nd) != ndim; len(%s) != %s' % (lo_nd,ndim)
    if len(hi_nd) != ndim:
        raise ValueError, 'len(hi_nd) != ndim; len(%s) != %s' % (hi_nd,ndim)
    # We must make a copy of hi_nd. If user passes in an ndarray, the
    # following "prep" operation will change the user's 'hi'.
    #hi_nd -= 1 # <----- don't do that!
    hi_nd = hi_nd-1 # prep hi for GA's inclusive indexing
    return lo_nd,hi_nd

cdef void* _convert_multiplier(int gtype, value,
        int *iv, long *lv, long long *llv,
        float *fv, double *dv, long double *ldv,
        SingleComplex *fcv, DoubleComplex *dcv):
    """Returns the address of an appropriately converted value.
    
    Functions which take an alpha/beta/value need to have the value
    appropriately converted from the (possible) Python type to a C type. Often
    the GA function takes a void* in this case, so the address of the
    converted value is returned.

    """
    cdef float complex pfcv=1.0
    cdef double complex pdcv=1.0
    if value is None:
        raise ValueError, "cannot convert None"
    if gtype == C_INT:
        iv[0] = value
        return iv
    elif gtype == C_LONG:
        lv[0] = value
        return lv
    elif gtype == C_LONGLONG:
        llv[0] = value
        return llv
    elif gtype == C_FLOAT:
        fv[0] = value
        return fv
    elif gtype == C_DBL:
        dv[0] = value
        return dv
    elif gtype == C_LDBL and float128_in_np:
        ldv[0] = value
        return ldv
    elif gtype == C_SCPL:
        pfcv = value
        fcv[0].real = pfcv.real
        fcv[0].imag = pfcv.imag
        return fcv
    elif gtype == C_DCPL:
        pdcv = value
        dcv[0].real = pdcv.real
        dcv[0].imag = pdcv.imag
        return dcv
    else:
        raise TypeError, "type of g_a not recognized"

def zip(lo, hi):
    """Transforms a GA lo,hi combination into a slice list."""
    return [slice(l,h) for l,h in __builtin__.zip(lo,hi)]

#############################################################################
# GA API
#############################################################################

def abs_value(int g_a, lo=None, hi=None):
    """Take element-wise absolute value of the array or patch.
    
    This is a collective operation.

    :Parameters:
        g_a : int
            the array handle
        lo : 1D array-like
            lower bound patch coordinates, inclusive
        hi : 1D array-like
            higher bound patch coordinates, exclusive
    
    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd
    if lo is None and hi is None:
        GA_Abs_value(g_a)
    else:
        lo_nd,hi_nd = _lohi(g_a,lo,hi)
        GA_Abs_value_patch64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data)

def acc(int g_a, buffer, lo=None, hi=None, alpha=None):
    """Combines data from buffer with data in the global array patch.
    
    The buffer array is assumed to be have the same number of
    dimensions as the global array.  If the buffer is not contiguous, a
    contiguous copy will be made.
    
        global array section (lo[],hi[]) += alpha * buffer

    This is a one-sided and atomic operation.

    :Parameters:
        g_a : int
            the array handle
        buffer : array-like
            must be contiguous and have same number of elements as patch
        lo : 1D array-like
            lower bound patch coordinates, inclusive
        hi : 1D array-like
            higher bound patch coordinates, exclusive
        alpha : object
            multiplier (converted to appropriate type)

    """
    _acc_common(g_a, buffer, lo, hi, alpha)

cdef _acc_common(int g_a, buffer, lo=None, hi=None, alpha=None,
        bint nb=False, bint periodic=False, skip=None):
    """Combines data from buffer with data in the global array patch.
    
    The local array is assumed to have the same shape as the requested region,
    or the local array can be 1-dimensional so long as it has the same number
    of elements as the requested region. Any detected inconsitencies raise a
    ValueError.
    
        global array section (lo[],hi[]) += alpha * buffer

    This is a one-sided and atomic operation.

    :Parameters:
        g_a : int
            the array handle
        buffer : array-like
            the data to put;
            should either be 1D and len(buffer)==np.prod(hi-lo), or
            np.all(buffer.shape == hi-lo) i.e. buffer is 1D and same size as
            requested region or buffer is the same shape as requested region
        lo : 1D array-like
            lower bound patch coordinates, inclusive
        hi : 1D array-like
            higher bound patch coordinates, exclusive
        alpha : object
            multiplier (converted to appropriate type)
        nb : bool
            whether the call is non-blocking
        periodic : bool
            whether the call is periodic
        skip : 1D array-like
            strides for each dimension

    :see: nbacc
    :see: periodic_acc
    :see: strided_acc

    :returns: None, however if nb=True, the nonblocking handle is returned.

    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd, ld_nd, shape, skip_nd
    cdef np.ndarray buffer_nd
    cdef int gtype=inquire_type(g_a)
    cdef int            ialpha
    cdef long           lalpha
    cdef long long      llalpha
    cdef float          falpha
    cdef double         dalpha
    cdef long double    ldalpha
    cdef SingleComplex  fcalpha
    cdef DoubleComplex  dcalpha
    cdef void          *valpha=NULL
    cdef ga_nbhdl_t     nbhandle
    dtype = _to_dtype[gtype]
    lo_nd,hi_nd = _lohi(g_a,lo,hi)
    shape = hi_nd-lo_nd+1
    if skip is None:
        skip_nd = None
    else:
        skip_nd = _inta64(skip)
        shape = (hi_nd-lo_nd)/skip_nd+1
    buffer_nd = np.asarray(buffer, dtype=dtype)
    if buffer_nd.dtype != dtype:
        raise ValueError, "buffer is wrong type :: buffer=%s != %s" % (
                buffer.dtype, dtype)
    # Due to GA restrictions, buffer must not have negative strides
    # and buffer's last stride must be same as itemsize
    strides = [buffer_nd.strides[i]/buffer_nd.itemsize
            for i in range(buffer_nd.ndim)]
    if (strides and (strides[-1] != 1 or np.any(np.asarray(strides) < 0))):
        buffer_nd = np.ascontiguousarray(buffer_nd)
    # we allow 1-d "flat" buffers in addition to buffers matching the shape of
    # the requested region
    if buffer_nd.ndim == 1:
        if buffer_nd.size != np.prod(shape):
            raise ValueError, ('buffer size does not match shape :: '
                    'buffer.size=%s != np.prod(shape)=%s' % (
                    buffer_nd.size, np.prod(shape)))
        ld_nd = shape[1:]
    else:
        buffer_shape = [buffer_nd.shape[i] for i in range(buffer_nd.ndim)]
        if not np.all(buffer_shape == shape):
            raise ValueError, ('buffer shape does not match request shape :: '
                    'buffer_shape=%s != shape=%s' % (
                    buffer_shape, shape))
        ld_nd = np.asarray([strides[i]/strides[i+1]
                for i in range(buffer_nd.ndim-1)], dtype=np.int64)
    if alpha is None:
        alpha = 1
    valpha = _convert_multiplier(gtype, alpha,
            &ialpha,  &lalpha,  &llalpha,
            &falpha,  &dalpha,  &ldalpha,
            &fcalpha, &dcalpha)
    if nb:
        NGA_NbAcc64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <void*>buffer_nd.data, <int64_t*>ld_nd.data, valpha, &nbhandle)
        return nbhandle
    elif periodic:
        NGA_Periodic_acc64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <void*>buffer_nd.data, <int64_t*>ld_nd.data, valpha)
    elif skip is not None:
        NGA_Strided_acc64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <int64_t*>skip_nd.data,
                <void*>buffer_nd.data, <int64_t*>ld_nd.data, valpha)
    else:
        NGA_Acc64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <void*>buffer_nd.data, <int64_t*>ld_nd.data, valpha)

def access(int g_a, lo=None, hi=None, int proc=-1):
    """Returns local array patch.
    
    This routine allows to access directly, in place elements in the local
    section of a global array. It useful for writing new GA operations.
    If no patch is specified, the entire local patch is returned.  If this
    process does not own any data, None is returned.
    
    Note: The entire local data is always accessed, but if a smaller patch is
    requested, an appropriately sliced ndarray is returned.

    If proc is not specified, then ga.nodeid() is used.

    Each call to ga.access has to be followed by a call to either ga.release
    or ga.release_update. You can access in this fashion only local data.
    Since the data is shared with other processes, you need to consider issues
    of mutual exclusion.

    This operation is local. 

    :Parameters:
        g_a : int
            the array handle
        lo : 1D array-like
            lower bound patch coordinates, inclusive
        hi : 1D array-like
            higher bound patch coordinates, exclusive
        proc : int
            defaults to ga.nodeid(), but can specify a proc within the same
            SMP node to access its data instead
    
    :returns: ndarray representing local patch

    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd
    cdef np.ndarray[np.int64_t, ndim=1] ld_nd, lo_dst, hi_dst, dims_nd
    cdef int i, gtype=inquire_type(g_a)
    cdef int dimlen=GA_Ndim(g_a)
    cdef np.dtype dtype = _to_dtype[gtype]
    cdef void *ptr
    cdef np.npy_intp *dims = NULL
    cdef np.npy_intp *strides = NULL
    if proc < 0:
        proc = GA_Nodeid()
    # first things first, if no data is owned, return None
    lo_dst,hi_dst = distribution(g_a, proc)
    if lo_dst[0] < 0 or hi_dst[0] < 0:
        return None
    # always access the entire local data
    ld_nd = np.ones(dimlen, dtype=np.int64)
    # put hi_dst back to GA inclusive indexing convention
    hi_dst -= 1
    NGA_Access64(g_a, <int64_t*>lo_dst.data, <int64_t*>hi_dst.data, &ptr,
            <int64_t*>ld_nd.data)
    # put hi_dst back to Python exclusive indexing convention
    hi_dst += 1
    dims_nd = hi_dst-lo_dst
    # fix the strides, for example if ghost cells are in use
    # or if memory is allocated in a non-contiguous way
    if dimlen >= 3:
        for i in range(-3, -dimlen-1, -1):
            ld_nd[i] *= ld_nd[i+1]
    ld_nd *= _to_dtype[gtype].itemsize
    # must convert int64_t ndarray shape to npy_intp array
    dims = <np.npy_intp*>malloc(dimlen*sizeof(np.npy_intp))
    strides = <np.npy_intp*>malloc(dimlen*sizeof(np.npy_intp))
    for i in range(dimlen):
        dims[i] = dims_nd[i]
        strides[i] = ld_nd[i]
    Py_INCREF(dtype)
    array = PyArray_NewFromDescr(<PyTypeObject*>np.ndarray,
            dtype, dimlen, dims, strides, ptr, np.NPY_WRITEABLE, None)
    free(dims)
    free(strides)
    if lo is not None or hi is not None:
        if lo is not None:
            lo_nd = _inta64(lo)
        else:
            lo_nd = lo_dst
        if hi is not None:
            hi_nd = _inta64(hi)
        else:
            hi_nd = hi_dst
        # sanity checks
        if np.sometrue(lo_nd>hi_nd):
            raise ValueError,"lo>hi lo=%s hi=%s"%(lo_nd,hi_nd)
        if np.sometrue(lo_nd<lo_dst):
            raise ValueError,"lo out of bounds lo_dst=%s lo=%s"%(lo_dst,lo_nd)
        if np.sometrue(hi_nd>hi_dst):
            raise ValueError,"hi out of bounds hi_dst=%s hi=%s"%(hi_dst,hi_nd)
        slices = []
        for i in range(dimlen):
            slices.append(slice(lo_nd[i]-lo_dst[i],hi_nd[i]-lo_dst[i]))
        return array[slices]
    return array

def access_block(int g_a, int idx):
    """Returns local array patch for a block-cyclic distribution.
    
    This routine allows to access directly, in place elements in the local
    section of a global array. It useful for writing new GA operations.
    
    Each call to ga.access_block has to be followed by a call to either
    ga.release_block or ga.release_update_block. You can access in this
    fashion only local data.  Since the data is shared with other processes,
    you need to consider issues of mutual exclusion.

    This operation is local. 

    :Parameters:
        g_a : int
            the array handle
        idx : int
            the block index

    :returns: ndarray representing local block

    """
    cdef np.ndarray[np.int64_t, ndim=1] ld_nd, lo_dst, hi_dst, dims_nd
    cdef int i, gtype=inquire_type(g_a)
    cdef int dimlen=GA_Ndim(g_a), typenum=_to_dtype[gtype].num
    cdef void *ptr
    cdef np.npy_intp *dims = NULL
    # first things first, if no data is owned, return None
    lo_dst,hi_dst = distribution(g_a, idx)
    if lo_dst[0] < 0 or hi_dst[0] < 0:
        return None
    # put hi_dst back to GA inclusive indexing convention
    hi_dst -= 1
    # always access the entire local data
    ld_nd = np.zeros(dimlen-1, dtype=np.int64)
    NGA_Access_block64(g_a, idx, &ptr, <int64_t*>ld_nd.data)
    dims_nd = hi_dst-lo_dst+1
    # must convert int64_t ndarray shape to npy_intp array
    dims = <np.npy_intp*>malloc(dimlen*sizeof(np.npy_intp))
    for i in range(dimlen):
        dims[i] = dims_nd[i]
    array = np.PyArray_SimpleNewFromData(dimlen, dims, typenum, ptr)
    free(dims)
    return array

def access_block_grid(int g_a, subscript):
    """Returns local array patch for a SCALAPACK block-cyclic distribution.

    The subscript array contains the subscript of the block in the array of
    blocks. This subscript is based on the location of the block in a grid,
    each of whose dimensions is equal to the number of blocks that fit along
    that dimension.

    Each call to ga.access_block_grid has to be followed by a call to either
    ga.release_block_grid or ga.release_update_block_grid. You can access in
    this fashion only local data.  Since the data is shared with other
    processes, you need to consider issues of mutual exclusion.

    This operation is local. 

    :Parameters:
        g_a : int
            the array handle
        subscript : 1D array-like
            subscript of the block in the array

    :returns: ndarray representing local block

    """
    raise NotImplementedError

def access_block_segment(int g_a, int proc):
    """This function can be used to gain access to the all the locally held
    data on a particular processor that is associated with a block-cyclic
    distributed array.

    The data  inside this segment has a lot of additional structure so this
    function is not generally useful to developers. It is primarily used
    inside the GA library to implement other GA routines. Each call to
    ga.access_block_segment should be followed by a call to either
    ga.release_block_segment or ga.release_update_block_segment.

    This is a local operation.

    :Parameters:
        g_a : int
            the array handle
        proc : int
            processor ID

    :returns: ndarray representing local block

    """
    cdef int64_t elems
    cdef int gtype=inquire_type(g_a)
    cdef int typenum=_to_dtype[gtype].num
    cdef void *ptr
    cdef np.npy_intp *dims = NULL
    # always access the entire local data
    NGA_Access_block_segment64(g_a, proc, &ptr, &elems)
    # must convert int64_t ndarray shape to npy_intp array
    dims = <np.npy_intp*>malloc(sizeof(np.npy_intp))
    dims[0] = elems
    array = np.PyArray_SimpleNewFromData(1, dims, typenum, ptr)
    free(dims)
    return array

def access_ghost_element(int g_a, subscript):
    """Returns a scalar ndarray representing the requested ghost element.

    This function can be used to return a pointer to any data element in the
    locally held portion of the global array and can be used to directly
    access ghost cell data. The array subscript refers to the local index of
    the  element relative to the origin of the local patch (which is assumed
    to be indexed by (0,0,...)).

    This is a  local operation. 

    :Parameters:
        g_a : int
            the array handle
        subscript : 1D array-like of integers
            index of the desired element

    :returns: ndarray scalar representing local block

    """
    raise NotImplementedError, "use access_ghosts(g_a) instead"

def access_ghosts(int g_a):
    """Returns ndarray representing local patch with ghost cells.

    This routine will provide access to the ghost cell data residing on each
    processor. Calls to NGA_Access_ghosts should normally follow a call to
    NGA_Distribution  that returns coordinates of the visible data patch
    associated with a processor. You need to make sure that the coordinates of
    the patch are valid (test values returned from NGA_Distribution).

    You can only access local data.

    This operation is local.

    :Parameters:
        g_a : int
            the array handle

    :returns: ndarray scalar representing local block with ghost cells

    """
    cdef np.ndarray[np.int64_t, ndim=1] ld_nd, lo_dst, hi_dst, dims_nd
    cdef int i, gtype=inquire_type(g_a)
    cdef int dimlen=GA_Ndim(g_a), typenum=_to_dtype[gtype].num
    cdef void *ptr
    cdef np.npy_intp *dims = NULL
    # first things first, if no data is owned, return None
    lo_dst,hi_dst = distribution(g_a)
    if lo_dst[0] < 0 or hi_dst[0] < 0:
        return None
    # always access the entire local data
    dims_nd = np.zeros(dimlen, dtype=np.int64)
    ld_nd = np.zeros(dimlen-1, dtype=np.int64)
    NGA_Access_ghosts64(g_a, <int64_t*>dims_nd.data, &ptr, <int64_t*>ld_nd.data)
    # must convert int64_t ndarray shape to npy_intp array
    dims = <np.npy_intp*>malloc(dimlen*sizeof(np.npy_intp))
    for i in range(dimlen):
        dims[i] = dims_nd[i]
    array = np.PyArray_SimpleNewFromData(dimlen, dims, typenum, ptr)
    free(dims)
    return array

def add(int g_a, int g_b, int g_c, alpha=None, beta=None, alo=None, ahi=None,
        blo=None, bhi=None, clo=None, chi=None):
    """Element-wise addition of two arrays.

    The arrays must be the same shape and identically aligned.
    The result (c) may replace one of the input arrays (a/b).
    Patches of arrays (which must have the same number of elements) may also
    be added together elementw=-wise, if patch coordinates are specified.
    c = alpha*a + beta*b

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle
        g_b : int
            the array handle
        g_c : int
            the array handle
        alpha : object
            multiplier (converted to appropriate type)
        beta : object
            multiplier (converted to appropriate type)
        alo : 1D array-like of integers
            lower bound patch coordinates of g_a, inclusive
        ahi : 1D array-like of integers
            higher bound patch coordinates of g_a, exclusive
        blo : 1D array-like of integers
            lower bound patch coordinates of g_b, inclusive
        bhi : 1D array-like of integers
            higher bound patch coordinates of g_b, exclusive
        clo : 1D array-like of integers
            lower bound patch coordinates of g_c, inclusive
        chi : 1D array-like of integers
            higher bound patch coordinates of g_c, exclusive


    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    cdef np.ndarray[np.int64_t, ndim=1] clo_nd, chi_nd
    cdef int gtype=inquire_type(g_a)
    cdef int            ialpha,  ibeta
    cdef long           lalpha,  lbeta
    cdef long long      llalpha, llbeta
    cdef float          falpha,  fbeta
    cdef double         dalpha,  dbeta
    cdef long double    ldalpha, ldbeta
    cdef SingleComplex  fcalpha, fcbeta
    cdef DoubleComplex  dcalpha, dcbeta
    cdef void          *valpha, *vbeta
    if alpha is None:
        alpha = 1
    valpha = _convert_multiplier(gtype, alpha,
            &ialpha,  &lalpha,  &llalpha,
            &falpha,  &dalpha,  &ldalpha,
            &fcalpha, &dcalpha)
    if beta is None:
        beta = 1
    vbeta = _convert_multiplier(gtype, beta,
            &ibeta,  &lbeta,  &llbeta,
            &fbeta,  &dbeta,  &ldbeta,
            &fcbeta, &dcbeta)
    if (alo is None and ahi is None
            and blo is None and bhi is None
            and clo is None and chi is None):
        GA_Add(valpha, g_a, vbeta, g_b, g_c)
    else:
        alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
        blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
        clo_nd,chi_nd = _lohi(g_c,clo,chi)
        NGA_Add_patch64(
                valpha, g_a, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                 vbeta, g_b, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data,
                        g_c, <int64_t*>clo_nd.data, <int64_t*>chi_nd.data)

def add_constant(int g_a, alpha, lo=None, hi=None):
    """Adds the constant alpha to each element of the array. 

    This operation is collective.

    :Parameters:
        g_a : int
            the array handle
        alpha : object
            the constant to add (converted to appropriate type)
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : 1D array-like of integers
            higher bound patch coordinates, exclusive

    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd
    cdef int gtype=inquire_type(g_a)
    cdef int            ialpha
    cdef long           lalpha
    cdef long long      llalpha
    cdef float          falpha
    cdef double         dalpha
    cdef long double    ldalpha
    cdef SingleComplex  fcalpha
    cdef DoubleComplex  dcalpha
    cdef void          *valpha
    valpha = _convert_multiplier(gtype, alpha,
            &ialpha,  &lalpha,  &llalpha,
            &falpha,  &dalpha,  &ldalpha,
            &fcalpha, &dcalpha)
    if lo is None and hi is None:
        GA_Add_constant(g_a, valpha)
    else:
        lo_nd,hi_nd = _lohi(g_a,lo,hi)
        GA_Add_constant_patch64(
                g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data, valpha)

def add_diagonal(int g_a, int g_v):
    """Adds the elements of the vector g_v to the diagonal of matrix g_a.

    This operation is collective.

    :Parameters:
        g_a : int
            the array handle
        g_v : int
            the vector handle

    """
    GA_Add_diagonal(g_a, g_v)

def allocate(int g_a):
    """Allocates memory for the handle obtained using ga.create_handle.

    At a minimum, the ga.set_data function must be called before the memory is
    allocated. Other ga.set_xxx functions can also be called before invoking
    this function.

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle

    :returns: True if allocation of g_a was successful.

    """
    if GA_Allocate(g_a) == 1:
        return True
    return False

def brdcst(buffer, int root=0):
    """Broadcast from process root to all other processes.

    If the buffer is not contiguous, an error is raised.  This operation is
    provided only for convenience purposes: it is available regardless of the
    message-passing library that GA is running with.

    This is a collective operation. 

    :Parameters:
        buffer : 1D array-like of objects
            the ndarray message (converted to the appropriate type)
        root : int
            the process which is sending

    :returns: The buffer in case a temporary was passed in.

    """
    cdef np.ndarray buffer_nd
    buffer_nd = np.asarray(buffer)
    if not buffer_nd.flags['C_CONTIGUOUS']:
        raise ValueError, "the buffer must be contiguous"
    #if buffer_nd.ndim != 1:
    #    raise ValueError, "the buffer must be one-dimensional"
    GA_Brdcst(buffer_nd.data, buffer_nd.size*buffer_nd.itemsize, root)
    return buffer_nd

def check_handle(int g_a, char *message):
    """Checks that the array handle g_a is valid.
    
    If not, calls ga.error withe the provided string.

    This operation is local.

    """
    GA_Check_handle(g_a, message)

def cluster_nnodes():
    """Returns the total number of nodes that the program is running on.

    On SMP architectures, this will be less than or equal to the total number
    of processors.

    This is a  local operation.

    """
    return GA_Cluster_nnodes()

def cluster_nodeid(int proc=-1):
    """Returns the node ID of this process or the given process.

    On SMP architectures with more than one processor per node, several
    processes may return the same node id.

    This is a local operation.

    :Parameters:
        proc : int
            process ID to lookup

    """
    if proc >= 0:
        return GA_Cluster_proc_nodeid(proc)
    return GA_Cluster_nodeid()

def cluster_proc_nodeid(int proc):
    """Returns the node ID of the specified process.

    On SMP architectures with more than one processor per node, several
    processors may return the same node id.

    This is a local operation.

    """
    return GA_Cluster_proc_nodeid(proc)

def cluster_nprocs(int node):
    """Returns the number of processors available on the given node.

    This is a local operation.

    """
    return GA_Cluster_nprocs(node)

def cluster_procid(int node, int proc):
    """Returns the proc ID associated with node and local proc ID.

    If node has N processors, then the value of proc lies between 0 and
    N-1.

    This is a  local operation. 

    """
    return GA_Cluster_procid(node, proc)

def compare_distr(int g_a, int g_b):
    """Compares the distributions of two global arrays.

    This is a collective operation.

    :returns: True if distributions are identical and False when they are not

    """
    if GA_Compare_distr(g_a, g_b) == 0:
        return True
    return False

def copy(int g_a, int g_b, alo=None, ahi=None, blo=None, bhi=None,
        bint trans=False):
    """Copies elements from array g_a into array g_b.

    For the operation over the entire arrays, the arrays must be the same
    type, shape, and identically aligned.  No transpose is allowed in this
    case.

    For patch operations, the patches of arrays may be of different shapes but
    must have the same number of elements. Patches must be nonoverlapping (if
    g_a=g_b).  Transposes are allowed for patch operations.

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle copying from
        g_b : int
            the array handle copying to
        alo : 1D array-like of integers
            lower bound patch coordinates of g_a, inclusive
        ahi : 1D array-like of integers
            higher bound patch coordinates of g_a, exclusive
        blo : 1D array-like of integers
            lower bound patch coordinates of g_b, inclusive
        bhi : 1D array-like of integers
            higher bound patch coordinates of g_b, exclusive
        trans : bool
            whether the transpose operator should be applied True=applied
             
    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    cdef char trans_c
    if alo is None and ahi is None and blo is None and bhi is None:
        GA_Copy(g_a, g_b)
    else:
        alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
        blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
        if trans:
            trans_c = "T"
        else:
            trans_c = "N"
        NGA_Copy_patch64(trans_c,
                g_a, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                g_b, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data)

def create(int gtype, dims, char *name="", chunk=None, int pgroup=-1):
    """Creates an n-dimensional array using the regular distribution model.

    The array can be distributed evenly or not. The control over the
    distribution is accomplished by specifying chunk (block) size for all or
    some of array dimensions. For example, for a 2-dimensional array, setting
    chunk[0]=dim[0] gives distribution by vertical strips (chunk[0]*dims[0]);
    setting chunk[1]=dim[1] gives distribution by horizontal strips
    (chunk[1]*dims[1]). Actual chunks will be modified so that they are at
    least the size of the minimum and each process has either zero or one
    chunk. Specifying chunk[i] as <1 will cause that dimension to be
    distributed evenly.

    As a convenience, when chunk is omitted or None, the entire array is
    distributed evenly.

    :Parameters:
        gtype : int
            the type of the array
        dims : 1D array-like of integers
            shape of the array
        name : string
            the name of the array
        chunk : 1D array-like of integers
            see above
        pgroup : int
            create array only as part of this processor group

    :returns: a non-zero array handle means the call was succesful.

    This is a collective operation. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] dims_nd, chunk_nd=None
    dims_nd = _inta64(dims)
    if pgroup < 0:
        pgroup = pgroup_get_default()
    if chunk:
        chunk_nd = _inta64(chunk)
        return NGA_Create_config64(gtype, len(dims_nd), <int64_t*>dims_nd.data,
                name, <int64_t*>chunk_nd.data, pgroup)
    else:
        return NGA_Create_config64(gtype, len(dims_nd), <int64_t*>dims_nd.data,
                name, NULL, pgroup)

def create_ghosts(int gtype, dims, width, char *name="", chunk=None,
        int pgroup=-1):
    """Creates an array with a layer of ghost cells around the visible data.

    The array can be distributed evenly or not evenly. The control over the
    distribution is accomplished by specifying chunk (block) size for all or
    some of the array dimensions. For example, for a 2-dimensional array,
    setting chunk(1)=dim(1) gives distribution by vertical strips
    (chunk(1)*dims(1)); setting chunk(2)=dim(2) gives distribution by
    horizontal strips (chunk(2)*dims(2)). Actual chunks will be modified so
    that they are at least the size of the minimum and each process has either
    zero or one chunk. Specifying chunk(i) as <1 will cause that dimension
    (i-th) to be distributed evenly. The  width of the ghost cell layer in
    each dimension is specified using the array width().  The local data of
    the global array residing on each processor will have a layer width[n]
    ghosts cells wide on either side of the visible data along the dimension
    n.

    :Parameters:
        gtype : int
            the type of the array
        dims : 1D array-like of integers
            shape of the array
        width : 1D array-like of integers
            ghost cell widths
        name : string
            the name of the array
        chunk : 1D array-like of integers
            see above
        pgroup : int
            create array only as part of this processor group

    :returns: a non-zero array handle means the call was successful.

    This is a collective operation. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] dims_nd, chunk_nd, width_nd
    dims_nd = _inta64(dims)
    width_nd = _inta64(width)
    if pgroup < 0:
        pgroup = pgroup_get_default()
    if chunk:
        chunk_nd = _inta64(chunk)
        return NGA_Create_ghosts_config64(gtype, len(dims_nd),
                <int64_t*>dims_nd.data, <int64_t*>width_nd.data, name,
                <int64_t*>chunk_nd.data, pgroup)
    else:
        return NGA_Create_ghosts_config64(gtype, len(dims_nd),
                <int64_t*>dims_nd.data, <int64_t*>width_nd.data, name,
                NULL, pgroup)

def create_handle():
    """Returns a global array handle that can be used to create a new array.
    
    The sequence of operations is to begin with a call to ga.create_handle to
    get a new array handle. The attributes of the array, such as dimension,
    size, type, etc. can then be set using successive calls to the ga.set_xxx
    subroutines. When all array attributes have been set, the ga.allocate
    subroutine is called and the global array is actually created and memory
    for it is allocated.

    This is a collective operation.

    """
    return GA_Create_handle()

def create_irreg(int gtype, dims, block, map, char *name="", int pgroup=-1):
    """Creates an array by following the user-specified distribution.

    The distribution is specified as a Cartesian product of distributions for
    each dimension. The array indices start at 0. For example, the following
    figure demonstrates distribution of a 2-dimensional array 8x10 on 6 (or
    more) processors. nblock[2]=[3,2], the size of map array is s=5 and array
    map contains the following elements map=[0,2,6, 0, 5]. The distribution is
    nonuniform because, P1 and P4 get 20 elements each and processors
    P0,P2,P3, and P5 only 10 elements each.

    This is a collective operation.

    :Parameters:
        gtype : int
            the type of the array
        dims : 1D array-like of integers
            shape of the array
        block : 1D array-like of integers
            the number of blocks each dimension is divided into
        map : 1D array-like of integers
            starting index for each block
            len(map) == sum of all elements of nblock array
        name : string
            the name of the array
        pgroup : int
            create array only as part of this processor group
    
    :returns: integer handle representing the array; a non-zero value indicates success

    """
    cdef np.ndarray[np.int64_t, ndim=1] dims_nd, block_nd, map_nd
    dims_nd = _inta64(dims)
    block_nd = _inta64(block)
    map_nd = _inta64(map)
    if pgroup < 0:
        pgroup = pgroup_get_default()
    return NGA_Create_irreg_config64(gtype, len(dims_nd),
            <int64_t*>dims_nd.data, name,
            <int64_t*>block_nd.data, <int64_t*>map_nd.data, pgroup)

def create_ghosts_irreg(int gtype, dims, width, block, map, char *name="",
        int pgroup=-1):
    """Creates an array with a layer of ghost cells around the visible data.

    The distribution is specified as a Cartesian product of distributions for
    each dimension. For example, the following figure demonstrates
    distribution of a 2-dimensional array 8x10 on 6 (or more) processors.
    nblock(2)=[3,2], the size of map array is s=5 and array map contains the
    following elements map=[1,3,7, 1, 6]. The distribution is nonuniform
    because, P1 and P4 get 20 elements each and processors P0,P2,P3, and P5
    only 10 elements each. 

    The array width[] is used to control the width of the ghost cell boundary
    around the visible data on each processor. The local data of the global
    array residing on each processor will have a layer width[n] ghosts cells
    wide on either side of the visible data along the dimension n.

    This is a collective operation. 

    :Parameters:
        gtype : int
            the type of the array
        dims : 1D array-like of integers
            shape of the array
        width : 1D array-like of integers
            ghost cell widths
        block : 1D array-like of integers
            number of blocks each dimension is divided into
        map : 1D array-like of integers
            starting index for each block
            len(map) == sum of all elements of nblock array
        name : string
            the name of the array
        pgroup : int
            create array only as part of this processor group
    
    :returns: a non-zero array handle means the call was succesful

    """
    cdef np.ndarray[np.int64_t, ndim=1] dims_nd, width_nd, block_nd, map_nd
    dims_nd = _inta64(dims)
    width_nd = _inta64(width)
    block_nd = _inta64(block)
    map_nd = _inta64(map)
    if pgroup < 0:
        pgroup = pgroup_get_default()
    return NGA_Create_ghosts_irreg_config64(gtype, len(dims_nd),
            <int64_t*>dims_nd.data, <int64_t*>width_nd.data, name,
            <int64_t*>block_nd.data, <int64_t*>map_nd.data, pgroup)

def create_mutexes(int number):
    """Creates a set containing the number of mutexes.

    Mutex is a simple synchronization object used to protect Critical
    Sections. Only one set of mutexes can exist at a time. Array of mutexes
    can be created and destroyed as many times as needed.

    Mutexes are numbered: 0, ..., number -1.

    This is a collective operation. 

    :Parameters:
        number : int
            the number of mutexes to create

    :returns: True on success, False on failure

    """
    if GA_Create_mutexes(number) == 1:
        return True
    return False

def deregister_type(int type):
    """Removes the data type previously registered using register_type.

    :Parameters:
        type : int
            the data type handle

    """
    return NGA_Deregister_type(type)
    
def destroy(int g_a):
    """Deallocates the array and frees any associated resources.

    This is a collective operation.

    """
    GA_Destroy(g_a)

def destroy_mutexes():
    """Destroys the set of mutexes created with ga_create_mutexes.
    
    :returns: True if the operation succeeded; False when failed

    This is a collective operation. 

    """
    if GA_Destroy_mutexes() == 1:
        return True
    return False

def diag(int g_a, int g_s, int g_v, evalues=None):
    """Solve the generalized eigen-value problem.

    The input matrices are not overwritten or destroyed.
    
    :Parameters:
        g_a : int
            the array handle of the matrix to diagonalize
        g_s : int
            the array handle of the metric
        g_v : int
            the array handle to return evecs

    :returns: All eigen-values as an ndarray in ascending order.

    This is a collective operation. 

    """
    if evalues is None:
        gtype,dims = inquire(g_a)
        evalues = np.ndarray((dims[0]), dtype=_to_dtype(gtype))
    else:
        evalues = np.asarray(evalues)
    GA_Diag(g_a, g_s, g_v, <void*>evalues.data)
    return evalues

def diag_reuse(int control, int g_a, int g_s, int g_v, evalues=None):
    """Solve the generalized eigen-value problem.

    Recommended for REPEATED calls if g_s is unchanged.
    The input matrices are not overwritten or destroyed.
    
    :Parameters:
        control : int
            0 indicates first call to the eigensolver;
            >0 consecutive calls (reuses factored g_s);
            <0 only erases factorized g_s; g_v and eval unchanged
            (should be called after previous use if another
            eigenproblem, i.e., different g_a and g_s, is to
            be solved) 
        g_a : int
            the array handle of the matrix to diagonalize
        g_s : int
            the array handle of the metric
        g_v : int
            the array handle to return evecs

    :returns: All eigen-values as an ndarray in ascending order.

    This is a collective operation. 

    """
    if evalues is None:
        gtype,dims = inquire(g_a)
        evalues = np.ndarray((dims[0]), dtype=_to_dtype(gtype))
    else:
        evalues = np.asarray(evalues)
    GA_Diag_reuse(control, g_a, g_s, g_v, <void*>evalues.data)
    return evalues

def diag_std(int g_a, int g_v, evalues=None):
    """Solve the standard (non-generalized) eigenvalue problem.

    The input matrix is neither overwritten nor destroyed.
    
    :Parameters:
        g_a : int
            the array handle of the matrix to diagonalize
        g_v : int
            the array handle to return evecs

    :returns: all eigenvectors via the g_v global array, and eigenvalues as an ndarray in ascending order

    This is a collective operation. 

    """
    if evalues is None:
        gtype,dims = inquire(g_a)
        evalues = np.ndarray((dims[0]), dtype=_to_dtype(gtype))
    else:
        evalues = np.asarray(evalues)
    GA_Diag_std(g_a, g_v, <void*>evalues.data)
    return evalues

cpdef distribution(int g_a, int proc=-1):
    """Return the distribution given to proc.

    If proc is not specified, then ga.nodeid() is used.  The range is
    returned as -1 for lo and -2 for hi if no elements are owned by
    proc.
    
    """
    cdef int ndim = GA_Ndim(g_a)
    cdef np.ndarray[np.int64_t, ndim=1] lo = np.zeros((ndim), dtype=np.int64)
    cdef np.ndarray[np.int64_t, ndim=1] hi = np.zeros((ndim), dtype=np.int64)
    if proc < 0:
        proc = GA_Nodeid()
    NGA_Distribution64(g_a, proc, <int64_t*>lo.data, <int64_t*>hi.data)
    # convert hi to python exclusive indexing convetion
    hi += 1
    return lo,hi

def dot(int g_a, int g_b, alo=None, ahi=None, blo=None, bhi=None,
        bint ta=False, bint tb=False):
    """Computes the element-wise dot product of two arrays.

    Arrays must be of the same type and same number of elements.
    Patch operation allows for possibly transposed patches.

    This is a collective operation.

    :Parameters:
        g_a : int
            the array handle
        g_b : int
            the array handle
        alo : 1D array-like of integers
            lower bound patch coordinates of g_a, inclusive
        ahi : 1D array-like of integers
            higher bound patch coordinates of g_a, exclusive
        blo : 1D array-like of integers
            lower bound patch coordinates of g_b, inclusive
        bhi : 1D array-like of integers
            higher bound patch coordinates of g_b, exclusive
        ta : bool
            whether the transpose operator should be applied to g_a True=applied
        tb : bool
            whether the transpose operator should be applied to g_b True=applied

    :returns: SUM_ij a(i,j)*b(i,j)

    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    cdef char ta_c, tb_c
    cdef int gtype=inquire_type(g_a)
    cdef float complex pfcv
    cdef double complex pdcv
    cdef SingleComplex gfcv
    cdef DoubleComplex gdcv
    if alo is None and ahi is None and blo is None and bhi is None:
        if gtype == C_INT:
            return GA_Idot(g_a, g_b)
        elif gtype == C_LONG:
            return GA_Ldot(g_a, g_b)
        elif gtype == C_LONGLONG:
            return GA_Lldot(g_a, g_b)
        elif gtype == C_FLOAT:
            return GA_Fdot(g_a, g_b)
        elif gtype == C_DBL:
            return GA_Ddot(g_a, g_b)
        elif gtype == C_SCPL:
            gfcv = GA_Cdot(g_a, g_b)
            pfcv.real = gfcv.real
            pfcv.imag = gfcv.imag
            return pfcv
        elif gtype == C_DCPL:
            gdcv = GA_Zdot(g_a, g_b)
            pdcv.real = gdcv.real
            pdcv.imag = gdcv.imag
            return pdcv
        else:
            raise TypeError
    else:
        alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
        blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
        if ta:
            ta_c = "T"
        else:
            ta_c = "N"
        if tb:
            tb_c = "T"
        else:
            tb_c = "N"
        if gtype == C_INT:
            return NGA_Idot_patch64(
                    g_a, ta_c, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                    g_b, tb_c, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data)
        elif gtype == C_LONG:
            return NGA_Ldot_patch64(
                    g_a, ta_c, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                    g_b, tb_c, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data)
        elif gtype == C_LONGLONG:
            return NGA_Lldot_patch64(
                    g_a, ta_c, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                    g_b, tb_c, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data)
        elif gtype == C_FLOAT:
            return NGA_Fdot_patch64(
                    g_a, ta_c, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                    g_b, tb_c, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data)
        elif gtype == C_DBL:
            return NGA_Ddot_patch64(
                    g_a, ta_c, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                    g_b, tb_c, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data)
        elif gtype == C_SCPL:
            gfcv = NGA_Cdot_patch64(
                    g_a, ta_c, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                    g_b, tb_c, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data)
            pfcv.real = gfcv.real
            pfcv.imag = gfcv.imag
            return pfcv
        elif gtype == C_DCPL:
            gdcv = NGA_Zdot_patch64(
                    g_a, ta_c, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                    g_b, tb_c, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data)
            pdcv.real = gdcv.real
            pdcv.imag = gdcv.imag
            return pdcv
        else:
            raise TypeError
    
def duplicate(int g_a, char *name=""):
    """Creates a new array by applying all the properties of another existing
    array.
    
    :Parameters:
        g_a : int
            the array handle
        name : string
            the new name of the created array

    :returns: a non-zero array handle means the call was succesful.

    This is a collective operation. 

    """
    return GA_Duplicate(g_a, name)

def elem_divide(int g_a, int g_b, int g_c, alo=None, ahi=None, blo=None,
        bhi=None, clo=None, chi=None):
    """Computes the element-wise quotient of the two arrays.

    Arrays or array patches must be of the same types and same number of
    elements. For two-dimensional arrays:

                c(i, j)  = a(i,j)/b(i,j)

    The result (c) may replace one of the input arrays (a/b).
    If one of the elements of array g_b is zero, the quotient for the element
    of g_c will be set to GA_NEGATIVE_INFINITY. 

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle
        g_b : int
            the array handle
        g_c : int
            the array handle
        alo : 1D array-like of integers
            lower bound patch coordinates of g_a, inclusive
        ahi : 1D array-like of integers
            higher bound patch coordinates of g_a, exclusive
        blo : 1D array-like of integers
            lower bound patch coordinates of g_b, inclusive
        bhi : 1D array-like of integers
            higher bound patch coordinates of g_b, exclusive
        clo : 1D array-like of integers
            lower bound patch coordinates of g_c, inclusive
        chi : 1D array-like of integers
            higher bound patch coordinates of g_c, exclusive

    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    cdef np.ndarray[np.int64_t, ndim=1] clo_nd, chi_nd
    if (alo is None and ahi is None
            and blo is None and bhi is None
            and clo is None and chi is None):
        GA_Elem_divide(g_a, g_b, g_c)
    else:
        alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
        blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
        clo_nd,chi_nd = _lohi(g_c,clo,chi)
        GA_Elem_divide_patch64(
                g_a, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                g_b, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data,
                g_c, <int64_t*>clo_nd.data, <int64_t*>chi_nd.data)

def elem_maximum(int g_a, int g_b, int g_c, alo=None, ahi=None, blo=None,
        bhi=None, clo=None, chi=None):
    """Computes the element-wise maximum of the two arrays.

    Arrays or array patches must be of the same types and same number of
    elements. For two-dimensional arrays::

        c(i,j) = max(a(i,j),b(i,j))

    If the data type is complex, then::

        c(i,j).real = max{ |a(i,j)|, |b(i,j)|} while c(i,j).image = 0

    The result (c) may replace one of the input arrays (a/b).

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle
        g_b : int
            the array handle
        g_c : int
            the array handle
        alo : 1D array-like of integers
            lower bound patch coordinates of g_a, inclusive
        ahi : 1D array-like of integers
            higher bound patch coordinates of g_a, exclusive
        blo : 1D array-like of integers
            lower bound patch coordinates of g_b, inclusive
        bhi : 1D array-like of integers
            higher bound patch coordinates of g_b, exclusive
        clo : 1D array-like of integers
            lower bound patch coordinates of g_c, inclusive
        chi : 1D array-like of integers
            higher bound patch coordinates of g_c, exclusive

    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    cdef np.ndarray[np.int64_t, ndim=1] clo_nd, chi_nd
    if (alo is None and ahi is None
            and blo is None and bhi is None
            and clo is None and chi is None):
        GA_Elem_maximum(g_a, g_b, g_c)
    else:
        alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
        blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
        clo_nd,chi_nd = _lohi(g_c,clo,chi)
        GA_Elem_maximum_patch64(
                g_a, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                g_b, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data,
                g_c, <int64_t*>clo_nd.data, <int64_t*>chi_nd.data)

def elem_minimum(int g_a, int g_b, int g_c, alo=None, ahi=None, blo=None,
        bhi=None, clo=None, chi=None):
    """Computes the element-wise minimum of the two arrays.

    Arrays or array patches must be of the same types and same number of
    elements. For two-dimensional arrays::

        c(i,j)  = min(a(i,j),b(i,j))

    If the data type is complex, then::

        c(i,j).real = min{ |a(i,j)|, |b(i,j)|} while c(i,j).image = 0

    The result (c) may replace one of the input arrays (a/b).

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle
        g_b : int
            the array handle
        g_c : int
            the array handle
        alo : 1D array-like of integers
            lower bound patch coordinates of g_a, inclusive
        ahi : 1D array-like of integers
            higher bound patch coordinates of g_a, exclusive
        blo : 1D array-like of integers
            lower bound patch coordinates of g_b, inclusive
        bhi : 1D array-like of integers
            higher bound patch coordinates of g_b, exclusive
        clo : 1D array-like of integers
            lower bound patch coordinates of g_c, inclusive
        chi : 1D array-like of integers
            higher bound patch coordinates of g_c, exclusive

    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    cdef np.ndarray[np.int64_t, ndim=1] clo_nd, chi_nd
    if (alo is None and ahi is None
            and blo is None and bhi is None
            and clo is None and chi is None):
        GA_Elem_minimum(g_a, g_b, g_c)
    else:
        alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
        blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
        clo_nd,chi_nd = _lohi(g_c,clo,chi)
        GA_Elem_minimum_patch64(
                g_a, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                g_b, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data,
                g_c, <int64_t*>clo_nd.data, <int64_t*>chi_nd.data)

def elem_multiply(int g_a, int g_b, int g_c, alo=None, ahi=None, blo=None,
        bhi=None, clo=None, chi=None):
    """Computes the element-wise product of the two arrays.

    Arrays or array patches must be of the same types and same number of
    elements. For two-dimensional arrays:

                c(i, j)  = a(i,j)*b(i,j)

    The result (c) may replace one of the input arrays (a/b).

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle
        g_b : int
            the array handle
        g_c : int
            the array handle
        alo : 1D array-like of integers
            lower bound patch coordinates of g_a, inclusive
        ahi : 1D array-like of integers
            higher bound patch coordinates of g_a, exclusive
        blo : 1D array-like of integers
            lower bound patch coordinates of g_b, inclusive
        bhi : 1D array-like of integers
            higher bound patch coordinates of g_b, exclusive
        clo : 1D array-like of integers
            lower bound patch coordinates of g_c, inclusive
        chi : 1D array-like of integers
            higher bound patch coordinates of g_c, exclusive

    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    cdef np.ndarray[np.int64_t, ndim=1] clo_nd, chi_nd
    if (alo is None and ahi is None
            and blo is None and bhi is None
            and clo is None and chi is None):
        GA_Elem_multiply(g_a, g_b, g_c)
    else:
        alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
        blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
        clo_nd,chi_nd = _lohi(g_c,clo,chi)
        GA_Elem_multiply_patch64(
                g_a, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                g_b, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data,
                g_c, <int64_t*>clo_nd.data, <int64_t*>chi_nd.data)

def error(char *message, int code=1):
    """Prints message and aborts safely with code."""
    GA_Error(message, code)

def fence():
    """Blocks the calling process until all the data transfers corresponding
    to GA operations called after ga.init_fence() complete.
    
    For example, since ga.put might return before the data reaches the final
    destination, ga_init_fence and ga_fence allow processes to wait until the
    data tranfer is fully completed:

        ga.init_fence()
        ga.put(g_a, ...)
        ga.fence()

    ga.fence() must be called after ga.init_fence(). A barrier, ga.sync(),
    assures completion of all data transfers and implicitly cancels all
    outstanding ga.init_fence() calls. ga.init_fence() and ga.fence() must be
    used in pairs, multiple calls to ga.fence() require the same number of
    corresponding ga.init_fence() calls. ga.init_fence()/ga_fence() pairs can
    be nested.

    ga.fence() works for multiple GA operations. For example:

        ga.init_fence()
        ga.put(g_a, ...)
        ga.scatter(g_a, ...)
        ga.put(g_b, ...)
        ga.fence()

    The calling process will be blocked until data movements initiated by two
    calls to ga_put and one ga_scatter complete.
    
    """
    GA_Fence()

def fill(int g_a, value, lo=None, hi=None):
    """Assign a single value to all elements in the array or patch.
    
    :Parameters:
        g_a : int
            the array handle
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : 1D array-like of integers
            higher bound patch coordinates, exclusive
    
    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd
    cdef int            ivalue
    cdef long           lvalue
    cdef long long      llvalue
    cdef float          fvalue
    cdef double         dvalue
    cdef long double    ldvalue
    cdef SingleComplex  fcvalue
    cdef DoubleComplex  dcvalue
    cdef void          *vvalue
    cdef int gtype=inquire_type(g_a)
    vvalue = _convert_multiplier(gtype, value,
            &ivalue,  &lvalue,  &llvalue,
            &fvalue,  &dvalue,  &ldvalue,
            &fcvalue, &dcvalue)
    if lo is None and hi is None:
        GA_Fill(g_a, vvalue)
    else:
        lo_nd,hi_nd = _lohi(g_a,lo,hi)
        NGA_Fill_patch64(
                g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data, vvalue)

def gather(int g_a, subsarray, np.ndarray values=None):
    """Gathers array elements from a global array into a local array.

    subsarray will be converted to an ndarray if it is not one already.  A
    two-dimensional array is allowed so long as its shape is (n,ndim) where n
    is the number of elements to gather and ndim is the number of dimensions
    of the target array.  Also, subsarray must be contiguous.

    For example, if the subsarray were two-dimensional::

        for k in range(n):
            v[k] = g_a[subsarray[k,0],subsarray[k,1],subsarray[k,2]...]

    For example, if the subsarray were one-dimensional::

        for k in range(n):
            base = n*ndim
            v[k] = g_a[subsarray[base+0],subsarray[base+1],subsarray[base+2]...]

    This is a one-sided operation. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] subsarray1_nd = None
    cdef np.ndarray[np.int64_t, ndim=2] subsarray2_nd = None
    cdef int gtype = inquire_type(g_a)
    cdef int ndim = GA_Ndim(g_a)
    cdef int64_t n
    # prepare subsarray
    try:
        subsarray1_nd = np.ascontiguousarray(subsarray, dtype=np.int64)
        n = len(subsarray1_nd) / ndim
    except ValueError:
        subsarray1_nd = None
        try:
            subsarray2_nd = np.ascontiguousarray(subsarray, dtype=np.int64)
            n = len(subsarray2_nd) # length of first dimension of subsarray2_nd
        except ValueError:
            raise ValueError, "subsarray must be either 1- or 2-dimensional"
    # prepare values array
    if values is None:
        values = np.ndarray(n, dtype=_to_dtype[gtype])
    else:
        if values.ndim != 1:
            raise ValueError, "values must be one-dimensional"
        if not values.flags['C_CONTIGUOUS']:
            raise ValueError, "values must be contiguous"
        if len(values) < n:
            raise ValueError, "values was not large enough"
    # call the wrapped function
    if subsarray1_nd is not None:
        NGA_Gather_flat64(g_a, <void*>values.data,
                <int64_t*>subsarray1_nd.data, n)
    elif subsarray2_nd is not None:
        NGA_Gather_flat64(g_a, <void*>values.data,
                <int64_t*>subsarray2_nd.data, n)
    else:
        raise ValueError, "how did this happen?"
    return values

def gemm(bint ta, bint tb, int64_t m, int64_t n, int64_t k,
        alpha, int g_a, int g_b, beta, int g_c):
    """Performs one of the matrix-matrix operations.
    
    C := alpha*op( A )*op( B ) + beta*C

    where op( X ) is one of op(X)=X or op(X) = X', alpha and beta are scalars,
    and A, B and C are matrices, with op(A) an m by k matrix, op(B) a k by n
    matrix, and C an m by n matrix.

    On entry, ta specifies the form of op( A ) to be used in the
    matrix multiplication as follows::

        ta = False, op(A) = A.
        ta = True,  op(A) = A'.

    This is a collective operation. 
    
    :Parameters:
        ta : bool
            transpose operator
        tb : bool
            transpose operator
        m : int
            number of rows of op(A) and of matrix C
        n : int
            number of columns of op(B) and of matrix C
        k : int
            number of columns of op(A) and rows of matrix op(B)
        alpha : object
            scale factor
        g_a : int
            handle to input array
        g_b : int
            handle to input array
        beta : object
            scale factor
        g_c : int
            handle to output array

    """
    cdef int gtype=inquire_type(g_a)
    #cdef int                 ialpha=1, ibeta=1
    #cdef long                lalpha=1, lbeta=1
    #cdef long long           llalpha=1, llbeta=1
    cdef float               falpha=1.0, fbeta=1.0
    cdef double              dalpha=1.0, dbeta=1.0
    #cdef long double         ldalpha=1.0, ldbeta=1.0
    cdef float complex       fcalpha=1.0, fcbeta=1.0
    cdef double complex      dcalpha=1.0, dcbeta=1.0
    #cdef long double complex ldalpha=1.0, ldbeta=1.0
    cdef SingleComplex       ga_fcalpha, ga_fcbeta
    cdef DoubleComplex       ga_dcalpha, ga_dcbeta
    cdef char ta_char = 'N'
    cdef char tb_char = 'N'
    if ta:
        ta_char = 'T'
    if tb:
        tb_char = 'T'
    if gtype == C_INT:
        raise TypeError, "C_INT not supported"
    elif gtype == C_LONG:
        raise TypeError, "C_LONG not supported"
    elif gtype == C_LONGLONG:
        raise TypeError, "C_LONGLONG not supported"
    elif gtype == C_FLOAT:
        falpha = alpha
        fbeta = beta
        GA_Sgemm64(ta_char, tb_char, m, n, k, falpha, g_a, g_b, fbeta, g_c)
    elif gtype == C_DBL:
        dalpha = alpha
        dbeta = beta
        GA_Dgemm64(ta_char, tb_char, m, n, k, dalpha, g_a, g_b, dbeta, g_c)
    elif gtype == C_LDBL:
        raise TypeError, "C_LDBL not supported"
    elif gtype == C_SCPL:
        fcalpha = alpha
        fcbeta = beta
        ga_fcalpha.real = fcalpha.real
        ga_fcalpha.imag = fcalpha.imag
        ga_fcbeta.real = fcbeta.real
        ga_fcbeta.imag = fcbeta.imag
        GA_Cgemm64(ta_char, tb_char, m, n, k, ga_fcalpha, g_a, g_b, ga_fcbeta, g_c)
    elif gtype == C_DCPL:
        dcalpha = alpha
        dcbeta = beta
        ga_dcalpha.real = dcalpha.real
        ga_dcalpha.imag = dcalpha.imag
        ga_dcbeta.real = dcbeta.real
        ga_dcbeta.imag = dcbeta.imag
        GA_Zgemm64(ta_char, tb_char, m, n, k, ga_dcalpha, g_a, g_b, ga_dcbeta, g_c)
    elif gtype == C_LDCPL:
        raise TypeError, "C_LDCPL not supported (yet)"
    else:
        raise TypeError

def get(int g_a, lo=None, hi=None, np.ndarray buffer=None):
    """Copies data from global array section to the local array buffer.
    
    The local array is assumed to be have the same number of dimensions as the
    global array. Any detected inconsitencies/errors in the input arguments
    are fatal.

    This is a one-sided operation.

    :Parameters:
        g_a : int
            the array handle
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : 1D array-like of integers
            higher bound patch coordinates, exclusive
        buffer : 
            an ndarray of the appropriate type, large enough to hold lo,hi

    :returns: The local array buffer.
    
    """
    return _get_common(g_a, lo, hi, buffer)

cdef _get_common(int g_a, lo=None, hi=None, np.ndarray buffer=None,
        bint nb=False, bint periodic=False, skip=None):
    """Copies data from global array section to the local array buffer.
    
    The local array is assumed to have the same shape as the requested region,
    or the local array can be 1-dimensional so long as it has the same number
    of elements as the requested region. Any detected inconsitencies raise a
    ValueError.

    This is a one-sided operation.

    :Parameters:
        g_a : int
            the array handle
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : 1D array-like of integers
            higher bound patch coordinates, exclusive
        buffer : ndarray
            should either be 1D and len(buffer)==np.prod(hi-lo), or
            np.all(buffer.shape == hi-lo) i.e. buffer is 1D and same size as
            requested region or buffer is the same shape as requested region
        nb : bool
            whether this call is non-blocking (see ga.nbget)
        periodic : bool
            whether this call is periodic (see ga.periodic_get)
        skip : 1D array-like of integers
            strides for each dimension

    :returns: The local array buffer (and the nonblocking handle if nb=True.)
    
    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd, ld_nd, shape, skip_nd
    cdef int gtype=inquire_type(g_a)
    cdef ga_nbhdl_t nbhandle
    dtype = _to_dtype[gtype]
    lo_nd,hi_nd = _lohi(g_a,lo,hi)
    shape = hi_nd-lo_nd+1
    if skip is None:
        skip_nd = None
    else:
        skip_nd = _inta64(skip)
        shape = (hi_nd-lo_nd)/skip_nd+1
    if buffer is None:
        buffer = np.ndarray(shape, dtype=dtype)
    elif buffer.dtype != dtype:
        raise ValueError, "buffer is wrong type :: buffer=%s != %s" % (
                buffer.dtype, dtype)
    # Due to GA restrictions, buffer must not have negative strides
    # and buffer's last stride must be same as itemsize
    strides = [buffer.strides[i]/buffer.itemsize for i in range(buffer.ndim)]
    if strides[-1] != 1:
        raise ValueError, "first dimension of buffer cannot be strided"
    if np.any(np.asarray(strides) < 0):
        raise ValueError, "buffer cannot have negative strides"
    # we allow 1-d "flat" buffers in addition to buffers matching the shape of
    # requested region
    if buffer.ndim == 1:
        if buffer.size != np.prod(shape):
            raise ValueError, ('buffer size does not match shape :: '
                    'buffer.size=%s != np.prod(shape)=%s' % (
                    buffer.size, np.prod(shape)))
        ld_nd = shape[1:]
    else:
        buffer_shape = [buffer.shape[i] for i in range(buffer.ndim)]
        if not np.all(buffer_shape == shape):
            raise ValueError, ('buffer shape does not match request shape :: '
                    'buffer_shape=%s != shape=%s' % (
                    buffer_shape, shape))
        ld_nd = np.asarray([strides[i]/strides[i+1]
                for i in range(buffer.ndim-1)], dtype=np.int64)
    if nb:
        NGA_NbGet64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <void*>buffer.data, <int64_t*>ld_nd.data, &nbhandle)
        return buffer,nbhandle
    elif periodic:
        NGA_Periodic_get64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <void*>buffer.data, <int64_t*>ld_nd.data)
        return buffer
    elif skip is not None:
        NGA_Strided_get64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <int64_t*>skip_nd.data,
                <void*>buffer.data, <int64_t*>ld_nd.data)
        return buffer
    else:
        NGA_Get64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <void*>buffer.data, <int64_t*>ld_nd.data)
        return buffer

def get_block_info(int g_a):
    """Returns information about the block-cyclic distribution.

    The number of blocks along each of the array axes are returned in the
    array num_blocks and the dimensions of the individual blocks, specified in
    the ga.set_block_cyclic or ga.set_block_cyclic_proc_grid subroutines, are
    returned in block_dims.

    This is a local function.

    :Parameters:
        g_a : int
            the array handle

    :returns: The number of blocks along each of the array axes and the dimensions of thet individual blocks, in that order, as ndarrays.

    """
    cdef np.ndarray[np.int32_t, ndim=1] num_blocks, block_dims
    cdef int ndim = GA_Ndim(g_a)
    num_blocks = np.zeros(ndim, dtype=np.int32)
    block_dims = np.zeros(ndim, dtype=np.int32)
    GA_Get_block_info(g_a, <int*>num_blocks.data, <int*>block_dims.data)
    return num_blocks,block_dims

def get_diag(int g_a, int g_v):
    """Inserts the diagonal elements of this matrix g_a into the vector g_v.

    This is a collective operation.
    
    :Parameters:
        g_a : int
            the array handle

    """
    GA_Get_diag(g_a, g_v)

def get_debug():
    """Returns the value of an internal flag in the GA library whose value can
    be set using the ga.set_debug() subroutine.

    This is a local operation.

    """
    return GA_Get_debug()

def gop(X, char *op):
    """Global operation.

    X(1:N) is a vector present on each process. gop 'sums' elements of X
    accross all nodes using the commutative operator op. The result is
    broadcast to all nodes. Supported operations include '+', '*', 'max',
    'min', 'absmax', 'absmin'. The use of lowerecase for operators is
    necessary.

    X must be a contiguous array-like.  X is not guaranteed to be modified
    in-place so use as:

    >>> value = ga.gop((1,2,3), "+")

    This operation is provided only for convenience purposes: it is available
    regardless of the message-passing library that GA is running with.

    This is a collective operation. 

    """
    cdef np.ndarray X_nd = np.asarray(X)
    cdef int size = 0
    if not X_nd.flags['C_CONTIGUOUS']:
        raise ValueError, "X must be contiguous"
    try:
        size = len(X_nd)
    except TypeError:
        size = 1
    if X_nd.dtype == np.dtype(np.intc):
        GA_Igop(<int*>X_nd.data, size, op)
    elif X_nd.dtype == np.dtype(np.long):
        GA_Lgop(<long*>X_nd.data, size, op)
    elif X_nd.dtype == np.dtype(np.longlong):
        GA_Llgop(<long long*>X_nd.data, size, op)
    elif X_nd.dtype == np.dtype(np.single):
        GA_Fgop(<float*>X_nd.data, size, op)
    elif X_nd.dtype == np.dtype(np.double):
        GA_Dgop(<double*>X_nd.data, size, op)
    elif X_nd.dtype == np.dtype(np.complex64):
        GA_Cgop(<SingleComplex*>X_nd.data, size, op)
    elif X_nd.dtype == np.dtype(np.complex128):
        GA_Zgop(<DoubleComplex*>X_nd.data, size, op)
    else:
        raise TypeError, "type not supported by ga.gop %s" % X_nd.dtype
    return X_nd

def gop_add(X):
    return gop(X, "+")

def gop_multiply(X):
    return gop(X, "*")

def gop_max(X):
    return gop(X, "max")

def gop_min(X):
    return gop(X, "min")

def gop_absmax(X):
    return gop(X, "absmax")

def gop_absmin(X):
    return gop(X, "absmin")

def has_ghosts(int g_a):
    """Determines whether any dimension of the given array has ghost cells.

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle

    :returns: True if the global array has some dimensions for which the ghost cell width is greater than zero, it returns False otherwise.

    """
    if GA_Has_ghosts(g_a) == 1:
        return True
    return False

def init_fence():
    """Initializes tracing of completion status of data movement operations.

    This operation is local.

    """
    GA_Init_fence()

def initialize():
    """Allocates and initializes internal data structures in Global Arrays.

    This is a collective operation.

    """
    import atexit
    global _initialized
    GA_Initialize()
    GA_Register_stack_memory(_gapy_malloc, _gapy_free)
    atexit.register(terminate)
    _initialized = True

def initialize_ltd(size_t limit):
    """Allocates and initializes internal data structures and sets limit for
    memory used in global arrays.
    
    The limit is per process: it is the amount of memory that the given
    processor can contribute to collective allocation of global arrays. It
    does not include temporary storage that GA might be allocating (and
    releasing) during execution of a particular operation.

    limit = 0 means "allow unlimited memory usage" in which case this
    operation is equivalent to GA_initialize.

    This is a collective operation. 

    """
    global _initialized
    GA_Initialize_ltd(limit)
    _initialized = True

def initialized():
    """Returns whether ga has been initialized."""
    return _initialized

cpdef inquire(int g_a):
    cdef int gtype
    cdef int ndim = GA_Ndim(g_a)
    cdef np.ndarray[np.int64_t, ndim=1] dims=np.zeros((ndim), dtype=np.int64)
    NGA_Inquire64(g_a, &gtype, &ndim, <int64_t*>dims.data)
    return gtype,dims

cpdef np.ndarray[np.int64_t, ndim=1] inquire_dims(int g_a):
    cdef int gtype
    cdef np.ndarray[np.int64_t, ndim=1] dims
    gtype,dims = inquire(g_a)
    return dims

def inquire_memory():
    """Returns amount of memory (in bytes) used in the allocated global arrays
    on the calling processor.

    This operation is local. 

    """
    return GA_Inquire_memory()

def inquire_name(int g_a):
    """Returns the name of an array represented by the handle g_a.

    This operation is local.

    :Parameters:
        g_a : int
            the array handle

    """
    return GA_Inquire_name(g_a)

cpdef int inquire_type(int g_a):
    cdef int gtype
    cdef np.ndarray[np.int64_t, ndim=1] dims
    gtype,dims = inquire(g_a)
    return gtype

def is_mirrored(int g_a):
    """Checks whether the array is mirrored.
    
    This is a  local  operation. 

    :returns: True if it is a mirrored array, else returns False.

    """
    if GA_Is_mirrored(g_a) == 1:
        return True
    return False

def llt_solve(int g_a, int g_b):
    """Solves a system of linear equations

        A * X = B

    using the Cholesky factorization of an NxN double precision symmetric
    positive definite matrix A (epresented by handle g_a). On successful exit
    B will contain the solution X.

    This is a collective operation. 

    :Parameters:
        g_a : int
            the coefficient matrix
        g_b : int
            the rhs/solution matrix

    :returns: 0 if successful; >0 if the leading minor of this order is not positive definite and the factorization could not be completed

    """
    return GA_Llt_solve(g_a, g_b)

def locate(int g_a, subscript):
    """Return the GA compute process id that 'owns' the data.
    
    If any element of subscript[] is out of bounds "-1" is returned.

    This operation is local.

    :Parameters:
        g_a : int
            the array handle
        subscript : 1D array-like of integers
            len(subscript) should be ndim

    """
    cdef np.ndarray[np.int64_t, ndim=1] subscript_nd
    subscript_nd = _inta64(subscript)
    return NGA_Locate64(g_a, <int64_t*>subscript_nd.data)

def locate_nnodes(int g_a, lo, hi):
    """Return the number of process which own the specified patch.

    This operation is local.

    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd,hi_nd
    cdef int result
    lo_nd,hi_nd = _lohi(g_a,lo,hi)
    return NGA_Locate_nnodes64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data)

def locate_region(int g_a, lo, hi):
    """Return the list of the GA processes id that 'own' the data.
    
    Parts of the specified patch might be actually 'owned' by several
    processes. If lo/hi are out of bounds "0" is returned, otherwise return
    value is equal to the number of processes that hold the data .
          
    map[i][0] - lo[ndim]
    map[i][1] - hi[ndim]
    procs[i]  - processor id that owns data in patch lo[i]:hi[i]

    This operation is local. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd
    cdef np.ndarray[np.int64_t, ndim=1] hi_nd
    cdef np.ndarray[np.int64_t, ndim=1] map
    cdef np.ndarray[np.int64_t, ndim=3] map_reshape
    cdef np.ndarray[np.int32_t, ndim=1] procs
    cdef int np_result
    cdef int np_guess
    cdef int ndim = GA_Ndim(g_a)
    lo_nd,hi_nd = _lohi(g_a,lo,hi)
    np_guess = NGA_Locate_nnodes64(
            g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data)
    map = np.ndarray(np_guess*ndim*2, dtype=np.int64)
    procs = np.ndarray(np_guess, dtype=np.int32)
    np_result = NGA_Locate_region64(g_a,
            <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
            <int64_t*>map.data, <int*>procs.data)
    assert(np_guess == np_result)
    # TODO then slice it and reshape to something useful?
    map_reshape = map.reshape(np_result,2,ndim)
    # need to add 1 to every 'hi' value
    map_reshape[:,1,:] += 1
    return map_reshape,procs

def lock(int mutex):
    """Locks a mutex object identified by the mutex number.
    
    It is a fatal error for a process to attempt to lock a mutex which was
    already locked by this process. 

    """
    GA_Lock(mutex)

def lu_solve(int g_a, int g_b, bint trans=False):
    """Solve the system of linear equations op(A)X = B based on the LU
    factorization.

    op(A) = A or A' depending on the parameter trans
    trans = False means that the transpose operator should not be applied.
    trans = True means that the transpose operator should be applied.

    Matrix A is a general real matrix. Matrix B contains possibly multiple rhs
    vectors. The array associated with the handle g_b is overwritten by the
    solution matrix X.

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle for the coefficient matrix
        g_b : int
            the array handle for the solution matrix
        trans : bool
            transpose (True) or not transpose (False)

    """
    cdef char ctrans = 'N'
    if trans:
        ctrans = 'T'
    GA_Lu_solve(ctrans, g_a, g_b)

def mask_sync(bint first, bint last):
    """This subroutine can be used to remove synchronization calls from around
    collective operations.
    
    Setting the parameter first=False removes the synchronization prior to the
    collective operation, setting last=False removes the synchronization call
    after the collective operation. This call is applicable to all collective
    operations. It most be invoked before each collective operation.
    
    This is a collective operation. 

    :Parameters:
        first : bool
            mask for prior internal synchronization
        last : bool
            mask for post internal synchronization

    """
    GA_Mask_sync(first,last)

def matmul_patch(bint transa, bint transb, alpha, beta,
        int g_a, alo, ahi,
        int g_b, blo, bhi,
        int g_c, clo, chi):
    """An n-dimensional patch version of ga_dgemm.

    C[clo[]:chi[]] := alpha * AA[alo[]:ahi[]] *
                                   BB[blo[]:bhi[]] ) + beta*C[clo[]:chi[]],

    where AA = op(A), BB = op(B), and op( X ) is one of

    op( X ) = X   or   op( X ) = X',

    It works for both double and DoubleComplex data tape.

    This is a collective operation. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    cdef np.ndarray[np.int64_t, ndim=1] clo_nd, chi_nd
    cdef int gtype=inquire_type(g_a)
    cdef int            ialpha,  ibeta
    cdef long           lalpha,  lbeta
    cdef long long      llalpha, llbeta
    cdef float          falpha,  fbeta
    cdef double         dalpha,  dbeta
    cdef long double    ldalpha, ldbeta
    cdef SingleComplex  fcalpha, fcbeta
    cdef DoubleComplex  dcalpha, dcbeta
    cdef void          *valpha, *vbeta
    cdef char char_transa = 'N'
    cdef char char_transb = 'N'
    if alpha is None:
        alpha = 1
    valpha = _convert_multiplier(gtype, alpha,
            &ialpha,  &lalpha,  &llalpha,
            &falpha,  &dalpha,  &ldalpha,
            &fcalpha, &dcalpha)
    if beta is None:
        beta = 1
    vbeta = _convert_multiplier(gtype, beta,
            &ibeta,  &lbeta,  &llbeta,
            &fbeta,  &dbeta,  &ldbeta,
            &fcbeta, &dcbeta)
    alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
    blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
    clo_nd,chi_nd = _lohi(g_c,clo,chi)
    if transa:
        char_transa = 'T'
    if transb:
        char_transb = 'T'
    NGA_Matmul_patch64(char_transa, char_transb, valpha, vbeta,
            g_a, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
            g_b, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data,
            g_c, <int64_t*>clo_nd.data, <int64_t*>chi_nd.data)

def median(int g_a, int g_b, int g_c, int g_m,
        alo=None, ahi=None, blo=None, bhi=None,
        clo=None, chi=None, mlo=None, mhi=None):
    """Computes the componentwise Median of three arrays or patches g_a, g_b,
    and g_c, and stores the result in this array or patch g_m.
    
    The result (m) may replace one of the input arrays (a/b/c).

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle
        g_b : int
            the array handle
        g_c : int
            the array handle
        g_m : int
            the array handle for the result
        alo : 1D array-like of integers
            lower bound patch coordinates of g_a, inclusive
        ahi : 1D array-like of integers
            higher bound patch coordinates of g_a, exclusive
        blo : 1D array-like of integers
            lower bound patch coordinates of g_b, inclusive
        bhi : 1D array-like of integers
            higher bound patch coordinates of g_b, exclusive
        clo : 1D array-like of integers
            lower bound patch coordinates of g_c, inclusive
        chi : 1D array-like of integers
            higher bound patch coordinates of g_c, exclusive
        mlo : 1D array-like of integers
            lower bound patch coordinates of g_m, inclusive
        mhi : 1D array-like of integers
            higher bound patch coordinates of g_m, exclusive

    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    cdef np.ndarray[np.int64_t, ndim=1] clo_nd, chi_nd
    cdef np.ndarray[np.int64_t, ndim=1] mlo_nd, mhi_nd
    if (alo is None and ahi is None
            and blo is None and bhi is None
            and clo is None and chi is None
            and mlo is None and mhi is None):
        GA_Median(g_a, g_b, g_c, g_m)
    else:
        alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
        blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
        clo_nd,chi_nd = _lohi(g_c,clo,chi)
        mlo_nd,mhi_nd = _lohi(g_m,mlo,mhi)
        GA_Median_patch64(
                g_a, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                g_b, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data,
                g_c, <int64_t*>clo_nd.data, <int64_t*>chi_nd.data,
                g_m, <int64_t*>mlo_nd.data, <int64_t*>mhi_nd.data)

def memory_avail():
    """Returns amount of memory (in bytes) left for allocation of new global
    arrays on the calling processor.

    Note: If ga.uses_ma() returns True, then ga.memory_avail() returns the
    lesser of the amount available under the GA limit and the amount available
    from MA (according to ma.inquire_avail() operation). If no GA limit has
    been set, it returns what MA says is available.

    If ( ! ga.uses_ma() && ! ga.memory_limited() ) returns < 0, indicating
    that the bound on currently available memory cannot be determined.

    This operation is local. 
    
    """
    return GA_Memory_avail()

def memory_limited():
    """Indicates if limit is set on memory usage in Global Arrays on the
    calling processor.
    
    This operation is local. 

    :returns: True for "yes", False for "no"

    """
    if 1 == GA_Memory_limited():
        return True
    return False

def merge_distr_patch(int g_a, alo, ahi, int g_b, blo, bhi):
    """This function merges all copies of a patch of a mirrored array (g_a)
    into a patch in a distributed array (g_b).

    This is a  collective  operation. 

    :Parameters:
        g_a : int
            array handle
        alo : 1D array-like of integers
            g_a patch coordinate
        ahi : 1D array-like of integers
            g_a patch coordinate
        g_b : int
            array handle
        blo : 1D array-like of integers
            g_b patch coordinate
        bhi : 1D array-like of integers
            g_b patch coordinate

    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
    blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
    NGA_Merge_distr_patch64(
            g_a, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
            g_b, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data)

def merge_mirrored(int g_a):
    """This subroutine merges mirrored arrays by adding the contents of each
    array across nodes.
    
    The result is that the each mirrored copy of the array represented by g_a
    is the sum of the individual arrays before the merge operation. After the
    merge, all mirrored arrays are equal.

    This is a  collective  operation. 

    :Parameters:
        g_a : int
            array handle

    """
    GA_Merge_mirrored(g_a)

def nbacc(int g_a, buffer, lo=None, hi=None, alpha=None):
    """Non-blocking version of ga.acc.

    The accumulate operation can be completed locally by making a call to the
    ga.nbwait() routine.
    
    Combines data from buffer with data in the global array patch.
    
    The buffer array is assumed to be have the same number of
    dimensions as the global array.  If the buffer is not contiguous, a
    contiguous copy will be made.
    
    global array section (lo[],hi[]) += alpha * buffer

    This is a non-blocking and one-sided and atomic operation.

    :Parameters:
        g_a : int
            the array handle
        buffer : array-like
            must be contiguous and have same number of elements as patch
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : 1D array-like of integers
            higher bound patch coordinates, exclusive
        alpha : object
            multiplier (converted to the appropriate type)

    :returns: The non-blocking request handle.

    """
    return _acc_common(g_a, buffer, lo, hi, alpha, True)

def nbget(int g_a, lo=None, hi=None, np.ndarray buffer=None):
    """Non-blocking version of the blocking ga.get operation.
    
    The get operation can be completed locally by making a call to the
    ga.nbwait() routine.
    
    Copies data from global array section to the local array buffer.
    
    The local array is assumed to be have the same number of dimensions as the
    global array. Any detected inconsitencies/errors in the input arguments
    are fatal.

    This is a non-blocking and one-sided operation.

    :Parameters:
        g_a : int
            the array handle
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : 1D array-like of integers
            higher bound patch coordinates, exclusive
        buffer : ndarray
            Fill this buffer instead of allocating a new one internally.
            Must be contiguous and have same number of elements as patch.

    :returns: The local array buffer.
    
    """
    return _get_common(g_a, lo, hi, buffer, True)

def nblock(int g_a):
    """Returns the number of partitions of each array dimension for g_a.

    This operation is local. 

    :Parameters:
        g_a : int
            array handle

    """
    cdef np.ndarray[np.int32_t, ndim=1] nblock_nd
    cdef int ndim = GA_Ndim(g_a)
    nblock_nd = np.zeros(ndim, dtype=np.intc)
    GA_Nblock(g_a, <int*>nblock_nd.data)
    return nblock_nd

def nbput(int g_a, buffer, lo=None, hi=None):
    """Non-blocking version of the blocking put operation.
    
    The put operation can be completed locally by making a call to the
    ga.nbwait() routine.
    
    Copies data from local array buffer to the global array section.
    
    The local array is assumed to be have the same number of dimensions as the
    global array.  Any detected inconsitencies/errors in input arguments are
    fatal.

    This is a one-sided operation. 

    :Parameters:
        g_a : int
            the array handle
        buffer : array-like
            the data to put
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : 1D array-like of integers
            higher bound patch coordinates, exclusive

    """
    return _put_common(g_a, buffer, lo, hi, True, False)
 
def nbwait(ga_nbhdl_t nbhandle):
    """This function completes a non-blocking one-sided operation locally.

    Waiting on a nonblocking put or an accumulate operation assures that data
    was injected into the network and the user buffer can be now be reused.
    Completing a get operation assures data has arrived into the user memory
    and is ready for use. Wait operation ensures only local completion. Unlike
    their blocking counterparts, the nonblocking operations are not ordered
    with respect to the destination. Performance being one reason, the other
    reason is that by ensuring ordering we incur additional and possibly
    unnecessary overhead on applications that do not require their operations
    to be ordered. For cases where ordering is necessary, it can be done by
    calling a fence operation. The fence operation is provided to the user to
    confirm remote completion if needed. 

    """
    NGA_NbWait(&nbhandle)

def ndim(int g_a):
    """Returns the number of dimensions in array represented by the handle g_a.

    This operation is local.
    
    :Parameters:
        g_a : int
            the array handle

    :returns: the number of dimensions in the array g_a

    """
    return GA_Ndim(g_a)

def nnodes():
    """Returns the number of the GA compute (user) processes.

    This operation is local.
    
    :returns: the number of GA compute (user) processes

    """
    return GA_Nnodes()

def nodeid():
    """Returns the GA process id (0, ..., ga.nnodes()-1) of the requesting
    compute process.

    This operation is local.
    
    :returns: the GA process id

    """
    return GA_Nodeid()

def norm1(int g_a):
    """Computes the 1-norm of the matrix or vector g_a.

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle

    :returns: the 1-norm of the matrix or vector g_a (as a float)

    """
    cdef double nm
    GA_Norm1(g_a, &nm)
    return nm

def norm_infinity(int g_a):
    """Computes the 1-norm of the matrix or vector g_a.

    This is a collective operation. 

    :Parameters:
        g_a : int
            the array handle

    :returns: the 1-norm of the matrix or vector g_a

    """
    cdef double nm
    GA_Norm_infinity(g_a, &nm)
    return nm

def enum(int g_a, lo=None, hi=None, start=None, inc=None):
    """This subroutine enumerates the values of an array between elements lo
    and hi starting with the value istart and incrementing each subsequent
    value by inc.
    
    This operation is only applicable to 1-dimensional arrays.

    An example of its use is shown below:

        ga.enum(g_a, 1, n, 7, 2)
        # g_a: 7  9 11 13 15 17 19 21 23 ...

    This is a collective operation.

    :Parameters:
        g_a : int
            the array handle
        lo : 1D array-like of integers
            patch coordinate
        hi : 1D array-like of integers
            patch coordinate
        start : object
            starting value of enumeration (converted to appropriate type)
        inc : object
            increment value (converted to appropriate type)

    """
    cdef np.ndarray[np.int64_t, ndim=1] hi_nd = inquire_dims(g_a)-1
    cdef int64_t c_lo=0, c_hi=hi_nd[0]
    cdef int gtype=inquire_type(g_a)
    cdef int            istart, iinc
    cdef long           lstart, linc
    cdef long long      llstart, llinc
    cdef float          fstart, finc
    cdef double         dstart, dinc
    cdef long double    ldstart, ldinc
    cdef SingleComplex  fcstart, fcinc
    cdef DoubleComplex  dcstart, dcinc
    cdef void          *vstart=NULL, *vinc=NULL
    if lo is not None:
        c_lo = lo
    if hi is not None:
        c_hi = hi
    if start is None:
        start = 1
    if inc is None:
        inc = 1
    vstart = _convert_multiplier(gtype, start,
            &istart,  &lstart,  &llstart,
            &fstart,  &dstart,  &ldstart,
            &fcstart, &dcstart)
    vinc = _convert_multiplier(gtype, inc,
            &iinc,  &linc,  &llinc,
            &finc,  &dinc,  &ldinc,
            &fcinc, &dcinc)
    GA_Patch_enum64(g_a, c_lo, c_hi, vstart, vinc)

def pack(int g_src, int g_dst, int g_msk, lo=None, hi=None):
    """The pack subroutine is designed to compress the values in the source
    vector g_src into a smaller destination array g_dst based on the values
    in an integer mask array g_msk. The values lo and hi denote the range of
    elements that should be compressed and the number of values placed in the
    compressed array is returned.  This operation is the complement of the
    ga.unpack operation. An example is shown below::

        icount = ga.pack(g_src, g_dst, g_msk, 1, n);
        # g_msk:   1  0  0  0  0  0  1  0  1  0  0  1  0  0  1  1  0
        # g_src:   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
        # g_dst:   1  7  9 12 15 16
        # icount:  6

    This is a collective operation.

    :Parameters:
        g_src : int
            handle for source arrray
        g_dst : int
            handle for destination array
        g_msk : int
            handle for integer array representing mask
        lo : 1D array-like of integers
            low value of range on which operation is performed
        hi : 1D array-like of integers
            hi value of range on which operation is performed

    """
    cdef np.ndarray[np.int64_t, ndim=1] hi_nd = inquire_dims(g_src)-1
    cdef int64_t c_lo=0, c_hi=hi_nd[0], icount
    if lo is not None:
        c_lo = lo
    if hi is not None:
        c_hi = hi
    GA_Pack64(g_src, g_dst, g_msk, lo, hi, &icount)
    return icount

def periodic_acc(int g_a, buffer, lo=None, hi=None, alpha=None):
    """Periodic version of ga.acc.

    The indices can extend beyond the array boundary/dimensions in which case
    the libray wraps them around.
    
    Combines data from buffer with data in the global array patch.
    
    The buffer array is assumed to be have the same number of
    dimensions as the global array.  If the buffer is not contiguous, a
    contiguous copy will be made.
    
        global array section (lo[],hi[]) += alpha * buffer

    This is a one-sided and atomic operation.

    :Parameters:
        g_a : int
            the array handle
        buffer : array-like
            must be contiguous and have same number of elements as patch
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : array-like of integers
            higher bound patch coordinates, exclusive
        alpha : object
            multiplier (converted to the appropriate type)

    """
    _acc_common(g_a, buffer, lo, hi, alpha, False, True)

def periodic_get(int g_a, lo=None, hi=None, np.ndarray buffer=None):
    """Periodic version of ga.get.
    
    The indices can extend beyond the array boundary/dimensions in which case
    the libray wraps them around.
    
    Copies data from global array section to the local array buffer.
    
    The local array is assumed to be have the same number of dimensions as the
    global array. Any detected inconsitencies/errors in the input arguments
    are fatal.

    This is a one-sided operation.

    :Parameters:
        g_a : int
            the array handle
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : array-like of integers
            higher bound patch coordinates, exclusive
        buffer : array-like
            must be contiguous and have same number of elements as patch

    :returns: The local array buffer.
    
    """
    return _get_common(g_a, lo, hi, buffer, False, True)

def periodic_put(int g_a, buffer, lo=None, hi=None):
    """Periodic version of ga.put.

    The indices can extend beyond the array boundary/dimensions in which case
    the libray wraps them around.
    
    Copies data from local array buffer to the global array section.
    
    The local array is assumed to be have the same number of dimensions as the
    global array.  Any detected inconsitencies/errors in input arguments are
    fatal.

    This is a one-sided operation. 

    :Parameters:
        g_a : int
            the array handle
        buffer : array-like
            the data to put
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : array-like of integers
            higher bound patch coordinates, exclusive

    """
    _put_common(g_a, buffer, lo, hi, False, True)

def pgroup_absolute_id(int pgroup, int pid):
    """TODO"""
    return GA_Pgroup_absolute_id(pgroup, pid)

def pgroup_brdcst(int pgroup, buffer, int root=0):
    """Broadcast from process root to all other processes in the same group.

    If the buffer is not contiguous, an error is raised.  This operation is
    provided only for convenience purposes: it is available regardless of the
    message-passing library that GA is running with.

    This is a collective operation. 

    :Parameters:
        pgroup : int
            processor group handle
        buffer : array-like
            the message
        root : int
            the process which is sending

    :returns: The buffer in case a temporary was passed in.

    """
    cdef np.ndarray buffer_nd
    buffer_nd = np.asarray(buffer)
    if not buffer_nd.flags['C_CONTIGUOUS']:
        raise ValueError, "the buffer must be contiguous"
    #if buffer_nd.ndim != 1:
    #    raise ValueError, "the buffer must be one-dimensional"
    GA_Pgroup_brdcst(pgroup, buffer_nd.data,
            buffer_nd.size*buffer_nd.itemsize, root)
    return buffer_nd

def pgroup_create(list):
    """Creates a processor group.
    
    At present, it must be invoked by all processors in the current default
    processor group. The list of processors use the indexing scheme of the
    default processor group.  If the default processor group is the world
    group, then these indices are the usual processor indices. This function
    returns a process group handle that can be used to reference this group by
    other functions.

    This is a collective operation on the default processor group. 

    """
    cdef np.ndarray[np.int32_t, ndim=1] list_nd
    list_nd = _inta32(list)
    return GA_Pgroup_create(<int*>list_nd.data, len(list_nd))

def pgroup_destroy(int pgroup):
    """Frees up a processor group handle.
    
    This is a collective operation on the default processor group.

    :returns: True if the handle was previously active.  False if the handle was not previously active.

    """
    if 0 == GA_Pgroup_destroy(pgroup):
        return False
    return True

def pgroup_get_default():
    """Returns a handle to the default processor group.

    The return value can then be used to create a global array using one of
    the ga.create or ga.set_pgroup calls.

    This is a local operation.  

    """
    return GA_Pgroup_get_default()

def pgroup_get_mirror():
    """Returns a handle to the mirrored processor group.

    The return value can then be used to create a global array using one of
    the ga.create or ga.set_pgroup calls.

    This is a local operation.  

    """
    return GA_Pgroup_get_mirror()

def pgroup_get_world():
    """Returns a handle to the world processor group.
    
    The return value can then be used to create a global array using one of
    the ga.create or ga.set_pgroup calls.

    This is a local operation. 

    """
    return GA_Pgroup_get_world()

def pgroup_gop(int pgroup, X, char *op):
    """Global operation.

    X(1:N) is a vector present on each process in the group. gop 'sums'
    elements of X accross all nodes using the commutative operator op. The
    result is broadcast to all nodes. Supported operations include '+', '*',
    'max', 'min', 'absmax', 'absmin'. The use of lowerecase for operators is
    necessary.

    X must be a contiguous array-like.  X is not guaranteed to be modified
    in-place so use as:

    >>> value = ga.gop((1,2,3), "+")

    This operation is provided only for convenience purposes: it is available
    regardless of the message-passing library that GA is running with.

    This is a collective operation. 

    """
    cdef np.ndarray X_nd = np.asarray(X)
    if not X_nd.flags['C_CONTIGUOUS']:
        raise ValueError, "X must be contiguous"
    if X_nd.dtype == np.intc:
        GA_Pgroup_igop(pgroup, <int*>X_nd.data, len(X_nd), op)
    elif X_nd.dtype == np.long:
        GA_Pgroup_lgop(pgroup, <long*>X_nd.data, len(X_nd), op)
    elif X_nd.dtype == np.longlong:
        GA_Pgroup_llgop(pgroup, <long long*>X_nd.data, len(X_nd), op)
    elif X_nd.dtype == np.single:
        GA_Pgroup_fgop(pgroup, <float*>X_nd.data, len(X_nd), op)
    elif X_nd.dtype == np.double:
        GA_Pgroup_dgop(pgroup, <double*>X_nd.data, len(X_nd), op)
    elif X_nd.dtype == np.complex64:
        GA_Pgroup_cgop(pgroup, <SingleComplex*>X_nd.data, len(X_nd), op)
    elif X_nd.dtype == np.complex128:
        GA_Pgroup_zgop(pgroup, <DoubleComplex*>X_nd.data, len(X_nd), op)
    else:
        raise TypeError, "type not supported by ga.pgroup_gop %s" % X_nd.dtype
    return X_nd

def pgroup_gop_add(int pgroup, X):
    return pgroup_gop(pgroup, X, "+")

def pgroup_gop_multiply(int pgroup, X):
    return pgroup_gop(pgroup, X, "*")

def pgroup_gop_max(int pgroup, X):
    return pgroup_gop(pgroup, X, "max")

def pgroup_gop_min(int pgroup, X):
    return pgroup_gop(pgroup, X, "min")

def pgroup_gop_absmax(int pgroup, X):
    return pgroup_gop(pgroup, X, "absmax")

def pgroup_gop_absmin(int pgroup, X):
    return pgroup_gop(pgroup, X, "absmin")

def pgroup_nnodes(int pgroup=-1):
    """Returns the number of processors contained in the group specified by
    pgroup.

    This is a local local operation. 

    :Parameters:
        pgroup : int
            the group handle

    """
    if pgroup < 0:
        pgroup = pgroup_get_default()
    return GA_Pgroup_nnodes(pgroup)

def pgroup_nodeid(int pgroup=-1):
    """Returns the relative index of the processor in the processor group
    specified by pgroup.
    
    This index will generally differ from the absolute processor index
    returned by ga.nodeid if the processor group is not the world group.

    This is a local operation. 

    :Parameters:
        pgroup : int
            the group handle

    """
    if pgroup < 0:
        pgroup = pgroup_get_default()
    return GA_Pgroup_nodeid(pgroup)

def pgroup_set_default(int pgroup=-1):
    """Resets the default processor group on a collection of processors.
    
    All processors in the group referenced by p_handle must make a call to
    this function. Any standard global array call that is made after resetting
    the default processor group will be restricted to processors in that
    group. Global arrays that are created after resetting the default
    processor group will only be defined on that group and global operations
    such as ga.sync or ga.gop will be restricted to processors in that group.
    The ga.pgroup_set_default call can be used to rapidly convert large
    applications, written with GA, into routines that run on processor groups.

    The default processor group can be overridden by using GA calls that
    require an explicit group handle as one of the arguments.

    This is a collective operation on the group represented by the handle
    pgroup. 

    """
    if pgroup < 0:
        pgroup = pgroup_get_world()
    GA_Pgroup_set_default(pgroup)

def pgroup_split(int pgroup, int num_group):
    """TODO"""
    return GA_Pgroup_split(pgroup, num_group)

def pgroup_split_irreg(int pgroup, int color):
    """TODO"""
    return GA_Pgroup_split_irreg(pgroup, color)

def pgroup_sync(int pgroup=-1):
    """Executes a synchronization group across the processors in the processor
    group specified by pgroup.
    
    Nodes outside this group are unaffected.

    This is a collective operation on the processor group specified by
    pgroup.  

    """
    if pgroup < 0:
        pgroup = pgroup_get_default()
    GA_Pgroup_sync(pgroup)

def print_distribution(int g_a):
    """Prints the array distribution.

    This is a collective operation. 

    """
    GA_Print_distribution(g_a)

def print_file(int g_a, file):
    """Prints an entire array to a file.

    This is a collective operation. 

    :Parameters:
        file : file-like
            file-like object which must implement fileno(), or a string
        g_a : int
            the array handle

    """
    #GA_Print_file(file.fileno(), g_a)
    raise NotImplementedError

def print_patch(int g_a, lo=None, hi=None, bint pretty=True):
    """Prints a patch of g_a array to the standard output.
    
    If pretty is False then output is printed in a dense fashion. If
    pretty is True then output is formatted and rows/columns labeled.

    This is a collective operation. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd
    cdef int apretty=0
    lo_nd,hi_nd = _lohi(g_a,lo,hi)
    if pretty:
        apretty = 1
    NGA_Print_patch64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data, apretty)

def print_stats():
    """This non-collective (MIMD) operation prints information about:

    * number of calls to the GA create/duplicate, destroy, get, put,
    * scatter, gather, and read_and_inc operations
    * total amount of data moved in the GA primitive operations
    * amount of data moved in GA primitive operations to logicaly
    * remote locations
    * maximum memory consumption in global arrays, and
    * number of requests serviced in the interrupt-driven
    * implementations by the calling process.

    This operation is local. 

    """
    GA_Print_stats()

def print_stdout(int g_a):
    """Prints an entire array to the standard output."""
    GA_Print(g_a)

def proc_topology(int g_a, int proc):
    """Based on the distribution of an array associated with handle g_a,
    determines coordinates of the specified processor in the virtual processor
    grid corresponding to the distribution of array g_a.
    
    The numbering starts from 0. The values of -1 means that the processor
    doesn't 'own' any section of array represented by g_a.

    This operation is local. 

    """
    cdef int ndim = GA_Ndim(g_a)
    cdef np.ndarray[np.int32_t, ndim=1] coord
    coord = np.zeros(ndim, dtype=np.intc)
    NGA_Proc_topology(g_a, proc, <int*>coord.data)
    return coord

def put(int g_a, buffer, lo=None, hi=None):
    """Copies data from local array buffer to the global array section.
    
    The local array is assumed to be have the same number of dimensions as the
    global array.  Any detected inconsitencies/errors in input arguments are
    fatal.

    This is a one-sided operation. 

    :Parameters:
        g_a : int
            the array handle
        buffer : array-like
            the data to put
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : array-like of integers
            higher bound patch coordinates, exclusive

    """
    _put_common(g_a, buffer, lo, hi)

cdef _put_common(int g_a, buffer, lo=None, hi=None,
        bint nb=False, bint periodic=False, skip=None):
    """Copies data from local array buffer to the global array section.
    
    The local array is assumed to have the same shape as the requested region,
    or the local array can be 1-dimensional so long as it has the same number
    of elements as the requested region. Any detected inconsitencies raise a
    ValueError.

    This is a one-sided operation. 

    :Parameters:
        g_a : int
            the array handle
        buffer : array-like
            the data to put;
            should either be 1D and len(buffer)==np.prod(hi-lo), or
            np.all(buffer.shape == hi-lo) i.e. buffer is 1D and same size as
            requested region or buffer is the same shape as requested region
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : array-like of integers
            higher bound patch coordinates, exclusive
        nb : bool
            whether this call is non-blocking (see ga.nbget)
        periodic : bool
            whether this call is periodic (see ga.periodic_get)
        skip : 1D array-like of integers
            strides for each dimension

    :returns: None, usually.  However if nb=True, the nonblocking handle is returned.

    """
    cdef np.ndarray buffer_nd
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd, ld_nd, shape, skip_nd
    cdef int gtype=inquire_type(g_a)
    cdef ga_nbhdl_t nbhandle
    dtype = _to_dtype[gtype]
    lo_nd,hi_nd = _lohi(g_a,lo,hi)
    shape = hi_nd-lo_nd+1
    if skip is None:
        skip_nd = None
    else:
        skip_nd = _inta64(skip)
        shape = (hi_nd-lo_nd)/skip_nd+1
    buffer_nd = np.asarray(buffer, dtype=dtype)
    if buffer_nd.dtype != dtype:
        raise ValueError, "buffer is wrong type :: buffer=%s != %s" % (
                buffer.dtype, dtype)
    # Due to GA restrictions, buffer must not have negative strides
    # and buffer's last stride must be same as itemsize
    strides = [buffer_nd.strides[i]/buffer_nd.itemsize
            for i in range(buffer_nd.ndim)]
    if (strides and (strides[-1] != 1 or np.any(np.asarray(strides) < 0))):
        buffer_nd = np.ascontiguousarray(buffer_nd)
    # we allow 1-d "flat" buffers in addition to buffers matching the shape of
    # the requested region
    if buffer_nd.ndim == 1:
        if buffer_nd.size != np.prod(shape):
            raise ValueError, ('buffer size does not match shape :: '
                    'buffer.size=%s != np.prod(shape)=%s' % (
                    buffer_nd.size, np.prod(shape)))
        ld_nd = shape[1:]
    else:
        buffer_shape = [buffer_nd.shape[i] for i in range(buffer_nd.ndim)]
        if not np.all(buffer_shape == shape):
            raise ValueError, ('buffer shape does not match request shape :: '
                    'buffer_shape=%s != shape=%s' % (
                    buffer_shape, shape))
        ld_nd = np.asarray([strides[i]/strides[i+1]
                for i in range(buffer_nd.ndim-1)], dtype=np.int64)
    if nb:
        NGA_NbPut64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <void*>buffer_nd.data, <int64_t*>ld_nd.data, &nbhandle)
        return nbhandle
    elif periodic:
        NGA_Periodic_put64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <void*>buffer_nd.data, <int64_t*>ld_nd.data)
    elif skip is not None:
        NGA_Strided_put64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <int64_t*>skip_nd.data,
                <void*>buffer_nd.data, <int64_t*>ld_nd.data)
    else:
        NGA_Put64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data,
                <void*>buffer_nd.data, <int64_t*>ld_nd.data)

def randomize(int g_a, val=None):
    """Fill array with random values in [0,val)."""
    cdef int gtype=inquire_type(g_a)
    cdef int            ival
    cdef long           lval
    cdef long long      llval
    cdef float          fval
    cdef double         dval
    cdef long double    ldval
    cdef SingleComplex  fcval
    cdef DoubleComplex  dcval
    cdef void          *vval=NULL
    if val is None:
        val = 1
    vval = _convert_multiplier(gtype, val,
            &ival,  &lval,  &llval,
            &fval,  &dval,  &ldval,
            &fcval, &dcval)
    GA_Randomize(g_a, vval)

def read_inc(int g_a, subscript, long inc=1):
    """Atomically read and increment an element in an integer array. 

    This is a one-sided and atomic operation.

    :Parameters:
        g_a : int
            the array handle
        subscript : 1D array-like of integers
            index for the referenced element
        inc : long
            the increment

    """
    cdef np.ndarray[np.int64_t, ndim=1] subscript_nd
    subscript_nd = _inta64(subscript)
    return NGA_Read_inc64(g_a, <int64_t*>subscript_nd.data, inc)

def register_dtype(dtype):
    """Creates a new data type based on the given dtype.

    :Parameters:
        dtype : dtype
            the numpy dtype to register

    """
    cdef int gatype
    dtype = np.dtype(dtype) # just in case it's not really a dtype instance
    gatype = NGA_Register_type(dtype.itemsize)
    _to_dtype[gatype] = dtype
    return gatype
    
def register_type(size_t bytes):
    """Creates a new data type of size bytes.

    :Parameters:
        bytes : size_t
            the size of the new data type

    """
    return NGA_Register_type(bytes)
    
def recip(int g_a, lo=None, hi=None):
    """Take element-wise reciprocal of the array or patch.

    This is a collective operation. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd
    if lo is None and hi is None:
        GA_Recip(g_a)
    else:
        lo_nd,hi_nd = _lohi(g_a,lo,hi)
        GA_Recip_patch64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data)

def release(int g_a, lo=None, hi=None):
    """Releases access to a global array when the data was read only.

    Your code should look like:

        array = ga.access(g_a)
        # <operate on the data referenced by ptr> 
        ga.release(g_a)

    NOTE: see restrictions specified for ga.access

    This operation is local. 
    
    """
    _release_common(g_a, lo, hi, False)

def release_block(int g_a, int index):
    """Releases access to the block of data specified by the integer index
    when data was accessed as read only.
    
    This is only applicable to block-cyclic data distributions created using
    the simple block-cyclic distribution. This is a local operation.

    """
    NGA_Release_block(g_a, index)

def release_block_grid(int g_a, subscript):
    """Releases access to the block of data specified by the subscript array
    when data was accessed as read only.
    
    This is only applicable to block-cyclic data distributions created using
    the SCALAPACK data distribution.
    
    This is a local operation.

    """
    cdef np.ndarray[np.int32_t, ndim=1] subscript_nd
    subscript_nd = _inta32(subscript)
    NGA_Release_block_grid(g_a, <int*>subscript_nd.data)

def release_block_segment(int g_a, int proc):
    """Releases access to the block of locally held data for a block-cyclic
    array, when data was accessed as read-only.
    
    This is a local operation.

    """
    NGA_Release_block_segment(g_a, proc)

cdef _release_common(int g_a, lo, hi, bint update):
    """TODO"""
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd, lo_dst, hi_dst
    # first things first, if no data is owned, return silently
    lo_dst,hi_dst = distribution(g_a)
    # convet hi_dst back to GA inclusive indexing convention
    hi_dst -= 1
    if lo_dst[0] < 0 or hi_dst[0] < 0:
        return
    if lo is not None:
        lo_nd = _inta64(lo)
    else:
        lo_nd = lo_dst
    if hi is not None:
        hi_nd = _inta64(hi)
    else:
        hi_nd = hi_dst
    # sanity checks
    if np.sometrue(lo_nd>hi_nd):
        raise ValueError,"lo>hi lo=%s hi=%s"%(lo_nd,hi_nd)
    if np.sometrue(lo_nd<lo_dst):
        raise ValueError,"lo out of bounds lo_dst=%s lo=%s"%(lo_dst,lo_nd)
    if np.sometrue(hi_nd>hi_dst):
        raise ValueError,"hi out of bounds hi_dst=%s hi=%s"%(hi_dst,hi_nd)
    if update:
        NGA_Release_update64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data)
    else:
        NGA_Release64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data)

def release_ghost_element(int g_a, subscript):
    """Releases access to the locally held data for an array with ghost
    elements, when data was accessed as read-only.
    
    This is a local operation.

    """
    cdef np.ndarray[np.int64_t, ndim=1] subscript_nd
    subscript_nd = _inta64(subscript)
    NGA_Release_ghost_element64(g_a, <int64_t*>subscript_nd.data)

def release_ghosts(int g_a):
    """Releases access to the locally held block of data containing ghost
    elements, when data was accessed as read-only.
    
    This is a local operation.

    """
    NGA_Release_ghosts(g_a)

def release_update(int g_a, lo=None, hi=None):
    """Releases access to the data.
    
    It must be used if the data was accessed for writing.
    NOTE: see restrictions specified for ga.access.
    
    This operation is local. 
    
    """
    _release_common(g_a, lo, hi, True)

def release_update_block(int g_a, int index):
    """Releases access to the block of data specified by the integer index
    when data was accessed in read-write mode.
    
    This is only applicable to block-cyclic data distributions created using
    the simple block-cyclic distribution.
    
    This is a local operation.

    """
    NGA_Release_update_block(g_a, index)

def release_update_block_grid(int g_a, subscript):
    """Releases access to the block of data specified by the subscript array
    when data was accessed in read-write mode.
    
    This is only applicable to block-cyclic data distributions created using
    the SCALAPACK data distribution.
    
    This is a local operation.

    """
    cdef np.ndarray[np.int32_t, ndim=1] subscript_nd
    subscript_nd = _inta32(subscript)
    NGA_Release_update_block_grid(g_a, <int*>subscript_nd.data)

def release_update_block_segment(int g_a, int proc):
    """Releases access to the block of locally held data for a block-cyclic
    array, when data was accessed as read-only.
    
    This is a local operation.

    """
    NGA_Release_update_block_segment(g_a, proc)

def release_update_ghost_element(int g_a, subscript):
    """Releases access to the locally held data for an array with ghost
    elements, when data was accessed in read-write mode.
    
    This is a local operation.

    """
    cdef np.ndarray[np.int64_t, ndim=1] subscript_nd
    subscript_nd = _inta64(subscript)
    NGA_Release_update_ghost_element64(g_a, <int64_t*>subscript_nd.data)
    pass

def release_update_ghosts(int g_a):
    """Releases access to the locally held block of data containing ghost
    elements, when data was accessed in read-write mode.
    
    This is a local operation. 

    """
    NGA_Release_update_ghosts(g_a)

def scale(int g_a, value, lo=None, hi=None):
    """Scales an array by the constant s.
    
    Note that the library is unable to detect errors when the pointed value is
    of different type than the array.

    This is a collective operation. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd
    cdef int gtype=inquire_type(g_a)
    cdef int            ivalue
    cdef long           lvalue
    cdef long long      llvalue
    cdef float          fvalue
    cdef double         dvalue
    cdef long double    ldvalue
    cdef SingleComplex  fcvalue
    cdef DoubleComplex  dcvalue
    cdef void          *vvalue
    vvalue = _convert_multiplier(gtype, value,
            &ivalue,  &lvalue,  &llvalue,
            &fvalue,  &dvalue,  &ldvalue,
            &fcvalue, &dcvalue)
    if lo is None and hi is None:
        GA_Scale(g_a, vvalue)
    else:
        lo_nd,hi_nd = _lohi(g_a,lo,hi)
        NGA_Scale_patch64(
                g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data, vvalue)

def scale_rows(int g_a, int g_v):
    """Scales the rows of this matrix g_a using the vector g_v.

    This is a collective operation. 

    """
    GA_Scale_rows(g_a, g_v)
            
def scale_cols(int g_a, int g_v):
    """Scales the columns of this matrix g_a using the vector g_v.

    This is a collective operation. 

    """
    GA_Scale_cols(g_a, g_v)

def scan_add(int g_src, int g_dst, int g_msk, lo=None, hi=None,
        bint excl=False):
    """Adds successive elements in a source vector g_src and put the results
    in a destination vector g_dst.
    
    The addition will restart based on the values of the integer mask vector
    g_msk. The scan is performed within the range specified by the integer
    values lo and hi. Note that this operation can only be applied to
    1-dimensional arrays. The excl flag determines whether the sum starts with
    the value in the source vector corresponding to the location of a 1 in the
    mask vector (excl=False) or whether the first value is set equal to 0
    (excl=True). Some examples of this operation are given below.

    ga.scan_add(g_src, g_dst, g_msk, 0, n, False);
    g_msk:   1  0  0  0  0  0  1  0  1  0  0  1  0  0  1  1  0
    g_src:   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
    g_dst:   1  3  6 10 16 21  7 15  9 19 30 12 25 39 15 16 33

    ga.scan_add(g_src, g_dst, g_msk, 0, n, True);
    g_msk:   1  0  0  0  0  0  1  0  1  0  0  1  0  0  1  1  0
    g_src:   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
    g_dst:   0  1  3  6 10 15  0  7  0  9 19  0 12 25  0  0 16

    This is a collective operation.

    :Parameters:
        g_src : int
            handle for source arrray
        g_dst : int
            handle for destination array
        g_msk : int
            handle for integer array representing mask
        lo : 1D array-like of integers
            low value of range on which operation is performed
        hi : 1D array-like of integers
            hi value of range on which operation is performed
        excl : bool
            whether the first value is set to 0 (see above)

    """
    cdef np.ndarray[np.int64_t, ndim=1] hi_nd = inquire_dims(g_src)-1
    cdef int64_t c_lo=0, c_hi=hi_nd[0]
    cdef int c_excl=0
    if lo is not None:
        c_lo = lo
    if hi is not None:
        c_hi = hi
    if excl:
        c_excl = 1
    GA_Scan_add64(g_src, g_dst, g_msk, c_lo, c_hi, c_excl)

def scan_copy(int g_src, int g_dst, int g_msk, lo=None, hi=None):
    """This subroutine does a segmented scan-copy of values in the source
    array g_src into a destination array g_dst with segments defined by
    values in the integer mask array g_msk. The scan-copy operation is only
    applied to the range between the lo and hi indices. This operation is
    restriced to 1-dimensional arrays. The resulting destination array will
    consist of segments of consecutive elements with the same value. An
    example is shown below

    GA_Scan_copy(g_src, g_dst, g_msk, 0, n);
    g_msk:   1  0  0  0  0  0  1  0  1  0  0  1  0  0  1  1  0
    g_src:   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
    g_dst:   1  1  1  1  1  1  7  7  9  9  9 12 12 12 15 16 16

    This is a collective operation.

    """
    cdef np.ndarray[np.int64_t, ndim=1] hi_nd = inquire_dims(g_src)-1
    cdef int64_t c_lo=0, c_hi=hi_nd[0]
    cdef int c_excl=0
    if lo is not None:
        c_lo = lo
    if hi is not None:
        c_hi = hi
    GA_Scan_copy64(g_src, g_dst, g_msk, c_lo, c_hi)

def scatter(int g_a, values, subsarray):
    """Scatters array elements from a global array into a local array.

    subsarray will be converted to an ndarray if it is not one already.  A
    two-dimensional array is allowed so long as its shape is (n,ndim) where n
    is the number of elements to gather and ndim is the number of dimensions
    of the target array.  Also, subsarray must be contiguous.

    For example, if the subsarray were two-dimensional::

        for k in range(n):
            v[k] = g_a[subsarray[k,0],subsarray[k,1],subsarray[k,2]...]

    For example, if the subsarray were one-dimensional::

        for k in range(n):
            base = n*ndim
            v[k] = g_a[subsarray[base+0],subsarray[base+1],subsarray[base+2]...]

    This is a one-sided operation. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] subsarray1_nd = None
    cdef np.ndarray[np.int64_t, ndim=2] subsarray2_nd = None
    cdef np.ndarray values_nd = None
    cdef int gtype = inquire_type(g_a)
    cdef int ndim = GA_Ndim(g_a)
    cdef int64_t n
    # prepare subsarray
    try:
        subsarray1_nd = np.ascontiguousarray(subsarray, dtype=np.int64)
        n = len(subsarray1_nd) / ndim
    except ValueError:
        subsarray1_nd = None
        try:
            subsarray2_nd = np.ascontiguousarray(subsarray, dtype=np.int64)
            n = len(subsarray2_nd) # length of first dimension of subsarray2_nd
        except ValueError:
            raise ValueError, "subsarray must be either 1- or 2-dimensional"
    # prepare values array
    values_nd = np.asarray(values, dtype=_to_dtype[gtype])
    if values_nd.ndim != 1:
        raise ValueError, "values must be one-dimensional"
    if not values_nd.flags['C_CONTIGUOUS']:
        raise ValueError, "values must be contiguous"
    if len(values_nd) < n:
        raise ValueError, "values was not large enough"
    # call the wrapped function
    if subsarray1_nd is not None:
        NGA_Scatter_flat64(g_a, <void*>values_nd.data,
                <int64_t*>subsarray1_nd.data, n)
    elif subsarray2_nd is not None:
        NGA_Scatter_flat64(g_a, <void*>values_nd.data,
                <int64_t*>subsarray2_nd.data, n)
    else:
        raise ValueError, "how did this happen?"

def scatter_acc(int g_a, values, subsarray, alpha=None):
    """Scatters array elements from a global array into a local array.

    Like scatter, but adds values to existing values in the global array after
    multiplying by alpha.

    subsarray will be converted to an ndarray if it is not one already.  A
    two-dimensional array is allowed so long as its shape is (n,ndim) where n
    is the number of elements to gather and ndim is the number of dimensions
    of the target array.  Also, subsarray must be contiguous.

    For example, if the subsarray were two-dimensional::

        for k in range(n):
            v[k] = g_a[subsarray[k,0],subsarray[k,1],subsarray[k,2]...]

    For example, if the subsarray were one-dimensional::

        for k in range(n):
            base = n*ndim
            v[k] = g_a[subsarray[base+0],subsarray[base+1],subsarray[base+2]...]

    This is a one-sided operation. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] subsarray1_nd = None
    cdef np.ndarray[np.int64_t, ndim=2] subsarray2_nd = None
    cdef np.ndarray values_nd = None
    cdef int gtype = inquire_type(g_a)
    cdef int ndim = GA_Ndim(g_a)
    cdef int64_t n
    cdef int            ialpha
    cdef long           lalpha
    cdef long long      llalpha
    cdef float          falpha
    cdef double         dalpha
    cdef long double    ldalpha
    cdef SingleComplex  fcalpha
    cdef DoubleComplex  dcalpha
    cdef void          *valpha=NULL
    # prepare subsarray
    try:
        subsarray1_nd = np.ascontiguousarray(subsarray, dtype=np.int64)
        n = len(subsarray1_nd) / ndim
    except ValueError:
        subsarray1_nd = None
        try:
            subsarray2_nd = np.ascontiguousarray(subsarray, dtype=np.int64)
            n = len(subsarray2_nd) # length of first dimension of subsarray2_nd
        except ValueError:
            raise ValueError, "subsarray must be either 1- or 2-dimensional"
    # prepare values array
    values_nd = np.asarray(values, dtype=_to_dtype[gtype])
    if values_nd.ndim != 1:
        raise ValueError, "values must be one-dimensional"
    if not values_nd.flags['C_CONTIGUOUS']:
        raise ValueError, "values must be contiguous"
    if len(values_nd) < n:
        raise ValueError, "values was not large enough"
    # prepare alpha
    if alpha is None:
        alpha = 1
    valpha = _convert_multiplier(gtype, alpha,
            &ialpha,  &lalpha,  &llalpha,
            &falpha,  &dalpha,  &ldalpha,
            &fcalpha, &dcalpha)
    # call the wrapped function
    if subsarray1_nd is not None:
        NGA_Scatter_acc_flat64(g_a, <void*>values_nd.data,
                <int64_t*>subsarray1_nd.data, n, valpha)
    elif subsarray2_nd is not None:
        NGA_Scatter_acc_flat64(g_a, <void*>values_nd.data,
                <int64_t*>subsarray2_nd.data, n, valpha)
    else:
        raise ValueError, "how did this happen?"

def select_elem(int g_a, char *op):
    """Returns the value and index for an element that is selected by the
    specified operator in a global array corresponding to g_a handle.

    This is a collective operation. 

    :returns: the selected element and the array index for the selected element

    """
    cdef np.ndarray[np.int64_t, ndim=1] index
    cdef int gtype=inquire_type(g_a)
    cdef int            ialpha
    cdef long           lalpha
    cdef long long      llalpha
    cdef float          falpha
    cdef double         dalpha
    cdef long double    ldalpha
    cdef SingleComplex  fcalpha
    cdef DoubleComplex  dcalpha
    cdef void          *valpha=NULL
    valpha = _convert_multiplier(gtype, 0,
            &ialpha,  &lalpha,  &llalpha,
            &falpha,  &dalpha,  &ldalpha,
            &fcalpha, &dcalpha)
    index = np.ndarray(GA_Ndim(g_a), dtype=np.int64)
    NGA_Select_elem64(g_a, op, valpha, <int64_t*>index.data)
    if gtype == C_INT:
        return ialpha,index
    elif gtype == C_LONG:
        return lalpha,index
    elif gtype == C_LONGLONG:
        return llalpha,index
    elif gtype == C_FLOAT:
        return falpha,index
    elif gtype == C_DBL:
        return dalpha,index
    elif gtype == C_LDBL:
        return ldalpha,index
    elif gtype == C_SCPL:
        # TODO explicitly convert GA complex to Python complex type
        return fcalpha,index
    elif gtype == C_DCPL:
        # TODO explicitly convert GA complex to Python complex type
        return dcalpha,index
    else:
        raise TypeError, "type of g_a not recognized"

def select_elem_min(int g_a):
    """Equivalent to ga.select_elem(g_a, "min")."""
    return select_elem(g_a, "min")

def select_elem_max(int g_a):
    """Equivalent to ga.select_elem(g_a, "max")."""
    return select_elem(g_a, "max")

def set_array_name(int g_a, char *name):
    """Assigns a unique character string name to a global array handle that
    was obtained using the GA_Create_handle function.

    This is a collective operation.

    """
    GA_Set_array_name(g_a, name)

def set_block_cyclic(int g_a, dims):
    """Creates a global array with a simple block-cyclic data distribution.

    The array is broken up into blocks of size dims and each block is numbered
    sequentially using a column major indexing scheme. The blocks are then
    assigned in a simple round-robin fashion to processors. This is
    illustrated in the figure below for an array containing 25 blocks
    distributed on 4 processors. Blocks at the edge of the array may be
    smaller than the block size specified in dims. In the example below,
    blocks 4,9,14,19,20,21,22,23, and 24 might be smaller thatn the remaining
    blocks. Most global array operations are insensitive to whether or not a
    block-cyclic data distribution is used, although performance may be slower
    in some cases if the global array is using a block-cyclic data
    distribution. Individual data blocks can be accessesed using the
    block-cyclic access functions.

    This is a collective operation.

    """
    cdef np.ndarray[np.int32_t, ndim=1] dims_nd
    dims_nd = _inta32(dims)
    GA_Set_block_cyclic(g_a, <int*>dims_nd.data)

def set_block_cyclic_proc_grid(int g_a, block, proc_grid):
    """Creates a global array with a SCALAPACK-type block cyclic data
    distribution.
    
    The user specifies the dimensions of the processor grid in the array
    proc_grid. The product of the processor grid dimensions must equal the
    number of total number of processors  and the number of dimensions in the
    processor grid must be the same as the number of dimensions in the global
    array. The data blocks are mapped onto the processor grid in a cyclic
    manner along each of the processor grid axes.  This is illustrated below
    for an array consisting of 25 data blocks disributed on 6 processors. The
    6 processors are configured in a 3 by 2 processor grid. Blocks at the edge
    of the array may be smaller than the block size specified in dims. Most
    global array operations are insensitive to whether or not a block-cyclic
    data distribution is used, although performance may be slower in some
    cases if the global array is using a block-cyclic data distribution.
    Individual data blocks can be accessesed using the block-cyclic access
    functions.

    This is a collective operation.

    """
    cdef np.ndarray[np.int32_t, ndim=1] block_nd, proc_grid_nd
    block_nd = _inta32(block)
    proc_grid_nd = _inta32(proc_grid)
    GA_Set_block_cyclic_proc_grid(g_a,
            <int*>block_nd.data,
            <int*>proc_grid_nd.data)

def set_chunk(int g_a, chunk):
    """This function is used to set the chunk array for a global array handle
    that was obtained using the GA_Create_handle function. The chunk array is
    used to determine the minimum number of array elements assigned to each
    processor along each coordinate direction.

    This is a collective operation.
    
    """
    cdef np.ndarray[np.int64_t, ndim=1] chunk_nd
    chunk_nd = _inta64(chunk)
    GA_Set_chunk64(g_a, <int64_t*>chunk_nd.data)

def set_data(int g_a, dims, int type):
    """Sets the array dimension, the coordinate dimensions, and the data type
    assigned to a global array handle obtained using the ga.create_handle
    function.

    This is a collective operation.

    """
    cdef np.ndarray[np.int64_t, ndim=1] dims_nd
    dims_nd = _inta64(dims)
    GA_Set_data64(g_a, len(dims_nd), <int64_t*>dims_nd.data, type)

def set_debug(bint debug):
    """Sets an internal flag in the GA library to either True or False.
    
    The value of this flag can be recovered at any time using the ga.get_debug
    function. The flag is set to false when the the GA library is initialized.
    This can be useful in a number of debugging situations, especially when
    examining the behavior of routines that are called in multiple locations
    in a code.

    This is a local operation.

    """
    GA_Set_debug(debug)

def set_diagonal(int g_a, int g_v):
    """Sets the diagonal elements of this matrix g_a with the elements of the
    vector g_v.

    This is a collective operation.

    """
    GA_Set_diagonal(g_a, g_v)

def set_ghosts(int g_a, width):
    """Sets the ghost cell widths for a global array handle that was obtained
    using the ga.create_handle function.
    
    The ghosts cells widths indicate how many ghost cells are used to pad the
    locally held array data along each dimension. The padding can be set
    independently for each coordinate dimension.

    This is a collective operation.

    """
    cdef np.ndarray[np.int64_t, ndim=1] width_nd
    width_nd = _inta64(width)
    GA_Set_ghosts64(g_a, <int64_t*>width_nd.data)

def set_irreg_distr(int g_a, mapc, nblock):
    """Partitions the array data among the individual processors for a global
    array handle obtained using the ga.create_handle function.

    The distribution is specified as a Cartesian product of distributions for
    each dimension. For example, the following figure demonstrates
    distribution of a 2-dimensional array 8x10 on 6 (or more) processors.
    nblock(2)={3,2}, the size of mapc array is s=5 and array mapc contains the
    following elements mapc={1,3,7, 1, 6}. The distribution is nonuniform
    because, P1 and P4 get 20 elements each and processors P0,P2,P3, and P5
    only 10 elements each.
     
    +----+----++--+
    |  5 |  5 ||  |
    +====+====++==+
    | P0 | P3 || 2|
    | P1 | P4 || 4|
    | P2 | P5 || 2|
    +----+----++--+

    The array width() is used to control the width of the ghost cell boundary
    around the visible data on each processor. The local data of the global
    array residing on each processor will have a layer width(n) ghosts cells
    wide on either side of the visible data along the dimension n.

    This is a collective operation.

    """
    cdef np.ndarray[np.int64_t, ndim=1] mapc_nd, nblock_nd
    mapc_nd = _inta64(mapc)
    nblock_nd = _inta64(nblock)
    GA_Set_irreg_distr64(g_a, <int64_t*>mapc_nd.data, <int64_t*>nblock_nd.data)

def set_memory_limit(size_t limit):
    """Sets the amount of memory to be used (in bytes) per process.

    This is a local operation. 

    :Parameters:
        limit : size_t
            the amount of memory in bytes per process

    """
    GA_Set_memory_limit(limit)

def set_pgroup(int g_a, int pgroup):
    """Sets the processor configuration assigned to a global array handle that
    was obtained using the ga.create_handle function.
    
    It can be used to create mirrored arrays by using the mirrored array
    processor configuration in this function call. It can also be used to
    create an array on a processor group by using a processor group handle in
    this call.

    This is a collective operation.

    """
    GA_Set_pgroup(g_a, pgroup)

def set_restricted(int g_a, list):
    """Restrict data in the global array g_a to only the processors listed in
    the array list.
    
    len(list) must be less than or equal to the number of available processors.
    If this call is used in conjunction with set_irreg_distr, then the
    decomposition in the set_irreg_distr call must be done assuming that the
    number of processors is nproc. The data that ordinarily would be mapped to
    process 0 is mapped to the process in list[0], the data that would be
    mapped to process 1 will be mapped to list[1], etc. This can be used to
    remap the data distribution to different processors, even if nproc equals
    the number of available processors.

    This is a collective operation.

    """
    cdef np.ndarray[np.int32_t, ndim=1] list_nd
    list_nd = _inta32(list)
    GA_Set_restricted(g_a, <int*>list_nd.data, len(list_nd))

def set_restricted_range(int g_a, int lo_proc, int hi_proc):
    """Restrict data in the global array to the given range of processors.

    Both lo_proc and hi_proc must be less than or equal to the total number of
    processors minus one (e.g., in the range [0,N-1], where N is the total
    number of processors) and lo_proc must be less than or equal to hi_proc. If
    lo_proc = 0 and hi_proc = N-1 then this call has no effect on the data
    distribution.

    This is a collective operation.

    """
    GA_Set_restricted_range(g_a, lo_proc, hi_proc)

def shift_diagoal(int g_a, value=None):
    """Adds this constant to the diagonal elements of the matrix.

    This is a collective operation.

    """
    cdef int            ivalue
    cdef long           lvalue
    cdef long long      llvalue
    cdef float          fvalue
    cdef double         dvalue
    cdef long double    ldvalue
    cdef SingleComplex  fcvalue
    cdef DoubleComplex  dcvalue
    cdef void          *vvalue
    cdef int gtype=inquire_type(g_a)
    if value is None:
        value = 1
    vvalue = _convert_multiplier(gtype, value,
            &ivalue,  &lvalue,  &llvalue,
            &fvalue,  &dvalue,  &ldvalue,
            &fcvalue, &dcvalue)
    GA_Shift_diagonal(g_a, vvalue)

def solve(int g_a, int g_b):
    """Solves a system of linear equations A * X = B.

    It first will call the Cholesky factorization routine and, if sucessfully,
    will solve the system with the Cholesky solver. If Cholesky will be not be
    able to factorize A, then it will call the LU factorization routine and
    will solve the system with forward/backward substitution. On exit B will
    contain the solution X.

    This is a collective operation.

    :returns: 0 if Cholesky factoriztion was succesful.  >0 if the leading minor of this order is not positive definite, Cholesky factorization could not be completed and LU factoriztion was used

    """
    return GA_Solve(g_a, g_b)

def spd_invert(int g_a):
    """Compute the inverse of a double precision using the Cholesky
    factorization of a NxN double precision symmetric positive definite matrix
    A stored in the global array represented by g_a. On successful exit, A
    will contain the inverse.

    This is a collective operation.

    :returns: 0 if successful exit; >0 if the leading minor of this order is not positive definite and the factorization could not be completed; <0 if it returns the index i of the (i,i) element of the factor L/U that is zero and the inverse could not be computed

    """
    return GA_Spd_invert(g_a)

def step_max(int g_a, int g_b, alo=None, ahi=None, blo=None, bhi=None):
    """Calculates the largest multiple of a vector g_b that can be added to
    this vector g_a while keeping each element of this vector non-negative.

    This is a collective operation. 

    """
    cdef np.ndarray[np.int64_t, ndim=1] alo_nd, ahi_nd
    cdef np.ndarray[np.int64_t, ndim=1] blo_nd, bhi_nd
    cdef double step
    if (alo is None and ahi is None
            and blo is None and bhi is None):
        GA_Step_max(g_a, g_b, &step)
    else:
        alo_nd,ahi_nd = _lohi(g_a,alo,ahi)
        blo_nd,bhi_nd = _lohi(g_b,blo,bhi)
        GA_Step_max_patch64(g_a, <int64_t*>alo_nd.data, <int64_t*>ahi_nd.data,
                g_b, <int64_t*>blo_nd.data, <int64_t*>bhi_nd.data, &step)
    return step

def strided_acc(int g_a, buffer, lo=None, hi=None, skip=None, alpha=None):
    """Strided version of ga.acc.
    
    The values corresponding to dimension n in buf are accumulated to every
    skip[n] values of the global array g_a.
    
    Combines data from buffer with data in the global array patch.
    
    The buffer array is assumed to be have the same number of dimensions as
    the global array.  If the buffer is not contiguous, a contiguous copy will
    be made.
    
        global array section (lo[],hi[]) += alpha * buffer

    This is a one-sided and atomic operation.

    :Parameters:
        g_a : int
            the array handle
        buffer : array-like
            must be contiguous and have same number of elements as patch
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : 1D array-like of integers
            higher bound patch coordinates, exclusive
        alpha : object
            multiplier (converted to the appropriate type)

    """
    _acc_common(g_a, buffer, lo, hi, alpha, False, False, skip)
       
def strided_get(int g_a, lo=None, hi=None, skip=None, np.ndarray buffer=None):
    """Strided version of ga.get.
    
    Copies data from global array section to the local array buffer.
    
    The local array is assumed to be have the same number of dimensions as the
    global array. Any detected inconsitencies/errors in the input arguments
    are fatal.

    This is a one-sided operation.

    :Parameters:
        g_a : int
            the array handle
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : 1D array-like of integers
            higher bound patch coordinates, exclusive
        skip : 1D array-like of integers
            strides for each dimension
        buffer : ndarray
            an ndarray of the appropriate type, large enough to hold lo,hi

    :returns: The local array buffer.
    
    """
    return _get_common(g_a, lo, hi, buffer, False, False, skip)

def strided_put(int g_a, buffer, lo=None, hi=None, skip=None):
    """Strided version of ga.put.
    
    Copies data from local array buffer to the global array section.
    
    The local array is assumed to be have the same number of dimensions as the
    global array.  Any detected inconsitencies/errors in input arguments are
    fatal.

    This is a one-sided operation. 

    :Parameters:
        g_a : int
            the array handle
        buffer : array-like
            the data to put
        lo : 1D array-like of integers
            lower bound patch coordinates, inclusive
        hi : array-like of integers
            higher bound patch coordinates, exclusive
        skip : 1D array-like of integers
            strides for each dimension

    """
    _put_common(g_a, buffer, lo, hi, False, False, skip)

def summarize(bint verbose):
    """Prints info about allocated arrays."""
    GA_Summarize(verbose)

def symmetrize(int g_a):
    """Symmetrizes matrix A represented with handle g_a: A:= .5 * (A+A').

    This is a collective operation.

    """
    GA_Symmetrize(g_a)

def sync():
    """Synchronize processes (a barrier) and ensure that all GA operations
    completed.

    This is a collective operation.

    """
    GA_Sync()

def terminate():
    """Delete all active arrays and destroy internal data structures.

    This is a collective operation. 

    """
    global _initialized
    _initialized = False
    GA_Terminate()

def total_blocks(int g_a):
    """Returns the total number of blocks contained in a global
    array with a block-cyclic data distribution.
    
    This is a local operation.

    """
    return GA_Total_blocks(g_a)

def transpose(int g_a, int g_b):
    """Transposes a matrix: B = A', where A and B are represented by handles
    g_a and g_b.

    This is a collective operation.

    """
    GA_Transpose(g_a, g_b)

def unlock(int mutex):
    """Unlocks a mutex object identified by the mutex number. It is a fatal
    error for a process to attempt to unlock a mutex which has not been locked
    by this process."""
    GA_Unlock(mutex)

def unpack(int g_src, int g_dst, int g_msk, lo=None, hi=None):
    """Expands the values in the source vector into a larger destination vector.

    The unpack subroutine is designed to expand the values in the source
    vector g_src into a larger destination array g_dst based on the values in
    an integer mask array g_msk. The values lo and hi denote the range of
    elements that should be uncompressed and icount is a variable that on
    output lists the number of values placed in the uncompressed array. This
    operation is the complement of the ga.pack operation. An example is shown
    below::

        ga.unpack(g_src, g_dst, g_msk, 1, n, &icount);
        g_src:   1  7  9 12 15 16
        g_msk:   1  0  0  0  0  0  1  0  1  0  0  1  0  0  1  1  0
        g_dst:   1  0  0  0  0  0  7  0  9  0  0 12  0  0 15 16  0
        icount:  6

    This is a collective operation.

    :Parameters:
        g_src : int
            handle for source arrray
        g_dst : int
            handle for destination array
        g_msk : int
            handle for integer array representing mask
        lo : 1D array-like of integers
            low value of range on which operation is performed
        hi : 1D array-like of integers
            hi value of range on which operation is performed

    """
    cdef np.ndarray[np.int64_t, ndim=1] hi_nd = inquire_dims(g_src)-1
    cdef int64_t c_lo=0, c_hi=hi_nd[0], icount
    if lo is not None:
        c_lo = lo
    if hi is not None:
        c_hi = hi
    GA_Unpack64(g_src, g_dst, g_msk, lo, hi, &icount)
    return icount

def update_ghosts(int g_a):
    """This call updates the ghost cell regions on each processor with the
    corresponding neighbor data from other processors.
    
    The operation assumes that all data is wrapped around using periodic
    boundary data so that ghost cell data that goes beyound an array boundary
    is wrapped around to the other end of the array. The ga.update_ghosts call
    contains two ga.sync calls before and after the actual update operation.
    For some applications these calls may be unecessary, if so they can be
    removed using the ga.mask_sync subroutine.

    This is a collective operation.

    """
    GA_Update_ghosts(g_a)

def update_ghost_dir(int g_a, int dimension, int dir, int flag):
    """This function can be used to update the ghost cells along individual
    directions. It is designed for algorithms that can overlap updates with
    computation. The variable dimension indicates which coordinate direction
    is to be updated (e.g. dimension = 1 would correspond to the y axis in a
    two or three dimensional system), the variable idir can take the values
    +/-1 and indicates whether the side that is to be updated lies in the
    positive or negative direction, and cflag indicates whether or not the
    corners on the side being updated are to be included in the update. The
    following calls would be equivalent to a call to GA_Update_ghosts  for a
    2-dimensional system:
     

    status = NGA_Update_ghost_dir(g_a,0,-1,1);
    status = NGA_Update_ghost_dir(g_a,0,1,1);
    status = NGA_Update_ghost_dir(g_a,1,-1,0);
    status = NGA_Update_ghost_dir(g_a,1,1,0);

    The variable cflag is set equal to 1 (or non-zero) in the first two calls
    so that the corner ghost cells are update, it is set equal to 0 in the
    second two calls to avoid redundant updates of the corners. Note that
    updating the ghosts cells using several independent calls to the
    nga_update_ghost_dir functions is generally not as efficient as using
    GA_Update_ghosts  unless the individual calls can be effectively
    overlapped with computation.

    """
    NGA_Update_ghost_dir(g_a, dimension, dir, flag)

def uses_ma():
    """TODO"""
    if GA_Uses_ma() == 1:
        return True
    return False

def wtime():
    """This function return a wall (or elapsed) time on the calling processor.
    Returns time in seconds representing elapsed wall-clock time since an
    arbitrary time in the past. Example:

    starttime = ga.wtime()
    # .... code snippet to be timed ....
    endtime   = ga.wtime()
    print "Time taken = %s seconds" % endtime-starttime

    This is a local operation.

    This function is only available in release 4.1 or greater.

    """
    return GA_Wtime()

def zero(int g_a, lo=None, hi=None):
    """Set all the elements in the array or patch to zero."""
    cdef np.ndarray[np.int64_t, ndim=1] lo_nd, hi_nd
    if lo is None and hi is None:
        GA_Zero(g_a)
    else:
        lo_nd,hi_nd = _lohi(g_a,lo,hi)
        NGA_Zero_patch64(g_a, <int64_t*>lo_nd.data, <int64_t*>hi_nd.data)

def zero_diagonal(int g_a):
    """Sets the diagonal elements of this matrix g_a with zeros.
    
    This is a collective operation. 

    """
    GA_Zero_diagonal(g_a)

initialize()
