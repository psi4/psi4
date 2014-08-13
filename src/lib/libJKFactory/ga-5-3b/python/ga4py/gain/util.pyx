# cython: profile=True
"""Contains index- and slice-related operations needed for bookkeeping.

Probably the most import set of classes are the Key-related ones including
Key, NoneKey, FixedKey, RangeKey, and MasterKey. These classes perform all of
the index arithmetic and index translation to make GAiN work.

The data distribution provided by Global Arrays is static. For a distributed
implementation of NumPy, this poses additional challenges on top of the
challenge of no longer having a single address space. The Key class contains
an attribute 'origin' which refers back to the GA dimension it was originally
associated with. This allows us to determine transpose operations.

"""
import numpy as np

cpdef list listify(thing):
    try:
        return list(thing)
    except:
        return [thing]

class Key(object):
    """An abstract base class for an index into an ndarray.

    You can index an array using an integer, a slice, or a None which are
    represented by the subclasses FixedKey, RangeKey, and NoneKey,
    respectively.
    
    Key instances may refer to an original dimension index using the 'origin'
    attribute. The 'origin' attribute only makes sense for the FixedKey and
    RangeKey subclasses.

    Subclasses should define __eq__, __getitem__, size, and pyobj.

    Attributes
    ----------
    origin : int, or None
        The original index this Key refers to.
    size : int
        shorthand for get_size() method
    
    """
    def __init__(self):
        """Create a new Key with origin=None."""
        self.origin = None

    def __eq__(self, other):
        """x.__eq__(y) <==> x==y

        Raises
        ------
        NotImplementedError :
            subclasses must define this

        """
        raise NotImplementedError, "subclasses must define this"

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]
        
        The given key is assumed to be a Python object.

        Raises
        ------
        NotImplementedError :
            subclasses must define this

        """
        raise NotImplementedError, "subclasses must define this"

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "Key()"

    def bound_by_lohi(self, lo, hi):
        """Return a new Key modified to fit between lo and hi.

        Raises
        ------
        NotImplementedError :
            subclasses must define this

        """
        raise NotImplementedError, "subclasses must define this"

    def get_size(self):
        """Returns number of elements represented by this Key.
        
        See Also
        --------
        size : Attribute getter for get_size method.
        
        """
        raise NotImplementedError, "subclasses must define this"
    size = property(get_size)

    def pyobj(self):
        """Returns a Python object based on this Key e.g. slice, int, None.
        
        Raises
        ------
        NotImplementedError :
            subclasses must define this
            
        """
        raise NotImplementedError, "subclasses must define this"

class NoneKey(Key):
    """A None with additional Key methods.
    
    Slicing a None is allowed. A None (newaxis) can be removed by slicing it
    with a value of 0. Slicing a None returns None which indicates the
    dimension should be removed.

    NoneKey instances should never have their origin attribute set since a
    None for an axis should never refer back to an original GA shape.
    
    """
    def __init__(self):
        Key.__init__(self)

    def __eq__(self, other):
        """Returns True of other is None or a NoneKey."""
        return isinstance(other, (type(None),NoneKey))

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]
        
        An index of 0 is allowed, or a slice starting from 0.
        
        Returns None if an index of 0 is passed because this effectively
        removes this Key from an index.
        
        """
        if isinstance(key, (RangeKey,slice)):
            if isinstance(key, slice):
                start,stop,step = key.indices(self.size)
            else:
                start,stop,step = key.start,key.stop,key.step
            if start != 0 or stop != 1 or step != 1:
                raise IndexError, "bad slice for NoneKey: %s" % key
            return self
        else:
            key = long(key)
            if key not in [0,0L]: # probably don't need '0' due to cast
                raise IndexError, "index out of bounds"
            return None # on purpose not self -- indicates index removal

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "NoneKey()"

    def bound_by_lohi(self, lo, hi):
        """Always raises NotImplementedError."""
        raise NotImplementedError, "bound_by_lohi nonsensical for NoneKey"

    def get_size(self):
        """Always returns 1."""
        return 1

    size = property(get_size)

    def pyobj(self):
        """Returns a Python object based on this Key (always None)."""
        return None

class FixedKey(Key):
    """An integer with additional Key methods.
    
    There must be a valid 'value' at all times. Anything that can be
    transformed into an int or long is allowed. long is prefered.

    FixedKey instances must always refer to an original GA dimension using the
    'origin' attribute. It does not make sense for a FixedKey to be without an
    'origin'.

    Attributes
    ----------
    origin : int, or None
        The original index this FixedKey refers to.
    size : int
        shorthand for get_size() method
    
    """
    def __init__(self, value, origin):
        """A integer value and an integer origin must be given.

        The given 'value' must be >= 0.
        The given 'origin' must be >= 0.

        """
        Key.__init__(self)
        self.value = long(value)
        self.origin = origin
        assert self.value >= 0
        assert self.origin >= 0

    def __eq__(self, other):
        """Compares the 'value' attributes."""
        if hasattr(other,'value'):
            return self.value == other.value
        return False

    def __getitem__(self, key):
        """You cannot index a FixedKey.
        
        Raises
        ------
        IndexError : 
            cannot index a FixedKey
            
        """
        raise IndexError, "cannot index a FixedKey"

    def __int__(self):
        return int(self.value)

    def __long__(self):
        return long(self.value)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "FixedKey(%s)" % self.value

    def bound_by_lohi(self, lo, hi):
        """Returns self but raises IndexError if value is not within range.
        
        Raises
        ------
        IndexError : 
            if not (lo <= self.value < hi)
            
        """
        if not (lo <= self.value < hi):
            raise IndexError, "lo/hi out of bounds %s <= %s < %s" % (
                    lo, self.value, hi)
        return self

    def get_size(self):
        """Always returns 1 even though 0 is more correct.

        FixedKey instances don't count towards the shape of an array because
        FixedKey instances remove dimensions. However, since it is often
        convenient to calculate the shape of multiple Key instances using a
        reduction of the sizes of the Key instances, we return 1 here instead
        of 0.

        """
        return 1
    size = property(get_size)

    def pyobj(self):
        """Returns a Python object based on this Key (long 'value' attr)."""
        return self.value

class RangeKey(Key):
    """A Python slice (w/o the indices function) with additional Key methods.
    
    RangeKey instances must always refer to an original GA dimension using the
    'origin' attribute. It does not make sense for a RangeKey to be without an
    'origin'.


    """
    def __init__(self, start, stop=None, step=None, origin=None):
        """Creates given either a slice or an explicit start, stop, and step.

        Asserts that the start, stop, and step are not None.
        An 'origin' is required. The given 'origin' must be >= 0.

        """
        Key.__init__(self)
        if isinstance(start, slice):
            self.start = start.start
            self.stop = start.stop
            self.step = start.step
            assert stop is None
            assert step is None
        else:
            self.start = start
            self.stop = stop
            self.step = step
        assert start is not None
        assert stop is not None
        assert step is not None
        assert origin is not None
        assert origin >= 0
        self.origin = origin

    def __eq__(self, other):
        """Returns true if all start, stop, and step values match other."""
        try:
            return (self.start == other.start
                    and self.stop == other.stop
                    and self.step == other.step)
        except AttributeError:
            return False

    def __getitem__(self, key):
        """You can index using either a slice, RangeKey, or integer."""
        if isinstance(key, (RangeKey,slice)):
            if isinstance(key, slice):
                _start,_stop,_step = key.indices(self.size)
            else:
                _start,_stop,_step = key.start,key.stop,key.step
            start = ((_start*self.step) + self.start)
            stop = ((_stop*self.step) + self.start)
            step = _step * self.step
            return RangeKey(start,stop,step,origin=self.origin)
        else:
            try:
                key = long(key)
            except TypeError:
                raise TypeError, ("long() argument must be a string or a "
                        "number, not %s" % key)
            if key < 0:
                key += self.size
            if key >= self.size or key < 0:
                raise IndexError, "invalid index"
            shifted = (key*self.step) + self.start
            if (self.step < 0 and shifted <= self.stop
                    or self.step > 0 and shifted >= self.stop):
                raise IndexError
            return FixedKey(shifted, self.origin)

    def __str__(self):
        return "RangeKey(%s,%s,%s)" % (self.start,self.stop,self.step)

    def __repr__(self):
        return str(self)

    def bound_by_lohi(self, lo, hi):
        """Return a new RangeKey modified to fit between lo and hi.

        Raises
        ------
        IndexError if start/stop/step is out of the bounds of lo/hi .

        Returns
        -------
        A new RangeKey between the values of lo and hi.

        Examples
        --------
        >>> RangeKey(1,10,2).bound_by_lohi(4,20)
        RangeKey(5,10,2)
        >>> RangeKey(1,10,2).bound_by_lohi(1,20)
        RangeKey(1,10,2)
        >>> RangeKey(1,10,2).bound_by_lohi(4,10)
        RangeKey(5,10,2)
        >>> RangeKey(1,10,2).bound_by_lohi(4,9)
        RangeKey(5,9,2)

        """
        new_start = 0
        new_stop = 0
        if self.step > 0:
            if self.start >= hi:
                raise IndexError, "start >= hi (out of bounds)"
            elif self.start >= lo: # start < hi is implied
                new_start = self.start
            else: # start < lo < hi is implied
                guess = (lo-self.start)//self.step
                new_start = guess*self.step + self.start
                while new_start < lo:
                    guess += 1
                    new_start = guess*self.step + self.start
            if self.stop <= lo:
                raise IndexError, "stop <= lo (out of bounds)"
            elif self.stop <= hi: # lo < stop is implied
                new_stop = self.stop
            else: # lo < hi < stop is implied
                new_stop = hi # this should be good enough
        else:
            if self.start < lo:
                raise IndexError, "negative step, start < lo (out of bounds)"
            elif self.start < hi:
                new_start = self.start
            else: # start >= hi >= lo
                guess = (hi-self.start)//self.step
                new_start = guess*self.step + self.start
                while new_start >= hi:
                    guess += 1
                    new_start = guess*self.step + self.start
            if self.stop >= hi:
                raise IndexError, "negative step, stop >= hi (out of bounds)"
            elif self.stop >= (lo-1):
                new_stop = self.stop
            else:
                new_stop = lo-1 # this should be good enough
        result = RangeKey(new_start,new_stop,self.step,self.origin)
        if result.size <= 0:
            raise IndexError, "bound_by_lohi resulted in 0 length"
        return result

    def get_size(self):
        """Returns the length of the range."""
        start,stop,step = self.start,self.stop,self.step # for brevity
        if (step < 0 and stop >= start) or (step > 0 and start >= stop):
            return 0
        elif step < 0:
            return (stop - start + 1) / (step) + 1
        else:
            return (stop - start - 1) / (step) + 1
    size = property(get_size)

    def pyobj(self):
        """Returns a Python object based on this Key (a slice)."""
        return slice(self.start,self.stop,self.step)

class MasterKey(object):
    """A list of Key instances with some extra methods.

    Represents one or more indices into an array. Contains methods to
    calculate shape, constrain the indices to lo/hi bounds of various types.

    Also contains attributes T and TT to represent axes to transpose as well
    as the inverse of the axes so that the transpose can be undone. The T and
    TT attributes will always be set in tandem -- both will be None or both
    will be set to a tuple of integers. Certain methods will query the T and
    TT attributes and change their behavior appropriately. Slicing a
    transposed MasterKey is one such example.

    """
    def __init__(self, shape=None, data=None, fixed=None, T=None, TT=None):
        if shape is None and data is None:
            raise ValueError, "specify either shape or data"
        elif shape is not None:
            self.data = [RangeKey(0,x,1,origin=i)
                    for i,x in enumerate(shape)]
        elif data is not None:
            for datum in data:
                assert isinstance(datum, Key)
            self.data = data
        else:
            raise ValueError, "specify either shape or data"
        self.T = T
        self.TT = TT
        if fixed is None:
            self.fixed = []
        else:
            self.fixed = fixed

    def __eq__(self, other):
        try:
            return (self.data == other.data) and (self.fixed == other.fixed)
        except AttributeError:
            return False

    def __iter__(self):
        return iter(self.data)

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return repr(self.data)

    def _count_instances(self, *args):
        """Return the number of instances of the given type."""
        return len(self._get_instances(args))

    def _get_instances(self, *args):
        return [item for item in self.data if isinstance(item, args)]

    def _get_origins_within_data(self):
        """Return list of keys which have their origin set."""
        return [item for item in self.data if item.origin is not None]

    def get_origin(self):
        """Return current list of FixedKey and RangeKey instances.

        The returned list is in the current (possibly transposed) order.

        See Also
        --------
        get_sorted_origin() : sorts the list by 'origin' property

        """
        origins = self._get_origins_within_data()
        result = [None]*(len(self.fixed)+len(origins))
        for fixed in self.fixed:
            result[fixed.origin] = fixed
        iter_origins = iter(origins)
        for i in range(len(result)):
            if result[i] is None:
                result[i] = iter_origins.next()
        return result

    def get_sorted_origin(self):
        """Return sorted list of FixedKey and RangeKey instances.

        The returned list is in the current (possibly transposed) order.

        See Also
        --------
        get_origin() : does not sort the list by 'origin' property

        """
        origins = self._get_origins_within_data()
        result = [None]*(len(self.fixed)+len(origins))
        for fixed in self.fixed:
            result[fixed.origin] = fixed
        for datum in origins:
            result[datum.origin] = datum
        return result

    def get_original_ndim(self):
        return len(self.get_origin())

    def get_shape(self):
        return tuple(
                [item.size for item in self._get_instances(NoneKey,RangeKey)])
    shape = property(get_shape)

    def get_size(self):
        return reduce(lambda x,y: x*y, self.shape, 1)
    size = property(get_size)

    def get_ndim(self):
        return self._count_instances(NoneKey,RangeKey)
    ndim = property(get_ndim)

    def pyobj(self):
        return [x.pyobj() for x in self.data]

    def replace_ellipsis(self, key):
        """Given key, replace one or more Ellipsis based on this MasterKey.

        Assumes key is a list, and if not, makes it a list.
        Raises IndexError if the key is too large.

        """
        if type(key) != type([]):
            key = [key]
        ndim = self.ndim
        count_real_keys = len(key)-key.count(None)-key.count(Ellipsis)
        if count_real_keys > ndim:
            raise IndexError, "invalid index %s[%s]" % (self,key)
        # implicit Ellipsis at end of key if one wasn't given
        if Ellipsis not in key:
            key.append(Ellipsis)
        # first Ellipsis replaced with as many slice(None,None,None) as needed
        ellipsis_index = key.index(Ellipsis)
        ellipsis_count = key.count(Ellipsis)
        key = (key[:ellipsis_index]
                + ([slice(None,None,None)] * (ndim-count_real_keys))
                + key[ellipsis_index+1:])
        # remove all remaining Ellipsis from key
        ellipsis_count -= 1
        while ellipsis_count > 0:
            key.remove(Ellipsis)
            ellipsis_count -= 1
        return key

    def __getitem__(self, key):
        """Slices this MasterKey and returns a new MasterKey."""
        key = self.replace_ellipsis(key)
        new_data = []
        new_fixed = self.fixed[:] # copy fixed dimensions
        key_iter = iter(key)
        for item in self.data:
            if isinstance(item, (RangeKey,NoneKey)):
                next_key = key_iter.next()
                while next_key is None:
                    new_data.append(NoneKey())
                    next_key = key_iter.next()
                result = item[next_key]
                if isinstance(item, RangeKey) and isinstance(result, FixedKey):
                    # went from range to fixed, add to fixed keys
                    new_fixed.append(result)
                elif result is not None:
                    new_data.append(result)
            elif isinstance(item, FixedKey):
                raise TypeError, "FixedKey found in MasterKey"
            else:
                raise TypeError, "unrecognized item in MasterKey"
        # there might be some items left in the key_iter
        # the only valid items that can remain are None
        for next_key in key_iter:
            assert next_key is None
            new_data.append(NoneKey())
        result = MasterKey(None, new_data, new_fixed, self.T, self.TT)
        return result

    def bound_by_lohi(self, lo, hi):
        """Calculate new MasterKey for the subarray denoted by lo and hi.
        
        returns a list of slice objects representing the piece of the subarray
        that the lo and hi bounds maintain after applying bounds to
        global_slice.

        raises IndexError if the lo and hi bounds don't get a piece of the
        subarray
        
        Examples
        --------
        >>> j = MasterKey((10,15,20,25))
        >>> k = j[1:10:2, 10:1:-3, 2, None, 10]
        >>> k.bound_by_lohi([4,4,0,5], [10,9,3,20])
        [RangeKey(5,10,2), RangeKey(7,3,-3)]

        """
        result = []
        lo = [l for l in lo]
        hi = [h for h in hi]
        assert len(lo) == len(hi) == self.get_original_ndim()
        # verify fixed dimensions within bounds; raises IndexError otherwise
        for item in self.fixed:
            item.bound_by_lohi(lo[item.origin], hi[item.origin])
        for item in self.data:
            if isinstance(item, RangeKey):
                result.append(
                        item.bound_by_lohi(lo[item.origin], hi[item.origin]))
            elif isinstance(item, NoneKey):
                pass
            elif isinstance(item, FixedKey):
                raise TypeError, "FixedKey found in MasterKey"
            else:
                raise TypeError, "unrecognized item in MasterKey"
        return result

    def lohi_T(self):
        """Return the current order of origin including fixed dimensions.

        This produces a tranpose suitable for a lo or hi.

        Examples
        --------
        >>> lo,hi = [0,0,0], [2,3,4]
        >>> k = MasterKey((2,3,4))
        >>> k_T = k.transpose((2,1,0))
        >>> k_T1 = k_T[:,1]
        >>> k_T1.lohi_T()

        """
        return [item.origin for item in self.get_origin()]

    def None_key(self):
        """Return a key based on current None's in this MasterKey.

        Examples
        --------
        >>> k = MasterKey((2,3,4))
        >>> k1 = k[None, :, None, :, None, :, None]
        >>> k1.None_key()
        [None, slice(None,None,None), 

        """
        return [None if isinstance(datum, NoneKey) else slice(None,None,None)
                for datum in self.data]

    def get_lo_hi_skip_adjust(self):
        """Return lo/hi/skip/adjust needed for successful ga.strided_get()."""
        sorted_origin = self.get_sorted_origin()
        lo = []
        hi = []
        sk = []
        ad = []
        need_strided = False
        for item in sorted_origin:
            if isinstance(item, FixedKey):
                lo.append(item.value)
                hi.append(item.value+1)
                sk.append(1)
                ad.append(slice(0,None,None))
            elif isinstance(item, RangeKey):
                if item.step > 1 or item.step < -1:
                    need_strided = True
                if item.step < 0:
                    length = item.size-1
                    lo.append(item.step*length + item.start)
                    hi.append(item.start+1)
                    sk.append(-item.step)
                    ad.append(slice(None,None,-1))
                else:
                    lo.append(item.start)
                    hi.append(item.stop)
                    sk.append(item.step)
                    ad.append(slice(0,None,None))
            else:
                raise TypeError, "unhandled piece of MasterKey"
        lo = np.asarray(lo, dtype=np.int64)
        hi = np.asarray(hi, dtype=np.int64)
        sk = np.asarray(sk, dtype=np.int64)
        return lo,hi,sk,ad,need_strided

    def access_key(self, lo, hi):
        """Converts from MasterKey to list appropriate for ga.access()[list].

        Does NOT return a new MasterKey or any Key subclasses. The returned
        list contains one or more of slice, int, long, etc which an
        access()ed ndarray can iterpret.

        """
        result = [None]*self.get_original_ndim()
        for fixed in self.fixed:
            result[fixed.origin] = fixed.value-lo[fixed.origin]
        iter_ranged_origin = iter(self.bound_by_lohi(lo, hi))
        for i in range(len(result)):
            if result[i] is None:
                ranged = iter_ranged_origin.next()
                nlo = lo[ranged.origin]
                start = ranged.start-nlo
                stop = ranged.stop-nlo
                if start < 0: start = None
                if stop < 0: stop = None
                result[i] = slice(start,stop,ranged.step)
        return result

    def get_key(self, lo, hi):
        """Converts from MasterKey to list appropriate for ga.get()[list]."""
        iter_ranged_origin = iter(self.bound_by_lohi(lo, hi))
        result = []
        for datum in self.data:
            if isinstance(datum, FixedKey):
                raise TypeError, "FixedKey found in MasterKey"
            elif isinstance(datum, RangeKey):
                ro = iter_ranged_origin.next()
                offset = RangeKey(
                        datum.start,ro.start,datum.step,datum.origin).size
                result.append(slice(offset, offset+ro.size, 1))
            elif isinstance(datum, NoneKey):
                result.append(slice(None,None,None))
            else:
                raise ValueError, "MasterKey contained unknown object"
        return result

    def transpose(self, axes):
        """Swap the axes of this MasterKey given axes.

        Returns the inverse of the axes in addition to the new MasterKey
        instance.

        Parameters
        ----------
        axes : list of ints
            `i` in the `j`-th place in the tuple means `a`'s `i`-th axis becomes
            `a.transpose()`'s `j`-th axis.
        
        Examples
        --------
        >>> k = MasterKey((40,30,20))
        >>> k = k[2:34:1, 4, None, 10:2:-1]
        >>> inv,new_k = k.transpose((1,2,0))
        >>> print inv
        [2, 0, 1]
        >>> print new_k

        """
        # normalize the inputs, making sure they're all iterable
        axes = listify(axes)
        # create the inverse of the given axes
        inverse = [None]*len(axes)
        for i,val in enumerate(axes):
            inverse[val] = i
        # modify the given axes based on this MasterKey
        new_data = [self.data[i] for i in axes]
        return MasterKey(None, new_data, self.fixed, axes, inverse)

def broadcast_shape(first, second):
    """Return the broadcasted version of shapes first and second."""
    def worker(x,y):
        if x is None: x = 1
        if y is None: y = 1
        if x != 1 and y != 1 and x != y:
            raise ValueError, ("shape mismatch:"
                    " objects cannot be broadcast to a single shape")
        return max(x,y)
    return tuple(reversed(map(worker, reversed(first), reversed(second))))

def broadcast_chomp(smaller_key, larger_key):
    """Return a key appropriate for the given shape."""
    new_key = []
    for s,l in zip(reversed(smaller_key),reversed(larger_key)):
        if s == 1:
            new_key.append(slice(0,1,1))
        else:
            new_key.append(l)
    new_key.reverse()
    return new_key
    
def unravel_index(x,dims):
    """Like np.unravel_index, but 'x' can be an integer array.

    Yeah, I know, numpy 1.6.0 has this already, but we're based on 1.5.1.
    I copied the code and modified it from 1.5.1.

    """
    import numpy as _nx
    x = _nx.asarray(x)
    if x.ndim == 0:
        return _nx.unravel_index(x,dims)
    max = _nx.prod(dims)-1
    if _nx.any(x>max) or _nx.any(x<0):
        raise ValueError("Invalid index, must be 0 <= x <= number of elements.")
    idx = _nx.empty_like(dims)

    # Take dimensions
    # [a,b,c,d]
    # Reverse and drop first element
    # [d,c,b]
    # Prepend [1]
    # [1,d,c,b]
    # Calculate cumulative product
    # [1,d,dc,dcb]
    # Reverse
    # [dcb,dc,d,1]
    dim_prod = _nx.cumprod([1] + list(dims)[:0:-1])[::-1]
    # Indices become [x/dcb % a, x/dc % b, x/d % c, x/1 % d]
    return tuple(x[:,None]//dim_prod % dims)
