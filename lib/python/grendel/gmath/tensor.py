from collections import Iterable, Callable
from copy import copy
import inspect
from itertools import product, islice, permutations
import math
import traceback
import types
from warnings import warn

import numpy as np
import operator

from grendel import sanity_checking_enabled
import sys

import grendel
from grendel.gmath.misc import permutation_is_even
from grendel.util.abstract_bases import MethodLike, FunctionLike
from grendel.util.decorators import typechecked
from grendel.util.overloading import get_kwarg, pop_kwarg, listify_args, pop_multikwarg
from grendel.util.strings import classname
#noinspection PyUnresolvedReferences
from grendel.util.units import isunit, IncompatibleUnitsError, hasunits, strip_units

__all__ = ["Tensor", "chop", "chopped", 'ComputableTensor']

def chop(tens, thresh=1e-14):
    """ Round off near-integer values (within `thresh`) to integers for each value in the Tensor.
    """
    for num in np.nditer(tens, op_flags=['readwrite']):
        if abs(num) < thresh: num[...] = 0.0
        if abs(num - round(num)) < thresh: num[...] = round(num)

def chopped(tens, thresh=1e-14):
    """ Same as `chop`, but does not modify the original and instead returns a copy.
    """
    it = np.nditer([tens, None])
    for src, dest in it:
        val = round(src) if abs(it[0] - round(it[0])) < thresh else src
        dest[...] = 0.0 if val == 0.0 else val
    return it.operands[1].view(type(tens))


class LightTensor(np.ndarray):
    """ A super-lightweight version of Tensor that doesn't have things like indices, units, etc.
    Primarily for internal and advanced use.  Light tensor implements "light" versions of some of the
    methods in `Tensor` that return a LightTensor rather than a Tensor.  These are labeled as "l_<method_name>"
    """

    def __new__(cls, iterable=None, **kwargs):
        if iterable:
            return np.array(iterable, **kwargs).view(cls)
        else:
            return np.array(**kwargs).view(cls)


    ###########
    # Methods #
    ###########


    #-----------------#
    # Inquiry methods #
    #-----------------#

    def is_zero(self, cutoff=None):
        """ Returns True if all elements of `self` have absolute values less than `cutoff`, which defaults
        to `Tensor.zero_cutoff`
        .. note::
           Since `zero_cutoff` is a pseudo-class attribute in Tensor, not LightTensor, individual LightTensor
           instances cannot set their own `zero_cutoff` attribute.
        """
        cutoff = cutoff if cutoff is not None else Tensor.zero_cutoff
        return (abs(self) < cutoff).all()

    def same_upto_phase_factor(self, *others, **kwargs):
        if sanity_checking_enabled:
            if not np.isrealobj(self) or any(not np.isrealobj(o) for o in others):
                raise NotImplementedError("phase factor detection complex Tensors is not yet implemented")
        cutoff = kwargs.pop('cutoff', Tensor.same_tensor_cutoff)
        for other in others:
            if (abs(other - self) > cutoff).any() and (abs(other + self) > cutoff).any():
                return False
        return True

    #------------------------------#
    # Methods that return a scalar #
    #------------------------------#

    def norm(self):
        """ Returns the square root of the sum of the element-wise product of `self` and `other`
        """
        return math.sqrt(np.einsum(self, range(len(self.shape)), self, range(len(self.shape)), []))

    def zero(self):
        """ Zeros all entries in the tensor
        """
        self.fill(0)

    #-----------#
    # Iterators #
    #-----------#

    def iter_vectors(self, with_indices=False):
        last_dim = self.shape[-1]
        flags = []
        if with_indices:
            flags.append('multi_index')
        it = np.nditer(self, flags=flags)
        curr_ary = []
        # TODO make this yield a view (?)
        for x in it:
            curr_ary.append(x)
            if len(curr_ary) == last_dim:
                if with_indices:
                    yield LightVector(curr_ary), x.multi_index[:-1]
                else:
                    yield LightVector(curr_ary)
                curr_ary = []

    def iter_with_indices(self):
        it = np.nditer(self, flags=['multi_index'])
        for val in it:
            yield val, it.multi_index

    #------------------------------------#
    # "Light" versions of Tensor methods #
    #------------------------------------#

    def l_flatten(self, order = 'C'):
        return super(LightTensor, self).flatten(order).view(LightVector)



# TODO transpose should transpose indices
# TODO different units for different entries in the Tensor (implemented by maintaining a sister array of units)
class Tensor(LightTensor):
    """
    Mostly a wrapper for the NumPy ndarray class.  This gives us a place to put tensor-related functionality that we
    need but is not available in NumPy.

    """

    ####################
    # Class Attributes #
    ####################

    same_tensor_cutoff = 1e-8

    # Used by is_zero to determine if all of the elements in the Tensor are 0.
    # note that elements smaller than this are still stored and otherwise
    # treated normally.
    zero_cutoff = 1e-8

    ##############
    # Attributes #
    ##############

    indices = None
    """ The einstein summation indices over which the tensor is defined.  Must be already 'declared' in to be part
     of some IndexRange object
    """

    index_range_set = None
    """
    The IndexRangeSet object that the indices are drawn from.
    """

    units = None
    """ If the entire tensor has the same units, then they can be stored here.
    Different units for individual elements is not yet implemented.
    """

    name = None
    """ An optional name for the tensor, used in output, etc.
    """

    ##################
    # Initialization #
    ##################

    def __new__(cls, *args, **kwargs):
        """
        TODO Move this to class level documentation, since it doesn't show up in Sphinx
        Takes several different forms.  Given a number of arguments (possibly nested lists), a Tensor is created as
        expected.  *This is the only form in which arguments may be specified without keywords.*
        `Tensor` can also be initialized by giving a list of dimensions for the `shape` (aliased as `dimension`)
        keyword argument, with a possible default value keyword argument `default_val`.
        Unless otherwise specified via the `dtype` keyword argument, it is assumed that the input data should be cast
        as a `numpy.float64`.

        """
        # Pop off any special args or kwargs and store their values for later
        has_indices = False
        ret_val = None
        indices = None
        units = pop_kwarg(kwargs, 'units')
        name = kwargs.pop('name', None)
        #--------------------------------------------------------------------------------#
        # pop off any kwargs that the subclass's __init__ takes...
        subclass_kwargs = {}
        for supercls in cls.__mro__:
            if supercls is Tensor:
                break
            if hasattr(supercls, '__tensor_init__') and callable(supercls.__tensor_init__):
                if hasattr(supercls.__tensor_init__, 'getargspec'):
                    argspec = supercls.__tensor_init__.getargspec()
                else:
                    argspec = inspect.getargspec(supercls.__tensor_init__)
                for arg in argspec.args[1:]:
                    subclass_kwargs[arg] = kwargs.pop(arg, None)
        #--------------------------------------------------------------------------------#
        # Check for indices...
        indices_kwarg = kwargs.pop('indices', None)
        if indices_kwarg is None and len(args) == 1 and isinstance(args[0], basestring):
            args = list(args)
            indices_kwarg = args.pop(0)
            args = tuple(args)
        index_range_set = pop_multikwarg(kwargs, 'index_range_set', 'in_set', 'set')
        if indices_kwarg is not None:
            has_indices = True
            indices = EinsumTensor.split_indices(indices_kwarg)
            if index_range_set is None:
                index_range_set = IndexRange.global_index_range_set
            shape = []
            for idx in indices:
                if idx not in index_range_set.known_ranges:
                    raise IndexError("unknown index '" + idx + "'")
                shape.append(index_range_set.known_ranges[idx].size)
            shape = tuple(shape)
            if "shape" in kwargs:
                kwshape = kwargs.pop('shape')
                if shape != kwshape:
                    raise TypeError("inconsistent shape:  indices '{}' indicate a shape of {}, but"
                                    " the keyword 'shape' was given with the shape {}".format(
                        ",".join(indices), shape, kwshape
                    ))
            kwargs['shape'] = shape
        #--------------------------------------------------------------------------------#
        # Now create a numpy.ndarray object...
        def_val = pop_kwarg(kwargs, 'default_val', 'default_value', 'default') or 0.0
        # Set the default data type to float, unless the user specifies
        dtype = get_kwarg(kwargs, 'dtype') or float
        if not callable(dtype) and grendel.show_warnings:
            warn("dtype given to {0} constructor is not Callable and thus cannot be used for casting." \
                 "  It is better to use callable types for dtype if possible; e.g. numpy.float64 instead" \
                 "of numpy.dtype('float64')".format(classname(cls)))
        kwargs['dtype'] = dtype
        # Typecast the default value
        if callable(dtype):
            def_val = dtype(def_val)
        # See if we have a shape...
        shape = pop_kwarg(kwargs, 'shape', 'dimension')
        # See if we have data...
        # This allows us to support the form Tensor(1, 2, 3, 4)
        if len(args) == 1 and isinstance(args[0], np.ndarray):
            data = args[0]
        else:
            data = listify_args(*args)
        if 'data' in kwargs:
            if data:
                raise TypeError("`data` may be specified as a keyword argument or as " \
                                "the regular arguments to {0} constructor, but not both.".format(classname(cls)))
            else:
                data = pop_kwarg('data')
        if shape and not isinstance(data, np.ndarray):
            has_content = False
            if len(args) != 0:
                has_content = True
            if not has_content:
                if def_val == 0.0:
                    try:
                        ret_val = np.zeros(shape=shape, **kwargs)
                    except:
                        raise
                else:
                    ret_val = (np.ones(shape=shape, **kwargs) * def_val)
            else:
                if grendel.sanity_checking_enabled:
                    # Check data length
                    tmp = np.array(data)
                    needed_data_size = 1
                    for dim in shape: needed_data_size *= dim
                    try:
                        tmp.reshape((needed_data_size,))
                    except ValueError:
                        raise ValueError("Data provided to {0} constructor is incompatible with the shape {1}".format(classname(cls), shape))
                    # Check data type
                ret_val = np.array(data, **kwargs).reshape(shape)
        else:
            # Just pass on the data and any surviving kwargs
            try:
                if isinstance(data, np.ndarray) and (
                        len(kwargs) == 0
                     or (len(kwargs) == 1 and 'dtype' in kwargs and kwargs['dtype'] == data.dtype)
                ):
                    # Just do a view
                    ret_val = data.view(cls)
                else:
                    # Otherwise, we need to call the numpy "constructor" (actually a factory function) of ndarray
                    ret_val = np.array(data, **kwargs)
            except:
                # debugging breakpoint hook
                raise
            if shape and ret_val.shape != shape:
                raise ValueError("Shape mismatch: data shape {0} does not match specified shape {1}".format(
                    data.shape, shape
                ))
        #--------------------------------------------------------------------------------#
        # View-cast the ret_val to the class in question, but only if we haven't already
        if not isinstance(ret_val, cls):
            ret_val = ret_val.view(cls)
        # Now assign stuff from any special args...
        if has_indices:
            ret_val.indices = indices
            ret_val.index_range_set = index_range_set
        else:
            ret_val.indices = None
        ret_val.units = units
        ret_val.name = name
        if name is None:
            ret_val.name = "(unnamed tensor)"
        # pass the remaining kwargs to the initializer...
        ret_val.__tensor_init__(**subclass_kwargs)
        return ret_val

    def __tensor_init__(self, **kwargs):
        if len(kwargs) > 0:
            raise TypeError("Unknown keyword arguments to {} constructor: ({})".format(
                self.__class__.__name__,
                ', '.join(kwargs.keys())
            ))

    ##############
    # Properties #
    ##############

    @property
    def index_ranges(self):
        """ List of pointers to IndexRange objects over which the tensor is defined.
        """
        return [self.index_range_set[i] for i in self.indices] if self.indices else None

    @property
    def value(self):
        """ Allows conformance with the `Unitized` protocol.
        """
        ret = self.view(self.__class__)
        ret.units = None
        return ret

    @property
    def diagonal(self):
        if not self.is_square():
            raise NotImplementedError
        ndim = len(self.shape)
        return Vector([self[(i,)*ndim] for i in xrange(self.shape[0])])

    @property
    def I(self):
        if not len(self.shape) == 2:
            raise NotImplementedError
        else:
            return self.view(np.matrix).I.view(type(self))


    ###################
    # Special Methods #
    ###################

    #---------------------#
    # Container Emulation #
    #---------------------#

    # TODO: Mixed index letters and ints?
    # TODO: partial index specification (e.g. T['a','b'] returns a matrix if T is a 4-index tensor)
    def __getitem__(self, args):
        """
        Returns the sub-``Tensor`` corresponding to the depth specified.  If the resulting sub-``Tensor`` is just an item
        (i.e. if the ``ndim`` attribute of self is equal to the number of arguments given), the ``numpy.ndarray``
        behavior is used (which just returns the entry, which is a special numpy subclass of ``float`` or ``int``).
        If the resulting sub-``Tensor`` has a ``ndim`` attribute of 1, a ``Vector`` object is returned.
        If the resulting sub-``Tensor`` has a ``ndim`` attribute of 2, a ``Matrix`` object is returned.

        :Examples:


        TODO:  Write example test cases
        TODO: Move this to the class documentation since it doesn't show up in sphinx

        """
        # if a single item is given, handle it as a length-1 tuple
        args = (args,) if not isinstance(args, tuple) else args
        #----------------------------------------#
        # Handle Einstein summation
        if all(isinstance(arg, basestring) for arg in args):
            argidxs = EinsumTensor.split_indices(args)
            return EinsumTensor(argidxs, self)
        #----------------------------------------#
        else:
            try:
                ret_val = np.ndarray.__getitem__(self, args)
            except:
                raise
            if np.isscalar(ret_val):
                return ret_val
            else:
                shp = ret_val.shape
                if len(shp) == 1:
                    return ret_val.view(Vector)
                elif len(shp) == 2:
                    return ret_val.view(Matrix)
                else:
                    return ret_val

    def __setitem__(self, key, value):
        key = (key,) if not isinstance(key, tuple) else key
        if all(isinstance(arg, basestring) for arg in key):
            argidxs = EinsumTensor.split_indices(key)
            #----------------------------------------#
            if isinstance(value, EinsumContraction):
                # Carry out the contraction and store the result...
                # TODO Sanity checking for correct shape
                value.contract(dest=EinsumTensor(argidxs, tensor=self))
            elif isinstance(value, EinsumTensor):
                idxs, subidxs = EinsumTensor.split_indices(argidxs, include_sub=True)
                eself = EinsumTensor(argidxs, tensor=self)
                if idxs == value.indices and subidxs == value.sub_indices:
                    if value.coeff == 1.0:
                        super(Tensor, eself.sliced_tensor).__setitem__(Ellipsis, value.sliced_tensor)
                    else:
                        super(Tensor, eself.sliced_tensor).__setitem__(Ellipsis, value.coeff * value.sliced_tensor)
                    #np.ndarray.__setitem__(self, Ellipsis, value.tensor)
                elif len(idxs) == len(value.indices):
                    # Just rearrange things...
                    rv = value.sort_to(dest=eself, multiplier=value.coeff)
                else:
                    # We're doing an internal contraction.  EinsumSum.sum_into handles this.
                    tmp = EinsumSum(value)
                    tmp.sum_into(eself)
            elif isinstance(value, EinsumSum):
                # TODO Sanity checking for correct shape
                value.sum_into(dest=EinsumTensor(argidxs, tensor=self))
            elif isinstance(value, Tensor):
                if grendel.sanity_checking_enabled:
                    # Check the shape of the tensor to be assigned to the block
                    # TODO Sanity checking for correct shape
                    pass
                #np.ndarray.__setitem__(dest_tens, tuple([Ellipsis]*self.ndim), value)
                # The user probably made a mistake
                raise ValueError
        else:
            super(Tensor, self).__setitem__(key, value)

    def __contains__(self, item):
        return super(Tensor, self).__contains__(item)

    #----------------------#
    # Copying and Pickling #
    #----------------------#

    def __copy__(self):
        return self.__class__(
            # Just the data...
            copy(self.view(np.ndarray)),
            # and the keyword arguments
            **self.__copy_kwargs__()
        )

    def __copy_kwargs__(self):
        return dict(
            indices=copy(self.indices),
            shape=self.shape,
            index_range_set=self.index_range_set,
            units=self.units,
            name=self.name
        )

    #-------------------------------#
    # Numpy ndarray "magic" methods #
    #-------------------------------#

    __array_finalize_attributes__ = ['units', 'indices', 'index_range_set']

    def __array_finalize__(self, obj):
        """
        See http://docs.scipy.org/doc/numpy/user/basics.subclassing.html#array-finalize for explanation
        Initialize any defaults that may need to be initialized in a case where the __new__ is not called explicitly
        (but remember that this also gets called when __new__ is called explicitly)
        Remember that both view casting (to Tensor) and new-from-template must be handled here.
        """
        # remember self is the *new* instance and obj is the old one!
        if obj is None:
            # The case of the call from __new__:  Nothing to do here...
            return
        # Now handle the other cases.  This is more difficult.
        #========================================#
        # Here we specify a general framework for easily inheriting from Tensor
        #   First, set the __array_finalize_attributes__ attribute of the subclass to
        #     a list of identifiers that are to be copied.
        #   Then, if anything special needs to be done for a given attribute, define
        #     a method named __array_finalize_{attribute}__ that takes the same
        #     arguments as __array_finalize__ and returns the new value to be assigned
        #     (it does _not_ need to do the assignment itself).  If all that needs to be
        #     done is the specification of a default value, you can set the value of
        #     __array_finalize_{attribute}__ to the default value (which defaults to None).
        #   Finally, if any sanity checking needs to be done, define a method named
        #     __array_finalize_{attribute}_sanity_check__ which takes the same arguments
        #     as __array_finalize__ and has a return value that is ignored.
        # TODO Benchmark this.  I think it is slowing things down a lot
        attributes_to_transfer = {}
        for cls in self.__class__.__mro__:
            if '__array_finalize_attributes__' in cls.__dict__:
                for attr in cls.__array_finalize_attributes__:
                    finalize_func_or_value = getattr(self, '__array_finalize_' + attr + "__", None)
                    if callable(finalize_func_or_value):
                        attributes_to_transfer[attr] = finalize_func_or_value(obj)
                    else:
                        attributes_to_transfer[attr] = getattr(obj, attr, finalize_func_or_value)
                    if sanity_checking_enabled:
                        sanity_func = getattr(
                            self,
                            '__array_finalize_' + attr + "_sanity_check__",
                            None)
                        if callable(sanity_func):
                            sanity_func(obj)
            if cls is Tensor:
                break
        for attr, val in attributes_to_transfer.iteritems():
            setattr(self, attr, val)
        #========================================#

    def __array_finalize_indices__(self, obj):
        # And the indices... for now only pass them on if we self and obj are the same exact shape.
        # It's possible that we could do more here...see numpy.ndarray.base
        if isinstance(obj, Tensor):
            if obj.indices:
                if self.shape == obj.shape:
                    return copy(obj.indices)
        return None

    def __array_wrap__(self, out_arr, context=None):
        """
        See http://docs.scipy.org/doc/numpy/user/basics.subclassing.html#array-wrap for explanation
        """
        if context is not None:
            if context[0] is np.multiply:
                if len(context[1]) == 2:
                    op1, op2 = context[1]
                else:
                    _, op1, op2 = context[1]
                if op1 is self:
                    tens = op1
                    other = op2
                elif op2 is self:
                    tens = op2
                    other = op1
                else:
                    raise Exception("something very enigmatic happened...")
                #----------------------------------------#
                if hasunits(other):
                    out_arr.units = self.__mul_units__(tens.units, other.units)
                else:
                    out_arr.units = self.units
        #--------------------------------------------------------------------------------#
        return np.ndarray.__array_wrap__(self, out_arr, context)

    #----------------------#
    # Comparison operators #
    #----------------------#

    def __eq__(self, other):
        """ True if and only if the tensors have the same shape and are element-wise equal.
        """
        if hasattr(other, 'units') and self.units != other.units:
            return False
        return np.ndarray.__eq__(self, other).all()

    def __ne__(self, other):
        """ True if and only if the tensors have the same shape and are element-wise equal.
        """
        if not hasattr(other, "shape"):
            # TODO Decide what to do with the 1x1x1x... case when compared with a number
            return True
        elif hasattr(other, 'units') and self.units != other.units:
            return True
        if self.shape == other.shape:
            return np.not_equal(self, other).any()
        else:
            return True

    #-----------------------#
    # Arithmatic operations #
    #-----------------------#

    def __mul__(self, other):
        # TODO Move this to unit and CompositeUnit __rmul__
        if isunit(other):
            ret_val = copy(self)
            if self.units:
                ret_val.units = self.units * other
            else:
                ret_val.units = other
            return ret_val
        else:
            return super(Tensor, self).__mul__(other)

    def __mul_units__(self, other):
        if self.units is None:
            if other.units is None:
                return None
            else:
                return other.units
        elif other.units is None:
            return None
        else:
            return self.units * other.units

    def __div__(self, other):
        if isunit(other):
            ret_val = copy(self)
            if self.units:
                ret_val.units = self.units / other
            else:
                ret_val.units = other
            return ret_val
        else:
            if hasattr(other, 'units'):
                ret_val = super(Tensor, self).__div__(other)
                ret_val.units = self.__div_units__(other)
                return ret_val
            else:
                return super(Tensor, self).__div__(other)

    def __truediv__(self, other):
        if isunit(other):
            ret_val = copy(self)
            if self.units:
                ret_val.units = self.units / other
            else:
                ret_val.units = other
            return ret_val
        else:
            if hasattr(other, 'units'):
                ret_val = super(Tensor, self).__truediv__(other)
                ret_val.units = self.__div_units__(other)
                return ret_val
            else:
                return super(Tensor, self).__truediv__(other)

    def __div_units__(self, other):
        if self.units is None:
            if other.units is None:
                return None
            else:
                return other.units**-1
        elif other.units is None:
            return None
        else:
            return self.units / other.units

    def __add__(self, other):
        if self.units:
            if hasattr(other, 'units'):
                return super(Tensor, self).__add__(other.in_units(self.units))
            else:
                # Some other array_like...assume the units work out for now...
                return super(Tensor, self).__add__(other)
        else:
            return super(Tensor, self).__add__(other)

    def __sub__(self, other):
        if self.units:
            if hasattr(other, 'units'):
                return super(Tensor, self).__sub__(other.in_units(self.units))
            else:
                # Some other array_like...assume the units work out for now...
                return super(Tensor, self).__sub__(other)
        else:
            return super(Tensor, self).__sub__(other)

    #------------------------#
    # Output Representations #
    #------------------------#

    def __str__(self):
        if self.units:
            return super(Tensor, self).__str__() + " " + str(self.units)
        else:
            return super(Tensor, self).__str__()

    def __repr__(self):
        try:
            if isinstance(self, Matrix):
                # Fix the weird spacing problem by the Matrix class
                ret_val = repr(self.view(np.matrix)).title()
            else:
                ret_val = super(Tensor, self).__repr__()
            if self.units:
                return ret_val[:-1] + ', units=' + repr(self.units) + ')'
            else:
                return ret_val
        except Exception as e:
            # be sure and always print something...
            try:
                return '<{}, {}; normal __repr__ raised {} at line {} of {}>'.format(
                    type(self).__name__,
                    'x'.join(str(dim) for dim in self.shape),
                    traceback.extract_tb(sys.exc_info()[2])[-1][0],
                    sys.exc_info()[0].__name__
                )
            except Exception as e2:
                return ('<' + self.__class__.__name__ + " with serious problems (printing of exception"
                                                        " {} raised further exception {}>".format(e, e2))


    ##################
    # Static Methods #
    ##################

    @staticmethod
    def remove_phase_factor(*args, **kwargs):
        """ Given any number of `Tensor` objects that are the same up to a phase factor,
        return a list of `Tensor` objects that have been.  Copies are only made if
        necessary.  The chosen phase is the one that makes the first non-zero element
        (in the `numpy.nditer(tensor)` iterator, using the `cutoff` keyword argument
        or Tensor.zero_cutoff if one is not given) positive.
        """
        cutoff = kwargs.pop('cutoff', Tensor.same_tensor_cutoff)
        force_list = kwargs.pop('force_list', False)
        args = listify_args(*args, ignore=np.ndarray)
        if len(args) == 0:
            return []
        if sanity_checking_enabled:
            if not args[0].same_upto_phase_factor(*args[1:], **{'cutoff': cutoff}):
                raise ValueError("can't align phases of tensors that are not the same to a phase factor")
        it = np.nditer(args[0], flags=['multi_index'])
        idx = None
        for i in it:
            if abs(i) < Tensor.zero_cutoff:
                idx = it.multi_index
                break
        if idx is None:
            # it's the zero vector (and so are all of the rest, assuming sanity
            #   checking is enabled and/or the user is being safe)
            if len(args) == 1 and not force_list:
                return args[0]
            else:
                return args
        ret_val = []
        for arg in args:
            if arg[idx] < 0:
                ret_val.append(-1.0 * arg)
            else:
                ret_val.append(arg)
        if len(ret_val) == 1 and not force_list:
            return ret_val[0]
        else:
            return ret_val

    ###########
    # Methods #
    ###########

    def max_abs(self):
        if all(s == 1 for s in self.shape):
            return abs(self).ravel()[0]
        return max(*abs(self).ravel())

    #-----------------#
    # Inquiry methods #
    #-----------------#

    def is_zero(self, cutoff=None):
        """ Returns True if all elements of `self` have absolute values less than `cutoff`, which defaults
        to `Tensor.zero_cutoff`
        .. note::
           `Tensor.zero_cutoff` is treated as a pseudo-class attribute for Tensor instances, meaning individual
            instances can also set a `zero_cutoff` attribute which will take precidence over the class-level
            default.

        """
        if cutoff:
            cutoff = strip_units(cutoff, self.units)
        else:
            cutoff = self.zero_cutoff
        return super(Tensor, self).is_zero(cutoff)

    def is_square(self):
        if len(self.shape) <= 1:
            return False
        return all(s == self.shape[0] for s in self.shape[1:])

    def is_symmetric(self, cutoff=1e-10):
        if not self.is_square():
            return False
        return all(np.all(abs(self - self.transpose(axes)) < cutoff) for axes in list(permutations(range(len(self.shape))))[1:])

    def is_antisymmetric(self, cutoff=1e-10):
        if not self.is_square():
            return False
        return all(
            np.all(abs(self - self.transpose(axes)) < cutoff) if permutation_is_even(axes) else
            np.all(abs(self + self.transpose(axes)) < cutoff)
                for axes in list(permutations(range(len(self.shape))))[1:])

    #-----------------------------#
    # Methods that return Tensors #
    #-----------------------------#

    def copy_shape(self, **kwargs):
        """
        Copy all properties of the tensor except for the data.
        """
        kw = self.__copy_kwargs__()
        kw.update(kwargs)
        return self.__class__(**kw)

    def flatten(self, order = 'C'):
        """ Same as `numpy.ndarray.flatten`, but modified to return a `Vector` object

        :Examples:


        >>> t = Tensor([[[1,2],[3,4]],[[5,6],[7,8]]])
        >>> t
        Tensor([[[ 1.,  2.],
                [ 3.,  4.]],
        <BLANKLINE>
               [[ 5.,  6.],
                [ 7.,  8.]]])
        >>> t.flatten()
        Vector([ 1.,  2.,  3., ...,  6.,  7.,  8.])

        """
        return super(Tensor, self).flatten(order).view(Vector)

    def in_units(self, other_units):
        if self.units is None and other_units is not None:
            raise IncompatibleUnitsError(self.units, other_units)
        else:
            ret_val = copy(self)
            ret_val.units = other_units
            conv = self.units.to(other_units)
            ret_val = ret_val * conv
            return ret_val

    def linearly_transformed(self, transmat, backwards=False):
        if backwards:
            einsum_args = [self, range(len(self.shape))]
            einsum_args += sum(([transmat, [i, i+len(self.shape)]] for i in xrange(len(self.shape))), [])
            einsum_args.append([i+len(self.shape) for i in xrange(len(self.shape))])
            return np.einsum(*einsum_args).view(self.__class__)
        else:
            einsum_args = [self, range(len(self.shape))]
            einsum_args += sum(([transmat, [i+len(self.shape), i]] for i in xrange(len(self.shape))), [])
            einsum_args.append([i+len(self.shape) for i in xrange(len(self.shape))])
            return np.einsum(*einsum_args).view(self.__class__)

    def reindexed(self, new_indices, dimensions=None, reverse=False):
        ret_val = copy(self)
        old_ret_val = self
        if dimensions is None:
            dimensions = range(len(self.shape))
        for dim in dimensions:
            idxs = [slice(None)] * len(self.shape)
            new_idxs = [slice(None)] * len(self.shape)
            for old_idx, new_idx in enumerate(new_indices):
                if reverse:
                    idxs[dim] = new_idx
                    new_idxs[dim] = old_idx
                else:
                    idxs[dim] = old_idx
                    new_idxs[dim] = new_idx
                ret_val[tuple(new_idxs)] = old_ret_val[tuple(idxs)]
            old_ret_val = copy(ret_val)
        return ret_val
        # Fail safe:
        #ret_val = copy(self)
        #if dimensions is None:
        #    for idxs in product(*map(xrange, self.shape)):
        #        new_idxs = tuple(new_indices[i] for i in idxs)
        #        if reverse:
        #            ret_val[idxs] = self[new_idxs]
        #        else:
        #            ret_val[new_idxs] = self[idxs]
        #else:
        #    for idxs in product(*map(xrange, self.shape)):
        #        new_idxs = tuple(new_indices[i] if d in dimensions else i for d, i in enumerate(idxs))
        #        if reverse:
        #            ret_val[idxs] = self[new_idxs]
        #        else:
        #            ret_val[new_idxs] = self[idxs]
        #return ret_val

    #-----------#
    # Iterators #
    #-----------#

    def iter_vectors(self, with_indices=False):
        it = super(Tensor, self).iter_vectors(True)
        for val, idxs in it:
            yieldval = Vector(val, units=self.units)
            if with_indices:
                yield yieldval, idxs
            else:
                yield yieldval


    #----------------#
    # Output methods #
    #----------------#

    def formatted_string(self, **kwargs):
        if 'name' not in kwargs and self.name is not None:
            kwargs.update(name=self.name)
        if len(self.shape) == 1:
            raise NotImplementedError
        elif len(self.shape) == 2:
            fmt = MatrixFormatter(**kwargs)
        else:
            fmt =TensorFormatter(**kwargs)
        return fmt.format(self)

    # TODO finish writing this for large tensors and tensors of dimension other than 2
    def zero_structure(self,
            max_width=120,
            row_label_width=5,
            one_based=True,
            positive_char='+',
            zero_char='0',
            one_char='1',
            negative_char='-',
            cutoff=1e-10):
        ret_val = ''
        off = 1 if one_based else 0
        if len(self.shape) <= 1:
            raise NotImplementedError
        for topidxs in product(*map(xrange, self.shape[:-2])):
            if len(topidxs) > 0:
                ret_val += "Zero structure for {}:\n".format(",".join(str(i+off) for i in topidxs))
            top_row_start = " " * (row_label_width + 1)
            top_labels = ''
            for i in xrange(self.shape[-1]):
                if (i + off) % 10 == 0:
                    top_labels += "|"
                else:
                    top_labels += str((i+off)%10)
            wdth = max_width - row_label_width - 1
            if len(top_labels) < wdth:
                ret_val += top_row_start + top_labels
            else:
                raise NotImplementedError
            for row in xrange(self.shape[-2]):
                ret_val += "\n" + ("{:>"+str(row_label_width)+"d} ").format(row+1)
                for col in xrange(self.shape[-1]):
                    idxs = topidxs + (row, col)
                    val = self[idxs]
                    if abs(1.0 - val) < cutoff:
                        ret_val += one_char
                    elif val < -cutoff:
                        ret_val += negative_char
                    elif val > cutoff:
                        ret_val += positive_char
                    else:
                        ret_val += zero_char
            if len(topidxs) > 0 and tuple(i+1 for i in topidxs) != self.shape[:-2]:
                ret_val += "\n\n"
        return ret_val




# Using this class is a very bad idea.  It makes development a nightmare.  Don't do it.
class ComputableTensor(Tensor):
    """ A tensor with elements that can be computed as needed.
    .. warning:
        Using this class is a very bad idea.  It makes development a nightmare.  Don't do it.
    """

    sentinal_value = None
    compute_function = None

    @typechecked(
        compute_function=(callable, None),
        uncomputed=(bool, None))
    def __tensor_init__(self,
            compute_function=None,
            uncomputed=False,
            sentinal_value=None,
            **kwargs):
        self.compute_function = compute_function
        self.sentinal_value = sentinal_value or float('nan')
        if uncomputed:
            # set all of the values in self to nan until they are computed
            #   do this by skipping all other __setitem__ methods up to that
            #   of numpy.ndarray.
            super(Tensor, self).__setitem__(Ellipsis, sentinal_value)
        if compute_function is None and uncomputed:
            raise ValueError("need compute_function to initialize an uncomputed tensor")
        super(ComputableTensor, self).__tensor_init__(**kwargs)
        self.__init_called = True

    def is_sentinal_value(self, test_val):
        if math.isnan(self.sentinal_value):
            return math.isnan(test_val)
        elif math.isinf(self.sentinal_value) and self.sentinal_value < 0:
            if self.sentinal_value < 0:
                # negative infinity
                return math.isinf(test_val) and test_val < 0
            else:
                # positive infinity
                return math.isinf(test_val) and test_val > 0
        else:
            return test_val == self.sentinal_value

    ###################
    # Special Methods #
    ###################

    def __array_wrap__(self, out_arr, context=None):
        if self.compute_function is not None:
            raise NotImplementedError(
                "some numpy ufuncs on ComputableTensor objects cause"
                " problems.  View-cast to a Tensor first."
            )
        else:
            return super(ComputableTensor, self).__array_wrap__(out_arr, context)

    #---------------------#
    # Container Emulation #
    #---------------------#

    def __getitem__(self, item):
        ret_val = super(ComputableTensor, self).__getitem__(item)
        if self.compute_function:
            if isinstance(ret_val, EinsumTensor):
                # TODO remove this code.  It shouldn't ever get called, since einsum tensor now always gets the appropriate slice before contracting
                ## make sure all of the values are filled in...
                #slices = [i.slice for i in ret_val.indices]
                #it = np.nditer(ret_val.tensor, flags=['multi_index'])
                #for val in it:
                #    if self.is_sentinal_value(val):
                #        offset_idxs = tuple(idx + s.start for idx, s in zip(it.multi_index, slices))
                #        val = self.compute_function(self, offset_idxs)
                #        if math.isnan(val):
                #            raise ValueError('compute function returned not-a-number for indices ({})'.format(
                #                ', '.join(str(i) for i in offset_idxs)
                #            ))
                #        self[offset_idxs] = val
                #return ret_val
                raise NotImplementedError
            else:
                # Regular case, not an EinsumTensor
                item = (item,) if not isinstance(item, tuple) else item
                if all(isinstance(i, (type(Ellipsis), slice, int, long)) for i in item):
                    if np.isscalar(ret_val):
                        if self.is_sentinal_value(ret_val):
                            if sanity_checking_enabled and not all(isinstance(i, int) for i in item):
                                raise NotImplementedError("{} can't get item '{}'".format(
                                    self.__class__.__name__,
                                    item
                                ))
                            ret_val = self.compute_function(self, item)
                            if self.is_sentinal_value(ret_val):
                                raise ValueError('compute function returned the sentinal value "{}"'
                                                 ' for indices ({})'.format(
                                    self.sentinal_value,
                                    ', '.join(str(i) for i in item)
                                ))
                            self[item] = ret_val
                        return ret_val
                    else:
                        # Fill in the needed parts.
                        def get_iter(part, idx):
                            if isinstance(part, (int, long)):
                                #noinspection PyRedundantParentheses
                                return (part,)
                            elif isinstance(part, slice):
                                return islice(xrange(self.shape[idx]), part.start, part.stop, part.step)
                            else: # part is Ellipsis
                                return xrange(self.shape[idx])
                        # iterate over the values needed, checking to see
                        needed = item + (Ellipsis,) * (len(self.shape) - len(item))
                        for idxs in product(*[get_iter(part, idx) for idx, part in enumerate(needed)]):
                            # Note the skip of superclasses!
                            #noinspection PyCallByClass,PyTypeChecker
                            val = np.ndarray.__getitem__(self, idxs)
                            if self.is_sentinal_value(val):
                                computed = self.compute_function(self, idxs)
                                if self.is_sentinal_value(computed):
                                    raise ValueError('compute function returned not-a-number'
                                                     ' for indices ({})'.format(
                                        ', '.join(str(i) for i in item)))
                                #noinspection PyCallByClass,PyTypeChecker
                                np.ndarray.__setitem__(self, idxs, computed)
                        return super(ComputableTensor, self).__getitem__(item)
                else:
                    raise NotImplementedError("{} can't get item '{}'".format(
                        self.__class__.__name__,
                        item
                    ))
        else:
            return ret_val

    def __setitem__(self, key, value):
        if self.compute_function is not None:
            # Set items individually to avoid infinite recursion.
            #   If a large chunk of a tensor is to be set, numpy.ndarray first
            #   calls __getitem__ with the chunk's indices.  But in the case of
            #   a ComputedTensor, the __getitem__ function itself may need to set the item
            #   (or the compute function defined by the user may itself set the value in chunks
            #   to avoid being called too many times, since some values like cross products
            #   are far more efficient to compute in chunks).
            #----------------------------------------#
            # if a single item is given, handle it as a length-1 tuple
            key = (key,) if not isinstance(key, tuple) else key
            # for now, only worry about cases where all the indices are integers
            if all(isinstance(k, int) for k in key):
                keydims = len(key)
                mydims = len(self.shape)
                if keydims < mydims:
                    # Make the shapes match by broadcasting to an array of the right shape
                    valary = np.zeros(self.shape[keydims:])
                    valary[...] = value
                    # construct the iterater over the broadcasted array
                    valit = np.nditer(valary, flags=['multi_index'])
                    # now assign the elements individually
                    for val in valit:
                        validxs = key + valit.multi_index
                        super(ComputableTensor, self).__setitem__(validxs, val)
                else:
                    # otherwise, we're not setting a chunk, so defer to super
                    super(ComputableTensor, self).__setitem__(key, value)
            elif any(isinstance(i, slice) for i in key):
                # TODO handle slice setting and other such complicated nonsense
                raise NotImplementedError
            else:
                # we're up to something else weird like Einstein summation.  Appeal to super!
                super(ComputableTensor, self).__setitem__(key, value)
        else:
            super(ComputableTensor, self).__setitem__(key, value)


    def __contains__(self, item):
        """Returns True if and only if the tensor has the item without having to compute it"""
        # Deal with the special case of Coordinate objects being passed in...
        if self.compute_function is None:
            # TODO Change/get rid of this.  It should act like the normal numpy.ndarray contains
            return True
        else:
            val = super(ComputableTensor, self).__getitem__(item)
            if np.isscalar(val):
                return not self.is_sentinal_value(val)
            else:
                return all(not self.is_sentinal_value(v) for v in val.flat)

    ###########
    # Methods #
    ###########

    def fill(self):
        if self.compute_function is not None:
            for idxs in product(*map(xrange, self.shape)):
                self.__getitem__(idxs)
            # now make it back into a regular tensor
            self.compute_function = None

    #-----------------#
    # Inquiry methods #
    #-----------------#

    def hasitem(self, item):
        """Returns True if and only if the tensor has the item without having to compute it.  Same as `item in self`."""
        return item in self


#####################
# Dependent Imports #
#####################

from grendel.gmath.einsum import EinsumTensor, EinsumContraction, EinsumSum
from grendel.gmath.einsum_indices import IndexRange
from grendel.gmath.vector import Vector, LightVector
from grendel.gmath.matrix import Matrix
from grendel.output.tensor_printer import MatrixFormatter, TensorFormatter

