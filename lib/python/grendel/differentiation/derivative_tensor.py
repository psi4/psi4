from copy import copy
from itertools import permutations
from numbers import Number

import numpy as np

from grendel import type_checking_enabled, sanity_checking_enabled
from grendel.coordinates.coordinate import Coordinate
from grendel.gmath.tensor import Tensor, ComputableTensor
from grendel.output.tensor_printer import MatrixFormatter, TensorFormatter
from grendel.representations.representation import Representation
from grendel.util.abstract_bases import NonstringIterable
from grendel.util.decorators import typechecked
from grendel.util.strings import shortstr, indent
from grendel.util.units.value_with_units import hasunits

class RepresentationDependentTensor(ComputableTensor):
    """
    """

    ##############
    # Attributes #
    ##############

    permutational_symmetry = None
    first_dimension_different = None

    ######################
    # Private Attributes #
    ######################

    _representation = None
    __init_called = False
    _secondary_representation = None

    ##################
    # Initialization #
    ##################

    # TODO accept a PermutationGroup object (or something like that) for the permutational_symmetry attribute
    @typechecked(
        representation=Representation,
        permutational_symmetry=(bool, None),
        first_dimension_different=(bool, None),
        secondary_representation=(Representation, None))
    def __tensor_init__(self,
            representation,
            permutational_symmetry=False,
            first_dimension_different=False,
            secondary_representation=None,
            **kwargs):
        self.first_dimension_different=first_dimension_different
        self.representation = representation
        self.permutational_symmetry = permutational_symmetry
        if sanity_checking_enabled:
            if self.first_dimension_different and len(self.shape) == 1:
                raise ValueError("can't have the first dimension different for a"
                                 " one-dimensional tensor (i.e. for a Vector)")
        if secondary_representation:
            self.secondary_representation = secondary_representation
        super(RepresentationDependentTensor, self).__tensor_init__(**kwargs)
        self.__init_called = True


    ##############
    # Properties #
    ##############

    @property
    def value(self):
        return self.view(Tensor)

    @property
    def representation(self):
        return self._representation

    @representation.setter
    @typechecked(new_rep=Representation)
    def representation(self, new_rep):
        if sanity_checking_enabled:
            if self.first_dimension_different:
                if self.shape[0] != len(new_rep):
                    raise ValueError("dimensions of Tensor ({}) are not compatible"
                                     " with representation (length {})".format(
                        self.shape[0], len(new_rep)
                    ))
            else:
                if not all(s == len(new_rep) for s in self.shape):
                    raise ValueError("dimensions of Tensor ({}) are not compatible"
                                     " with representation (length {})".format(
                        self.shape[0], len(new_rep)
                    ))
        self._representation = new_rep
        # Freeze the new representation to make sure it doesn't change on us
        self._representation.freeze()

    @property
    def secondary_representation(self):
        return self._secondary_representation

    @secondary_representation.setter
    @typechecked(new_rep=Representation)
    def secondary_representation(self, new_rep):
        if sanity_checking_enabled:
            if not self.first_dimension_different:
                raise ValueError("can't set second_representation for a"
                                 " RepresentationDependentTensor that"
                                 " depends on only one representation "
                                 " for all dimensions; first_dimension_different"
                                 " must be set to associate two representations"
                                 " with one tensor")
            else:
                if any(s != len(new_rep) for s in self.shape[1:]):
                    raise ValueError("dimensions of Tensor:\n{}\nare not compatible"
                                     " with representation:\n{}".format(
                        self, new_rep
                    ))
        self._secondary_representation = new_rep
        # Freeze the new representation to make sure it doesn't change on us
        self._secondary_representation.freeze()

    ###################
    # Special Methods #
    ###################

    #-------------------------------#
    # Numpy ndarray "magic" methods #
    #-------------------------------#

    def __array_finalize__(self, obj):
        if not self.__init_called:
            # then we're probably calling from Tensor.__new__ (or just someone who's looking for trouble)
            # just forget about it...
            return
        # remember self is the *new* instance and obj is the old one!
        if obj is None:
            # Just call super and be done with it.
            super(RepresentationDependentTensor, self).__array_finalize__(obj)
            return
        #========================================#
        # If we're coming from something that doesn't have a representation attribute,
        # then just throw our hands up and give up
        new_rep = getattr(obj, 'representation', None)
        if new_rep is None:
            raise TypeError("can't create a RepresentationDependentTensor"
                            " by casting an object that doesn't have a"
                            " representation associated with it.")
        # Other checking will happen in the representation setter...
        self.representation = new_rep
        #========================================#
        # Now see if the other has a second representation...
        first_different = getattr(obj, 'first_dimension_different', None)
        if first_different:
            sec_rep = getattr(obj, 'secondary_representation', None)
            if sec_rep is None:
                raise TypeError("can't create a RepresentationDependentTensor"
                                " with a different first dimension by casting"
                                " an object that doesn't have a secondary"
                                " representation associated with it.")
            # Other checking will happen in the representation setter...
            self.first_dimension_different = True
            self.secondary_representation = new_rep
        else:
            self.first_dimension_different = False
        #========================================#
        # Check for permutational symmetry...
        psym = getattr(obj, 'permutational_symmetry', None)
        if psym:
            self.permutational_symmetry = psym
        #========================================#
        # Finally, call super
        super(RepresentationDependentTensor, self).__array_finalize__(obj)


    __passthrough_ufuncs__ = [np.equal, np.isnan, np.isinf]

    def __array_prepare__(self, out_arr, context=None):
        fail_error = NotImplementedError(
            "some numpy ufuncs ({} was called) on RepresentationDependentTensor objects cause"
            " problems.  View-cast to a Tensor first.".format(
                context[0] if context is not None else "unknown ufunc"
            )
        )
        if context is not None:
            #========================================#
            if context[0] is np.multiply:
                if isinstance(context[1][0], Number) or isinstance(context[1][1], Number):
                    # we're just scaling the array.  call super...
                    return super(RepresentationDependentTensor, self).__array_prepare__(out_arr, context)
            elif context[0] is np.divide:
                if isinstance(context[1][1], Number):
                    # we're just scaling the array.  call super...
                    return super(RepresentationDependentTensor, self).__array_prepare__(out_arr, context)
            #========================================#
            elif context[0] is np.add:
                other = None
                if context[1][0] is self:
                    other = context[1][1]
                else:
                    other = context[1][0]
                #----------------------------------------#
                if isinstance(other, RepresentationDependentTensor) and other.representation is self.representation:
                    # all is well.  call on up.  We'll take care of the transfer in the __array_wrap__
                    return super(RepresentationDependentTensor, self).__array_prepare__(out_arr, context)
                else:
                    # don't do any casting...
                    # TODO this doesn't do what I think it does...
                    return super(RepresentationDependentTensor, self).__array_prepare__(out_arr, context)
            #========================================#
            elif context[0] in self.__passthrough_ufuncs__:
                # just call super
                return super(RepresentationDependentTensor, self).__array_prepare__(out_arr, context)
            #========================================#
        raise fail_error

    def __array_wrap__(self, out_arr, context=None):
        fail_error = NotImplementedError(
            "some numpy ufuncs ({} was called) on RepresentationDependentTensor objects cause"
            " problems.  View-cast to a Tensor first.".format(
                context[0] if context is not None else "unknown ufunc"
            )
        )
        #========================================#
        def copy_over_stuff():
            kwargs = self.__copy_kwargs__(include_parent=False)
            out_arr.__tensor_init__(**kwargs)
            out_arr.units = self.units
        if context is not None:
            #========================================#
            if context[0] is np.multiply:
                if isinstance(context[1][0], Number) or isinstance(context[1][1], Number):
                    # we're just scaling the array.  transfer important stuff and call super
                    copy_over_stuff()
                    return super(RepresentationDependentTensor, self).__array_wrap__(out_arr, context)
            elif context[0] is np.divide:
                if isinstance(context[1][1], Number):
                    # we're just scaling the array.  transfer important stuff and call super
                    copy_over_stuff()
                    return super(RepresentationDependentTensor, self).__array_wrap__(out_arr, context)
            #----------------------------------------#
            elif context[0] is np.add:
                if context[1][0] is self:
                    other = context[1][1]
                else:
                    other = context[1][0]
                if isinstance(other, RepresentationDependentTensor) and other.representation is self.representation:
                    # all is well.  call on up.  We'll take care of the transfer in the __array_wrap__
                    return super(RepresentationDependentTensor, self).__array_wrap__(out_arr, context)
                else:
                    return super(RepresentationDependentTensor, self).__array_wrap__(out_arr.view(Tensor), context)
            #----------------------------------------#
            elif context[0] in self.__passthrough_ufuncs__:
                # just call super
                return super(RepresentationDependentTensor, self).__array_wrap__(out_arr, context)
            #========================================#
            else:
                raise fail_error
        else:
            # if context is none, ignore it (it's probably something like scalar multiplication)
            return super(RepresentationDependentTensor, self).__array_wrap__(out_arr, context)

    #---------------------#
    # Container Emulation #
    #---------------------#

    def __getitem__(self, item):
        item = (item,) if not isinstance(item, tuple) else item
        #----------------------------------------#
        # Deal with the special case of Coordinate objects being passed in...
        if sanity_checking_enabled:
            if any(isinstance(c, Coordinate) and c.is_orphaned() for c in item):
                raise ValueError("cannot set item using index of an orphaned coordinate")
        indices = tuple(i.index if isinstance(i, Coordinate) else i for i in item)
        #----------------------------------------#
        return super(RepresentationDependentTensor, self).__getitem__(indices)

    def __setitem__(self, item, value):
        item = (item,) if not isinstance(item, tuple) else item
        #----------------------------------------#
        # Deal with the special case of Coordinate objects being passed in...
        if sanity_checking_enabled:
            if any(isinstance(c, Coordinate) and c.is_orphaned() for c in item):
                raise ValueError("cannot set item using index of an orphaned coordinate")
        indices = tuple(i.index if isinstance(i, Coordinate) else i for i in item)
        #----------------------------------------#
        # Einstein summation and other hoopla will be handled in the superclass
        if self.permutational_symmetry and isinstance(indices, NonstringIterable) and all(isinstance(i, int) for i in indices):
            if sanity_checking_enabled and not np.isscalar(value):
                raise NotImplementedError("permutationally symmetric assignments of non-scalar objects is not yet implemented")
            if self.first_dimension_different:
                # Assign all the permutations of the remaining dimensions...
                for idxs in permutations(indices[1:]):
                    super(RepresentationDependentTensor, self).__setitem__((indices[0],) + idxs, value)
            else:
                # Permute all dimensions...
                for idxs in permutations(indices):
                    super(RepresentationDependentTensor, self).__setitem__(idxs, value)
        else:
            super(RepresentationDependentTensor, self).__setitem__(indices, value)

    #----------------------#
    # Copying and Pickling #
    #----------------------#

    def __copy_kwargs__(self, include_parent=True):
        if include_parent:
            ret_val = super(RepresentationDependentTensor, self).__copy_kwargs__()
        else:
            ret_val = {}
        ret_val.update(
            representation=self.representation,
            permutational_symmetry=self.permutational_symmetry,
            first_dimension_different=self.first_dimension_different,
            secondary_representation=self.secondary_representation
        )
        return ret_val


    #------------------------#
    # Output Representations #
    #------------------------#

    def __str__(self):
        namestrs = []
        valstrs = []
        it = np.nditer(self, flags=['multi_index'])
        for x in it:
            namestrs.append(", ".join(shortstr(self.representation[i]) for i in it.multi_index))
            valstrs.append("{: .8f}{}".format(float(x), (" " + str(self.units) if hasunits(self) else '')))
        width=max(len(s) for s in namestrs) + 4
        ret_val = self.__class__.__name__ + ":\n"
        for name, val in zip(namestrs, valstrs):
            ret_val += ("{:>"+str(width)+"}: {}\n").format(name, val)
        ret_val += "  for representation:\n{}".format(indent(str(self.representation)))
        return ret_val

    def __repr__(self):
        try:
            return self.__short_repr__()  # For now...
        except:
            # Always print something, no matter what...
            return super(RepresentationDependentTensor, self).__repr__()

    def __short_repr__(self):
        return "<{shape} {cls}>".format(
            shape='x'.join(str(s) for s in self.shape),
            cls=self.__class__.__name__,
        )


    ###########
    # Methods #
    ###########

    def formatted_string(self, **kwargs):
        if len(self.shape) == 1:
            raise NotImplementedError
        elif len(self.shape) == 2:
            fmt = MatrixFormatter(**kwargs)
        else:
            fmt =TensorFormatter(**kwargs)
        if self.first_dimension_different:
            first_lbls = [coord.name for coord in self.representation]
            other_lbls = [coord.name for coord in self.secondary_representation]
            lbls = [first_lbls] + [other_lbls] * (len(self.shape) - 1)
        else:
            lbls = [coord.name for coord in self.representation]
        return fmt.format(self, labels=lbls)

    def in_representation(self, new_rep):
        return self.representation.transform_tensor(self, new_rep)


class DerivativeTensor(RepresentationDependentTensor):

    ##############
    # Attributes #
    ##############

    compute_function = None
    collection = None

    ######################
    # Private Attributes #
    ######################

    _value_computed = None

    ##################
    # Initialization #
    ##################

    @typechecked(
        collection=('DerivativeCollection', None))
    def __tensor_init__(self,
            collection=None,
            **kwargs):
        self.collection = collection
        super(DerivativeTensor, self).__tensor_init__(**kwargs)
        self.__init_called = True

    ###################
    # Special Methods #
    ###################

    #---------------------#
    # Container Emulation #
    #---------------------#

    def __getitem__(self, item):
        item = (item,) if not isinstance(item, tuple) else item
        #----------------------------------------#
        # Coordinate objects being passed in is handled by superclass
        # Einstein summation and other hoopla will be handled in the superclass
        # View-cast to prevent a DerivativeTensor object from being returned,
        #   since that would be confusing (for instance, how does the first_dimension_different
        #   get propagated?).  This could be done, but it would be difficult...
        ret_val = super(DerivativeTensor, self).__getitem__(item)
        if isinstance(ret_val, Tensor) and len(ret_val.shape) > 0:
            # TODO make this work with Einstein summation (right now the EinsumTensor would have a tensor attribute of type DerivativeTensor, which is confusing for the same reasons.)
            ret_val = ret_val.view(Tensor)
        # TODO Cast 1D to vector and 2D to Matrix?
        return ret_val


    def __setitem__(self, item, value):
        # Just a call to superclass for now...
        # Deal with Coordinate and other stuff objects in the parent class
        super(DerivativeTensor, self).__setitem__(item, value)

    def __contains__(self, item):
        item = (item,) if not isinstance(item, tuple) else item
        #----------------------------------------#
        # Deal with the special case of Coordinate objects being passed in...
        if sanity_checking_enabled:
            if any(isinstance(c, Coordinate) and c.is_orphaned() for c in item):
                raise ValueError("cannot set item using index of an orphaned coordinate")
        indices = tuple(i.index if isinstance(i, Coordinate) else i for i in item)
        #----------------------------------------#
        return super(DerivativeTensor, self).__contains__(indices)

    #----------------------#
    # Copying and Pickling #
    #----------------------#

    def __copy_kwargs__(self):
        ret_val = super(DerivativeTensor, self).__copy_kwargs__()
        ret_val.update(
            compute_function=self.compute_function,
            # when copying, we don't want to obliterate the data
            uncomputed=False
        )


#####################
# Dependent Imports #
#####################

# needed for dynamic typechecking, do not delete
if type_checking_enabled:
    from grendel.differentiation.derivative_collection import DerivativeCollection
