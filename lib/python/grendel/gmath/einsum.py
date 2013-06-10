# -*- coding: utf-8 -*-
""" Module containing classes to make Einstein summation straitforward and easy.
These shouldn't need to be accessed by the outside world directly.
"""
from __future__ import division, print_function
from copy import copy, deepcopy
import sys
from grendel import sanity_checking_enabled

import numpy as np
import re
from numbers import Number
from collections import Iterable

from grendel.util.decorators import forwardable, def_delegators

__all__ = []

@forwardable
class EinsumTensor(object):
    """ A Tensor with indices, for performing contractions.
    ..note::
      Don't store EinsumTensors!  The should go away by the end of the line, either by contraction
      or by summation into the destination.  Basically, things should work like you would expect if
      you don't do this.
    """

    ##############
    # Attributes #
    ##############

    indices = None
    sub_indices = None
    _tensor = None
    coeff = None

    ##################
    # Initialization #
    ##################

    def __new__(cls,
            indices,
            tensor=None,
            shape=None,
            coeff=None,
            subidxs=None,
            known_indices=False,
            index_range_set=None):
        """ Initialize a EinsumTensor with `indices` and a reference to `tensor`, or create a new
        Tensor to refer to with the give shape
        """
        self = object.__new__(cls)
        if coeff is None:
            self.coeff = 1.0
        else:
            self.coeff = coeff
        #----------------------------------------#
        # Indices determine slices, sub_indices determine einstein summation (if sub_indices
        #   are not used, then the two tuples are identical)
        if subidxs is None:
            self.indices, self.sub_indices = EinsumTensor.split_indices(indices, include_sub=True)
        else:
            # Mostly for calling from __copy__
            self.indices, self.sub_indices = indices, subidxs
        #----------------------------------------#
        if tensor is not None:
            self._tensor = tensor
        elif known_indices:
            idx_rng_set = index_range_set if index_range_set is not None else IndexRange.global_index_range_set
            self._tensor = Tensor(indices=self.sub_indices, index_range_set=idx_rng_set)
        elif shape is not None:
            self._tensor = Tensor(shape=shape)
        else:
            raise TypeError("EinsumTensor initialization must have either a tensor or a shape to proceed.")
        #----------------------------------------#
        # See if it's fully internally contracted.  If so, perform the internal contraction
        #   and return a float
        if all(self.sub_indices.count(i) == 2 for i in self.sub_indices):
            idxs = list(set(self.sub_indices))
            idxmap = {}
            for n, i in enumerate(idxs):
                idxmap[i] = n
            return np.einsum(
                self.sliced_tensor,
                [idxmap[i] for i in self.sub_indices],
                []
            )
        #----------------------------------------#
        if self.indices != self.sub_indices and not self._has_only_known_indices():
            raise ValueError('sub-indices can only be used with main indices from declared ranges')
        #----------------------------------------#
        return self


    ##############
    # Delegators #
    ##############

    def_delegators('sliced_tensor',
        'zero',
        'diagonal',
        'formatted_string',
        'zero_structure',
        'shape',
    )

    ##############
    # Properties #
    ##############

    @property
    def sliced_tensor(self):
        return self._tensor[self._get_slices()].view(Tensor)

    @property
    def t(self):
        return self.sliced_tensor

    @property
    def _(self):
        return self.sliced_tensor

    @property
    def name(self):
        if self._tensor.name is not None:
            return self._tensor.name
        else:
            return "(unnamed tensor)"

    ###################
    # Special Methods #
    ###################

    def __copy__(self):
        return EinsumTensor(indices=self.indices, tensor=self._tensor, coeff=self.coeff, subidxs=self.sub_indices)

    def __deepcopy__(self, memo):
        return EinsumTensor(indices=self.indices, tensor=deepcopy(self._tensor, memo), coeff=self.coeff, subidxs=self.sub_indices)

    #----------------------#
    # Arithmetic Operators #
    #----------------------#

    def __neg__(self):
        self.coeff *= -1.0
        return self

    def __str__(self):
        if self.indices == self.sub_indices:
            return self.name + '["' + ",".join(self.sub_indices) + '"]'
        else:
            # TODO str for EinsumTensors with subindices
            return repr(self) # for now

    def __repr__(self):
        return "<EinsumTensor indices = ['{0}'],\n"\
               "    tensor = {1}\n>".format("', '".join(self.indices), str(self._tensor))

    def __mul__(self, other):
        if isinstance(other, EinsumTensor):
            return EinsumContraction(self, other)
        elif isinstance(other, EinsumContraction):
            result = other.append_tensor(self)
            if result is not None:
                return result
            else:
                return other
        elif isinstance(other, Number):
            rv = copy(self)
            rv.coeff *= other
            return rv
        else: # pragma: no cover
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Number):
            return self.__mul__(other)
        else: # pragma: no cover
            return NotImplemented

    def __add__(self, other):
        if isinstance(other, (EinsumTensor, EinsumContraction)):
            return EinsumSum(self, other)
        elif isinstance(other, EinsumSum):
            other.summands.append(self)
            return other
        else: # pragma: no cover
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, (EinsumTensor, EinsumContraction, EinsumSum)):
            return EinsumSum(self, -other)
        else: # pragma: no cover
            return NotImplemented

    def __iadd__(self, other):
        if isinstance(other, (EinsumTensor, EinsumContraction)):
            return EinsumSum(other).sum_into(self, accumulate=True)
        elif isinstance(other, EinsumSum):
            return other.sum_into(self, accumulate=True)
        else: # pragma: no cover
            return NotImplemented

    def __isub__(self, other):
        if isinstance(other, (EinsumTensor, EinsumContraction)):
            return EinsumSum(-other).sum_into(self, accumulate=True)
        elif isinstance(other, EinsumSum):
            return (-other).sum_into(self, accumulate=True)
        else: # pragma: no cover
            return NotImplemented

    def __imul__(self, other):
        if isinstance(other, Number):
            t = self.sliced_tensor
            #noinspection PyUnusedLocal
            t *= other
        elif isinstance(other, (EinsumTensor, EinsumContraction, EinsumSum)):
            raise NotImplementedError
        else: # pragma: no cover
            return NotImplemented

    def __idiv__(self, other):
        if isinstance(other, Number):
            #noinspection PyTypeChecker
            self.__imul__(1.0 / other)
        elif isinstance(other, EinsumTensor):
            #noinspection PyPropertyAccess
            tmp = EinsumSum(other)
            self.__idiv__(tmp)
        elif isinstance(other, EinsumContraction):
            dest = EinsumTensor(self.sub_indices, shape=self.sliced_tensor.shape)
            other.contract(dest)
            t = self.sliced_tensor
            #noinspection PyUnusedLocal
            t /= dest.sliced_tensor
        elif isinstance(other, EinsumSum):
            dest = EinsumTensor(self.sub_indices, shape=self.sliced_tensor.shape)
            other.sum_into(dest)
            t = self.sliced_tensor
            #noinspection PyUnusedLocal
            t /= dest.sliced_tensor
        else: # pragma: no cover
            return NotImplemented
    __itruediv__ = __idiv__

    #---------------------#
    # Container Emulation #
    #---------------------#

    def __contains__(self, item):
        return item in self.sliced_tensor


    #-----------------------------#
    # Numpy "array-like" behavior #
    #-----------------------------#

    def __array__(self):
        """ Included for compatibility with numpy
        """
        return self.sliced_tensor

    #----------------------#
    # Comparison Operators #
    #----------------------#

    def __eq__(self, other):
        if isinstance(other, EinsumTensor):
            tmp = EinsumTensor(indices=self.indices, shape=self.shape, coeff=self.coeff, subidxs=self.sub_indices)
            other.sort_to(tmp, multiplier=other.coeff)
            return self.sliced_tensor == tmp.sliced_tensor
        elif not isinstance(other, (EinsumContraction, EinsumSum)):
            return self.sliced_tensor == other
        else:
            return NotImplemented

    ##################
    # Static Methods #
    ##################

    @staticmethod
    def split_indices(in_indices, include_sub=False):
        """ Split a string or list/tuple into a tuple of indices

        Examples
        --------
        >>> EinsumTensor.split_indices("i,j,k,l")
        ('i', 'j', 'k', 'l')
        >>> EinsumTensor.split_indices(("i,j","k,l"))
        ('i', 'j', 'k', 'l')

        """
        indices = []
        subindices = []
        def add_index(idx):
            if include_sub and '_' in idx:
                parts = idx.split('_')
                if len(parts) == 2 and len(parts[0]) > 0 and len(parts[1]) > 0:
                    indices.append(parts[0])
                    subindices.append(parts[1])
                else:
                    indices.append(idx)
                    subindices.append(idx)
            else:
                indices.append(idx)
                subindices.append(idx)
        if isinstance(in_indices, basestring):
            for idx in re.split(r'\s*,\s*', in_indices):
                add_index(idx)
        elif isinstance(in_indices, Iterable):
            for subpart in in_indices:
                for idx in re.split(r'\s*,\s*', subpart):
                    add_index(idx)
        if include_sub:
            return tuple(indices), tuple(subindices)
        else:
            return tuple(indices)


    ###########
    # Methods #
    ###########

    def sort_to(self, dest, multiplier = 1.0):
        idxmap = {}
        if sanity_checking_enabled and set(self.sub_indices) != set(dest.sub_indices):
            raise ValueError("can't sort tensor with indices '{}' into tensor with indices '{}'".format(
                self.sub_indices, dest.sub_indices
            ))
        for n, i in enumerate(set(self.sub_indices + dest.sub_indices)):
            idxmap[i] = n
        return np.einsum(
            multiplier * self.sliced_tensor,
            [idxmap[i] for i in self.sub_indices],
            [idxmap[i] for i in dest.sub_indices],
            out=dest.sliced_tensor)

    ###################
    # Private Methods #
    ###################

    def _get_slices(self):
        if self._tensor.indices is None:
            # Not using index ranges
            return tuple([Ellipsis] * len(self._tensor.shape))
        else:
            # make sure all of the indices are well-known...
            if len(self.indices) != len(self._tensor.indices):
                raise IndexError('dimension mismatch ({} != {})'.format(self.indices, self._tensor.indices))
            if self._has_only_known_indices():
                slices = []
                for aidx, myidx in zip(self.indices, self._tensor.indices):
                    arange = self._tensor.index_range_set.known_ranges[aidx]
                    myrange = self._tensor.index_range_set.known_ranges[myidx]
                    if arange is myrange:
                        slices.append(Ellipsis)
                    elif is_subrange(arange, myrange):
                        slices.append(arange.slice_in(myrange))
                return tuple(slices)
            else:
                badidx = [self.indices[i] for i in xrange(len(self.indices))
                    if not is_subrange(
                        self._tensor.index_range_set.known_ranges[self.indices[i]],
                        self._tensor.index_range_set.known_ranges[self._tensor.indices[i]]
                    )
                ]
                raise IndexError("Index '{0[0]}' unknown index range used on a tensor '{1}' with indices from"
                                 " well-defined ranges.  Probably a typo?  If not, check your Tensor and"
                                 " IndexRange constructors.".format(
                    badidx,
                    self._tensor.name if self._tensor.name is not None else "(unnamed tensor)"
                ))

    def _has_only_known_indices(self):
        if self._tensor.index_range_set is None:
            return False
        else:
            return all(
                is_subrange(
                    self._tensor.index_range_set.known_ranges[self.indices[i]],
                    self._tensor.index_range_set.known_ranges[self._tensor.indices[i]]
                ) for i in xrange(len(self.indices)))

##################################################################################
#                             EinsumContraction                                  #
##################################################################################

class EinsumContraction(object):
    """ A contraction between two EinsumTensors, an EinsumTensor and an EinsumContraction, or two EinsumContractions
    """

    ####################
    # Class Attributes #
    ####################

    factorize_contraction = True
    print_factorization = False

    ##############
    # Attributes #
    ##############

    tensors = None
    coeff = None

    ##################
    # Initialization #
    ##################

    def __new__(cls, left, right, coeff=None):
        ret_val = object.__new__(cls)
        if coeff is not None:
            # Only used when called from __copy__
            ret_val.coeff = coeff
            # ignore the coefficients from the left and right parts
            left.coeff = 1.0
            right.coeff = 1.0
        else:
            ret_val.coeff = 1.0
        ret_val.tensors = []
        ret_val.append_tensor(left)
        result = ret_val.append_tensor(right)
        if result is not None:
            # Coefficient multiplied in contract, don't need to do it here
            return result
        else:
            return ret_val


    ##############
    # Properties #
    ##############

    @property
    def external_indices(self):
        all_idxs = []
        for tens in self.tensors:
            all_idxs.extend(tens.sub_indices)
        return tuple(a for a in all_idxs if all_idxs.count(a) == 1)

    @property
    def internal_indices(self):
        all_idxs = []
        for tens in self.tensors:
            all_idxs.extend(tens.sub_indices)
        return tuple(sorted(list(set([a for a in all_idxs if all_idxs.count(a) == 2]))))

    ###################
    # Special Methods #
    ###################

    def __copy__(self):
        rv = EinsumContraction(self.tensors[0], self.tensors[1], coeff=self.coeff)
        for tensor in self.tensors[2:]:
            rv.tensors.append(tensor)
        return rv

    def __neg__(self):
        return self.__mul__(-1.0)

    def __add__(self, other):
        if isinstance(other, (EinsumContraction, EinsumTensor)):
            return EinsumSum(self, other)
        else: # pragma: no cover
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, (EinsumContraction, EinsumTensor)):
            return EinsumSum(self, -other)
        else: # pragma: no cover
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, EinsumTensor):
            result = self.append_tensor(other)
            if result is not None:
                return result
            else:
                return self
        elif isinstance(other, Number):
            rv = copy(self)
            rv.coeff *= other
            return rv
        else: # pragma: no cover
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Number):
            rv = copy(self)
            rv.coeff *= other
            return rv
        else: # pragma: no cover
            return NotImplemented

    def __imul__(self, other):
        # This should never get called from outside of this file...
        if isinstance(other, Number):
            self.coeff *= other
            return self
        else: # pragma: no cover
            return NotImplemented

    #------------------------#
    # Output Representations #
    #------------------------#

    def __str__(self):
        return " * ".join(str(t) for t in self.tensors)


    ###########
    # Methods #
    ###########

    def append_tensor(self, tens):
        self.tensors.append(tens)
        self.coeff *= tens.coeff
        if self.is_dot_product():
            return self.contract(None)
        else:
            return None


    def is_dot_product(self):
        return len(self.external_indices) == 0


    # TODO @hard optimal automatic factoring
    # TODO raise more readable errors
    def contract(self, dest=None):
        all_idxs = []
        for tens in self.tensors:
            all_idxs.extend(tens.sub_indices)
        #========================================#
        if sanity_checking_enabled:
            bad_idxs = [a for a in all_idxs if all_idxs.count(a) > 2]
            if len(bad_idxs) > 0:
                raise ValueError("too many appearances of index '{}' in contraction".format(bad_idxs[0]))
        #========================================#
        if dest is not None:
            # Copy so as to avoid overwriting data that might be on the right hand side.
            rv = dest.sliced_tensor.copy_shape()
            dot_product = False
        else:
            # for now, dest=None is only allowed in the case of a dot-product-like contraction
            rv = 0.0  # initialize rv so that the code inspector will be happy
            dot_product = True
            if sanity_checking_enabled:
                bad_idxs = [a for a in all_idxs if all_idxs.count(a) != 2]
                if len(bad_idxs) > 0:
                    raise ValueError("index '{}' does not appear exactly twice in dot-product-like"
                                     " contraction with indices {}.".format(bad_idxs[0], all_idxs))
        #========================================#
        # TODO automatic factorization of dot products
        if dot_product or not self.factorize_contraction or len(self.tensors) <= 2:
            indices = list(set(all_idxs))
            idxmap = {}
            for n, i in enumerate(indices):
                idxmap[i] = n
            esum_args = sum(
                ([t.sliced_tensor, [idxmap[i] for i in t.sub_indices]]
                    for t in self.tensors), [])
            if not dot_product:
                esum_args.append([idxmap[i] for i in dest.sub_indices])
            else:
                esum_args.append([])
            try:
                if not dot_product:
                    np.einsum(*esum_args, out=rv)
                else:
                    rv = np.einsum(*esum_args)
            except: # pragma: no cover
                # Debugging breakpoint line (since we can't break in the numpy source code)
                # Should never get here; any errors should be raised before this
                raise
            #----------------------------------------#
            if self.coeff != 1.0:
                #noinspection PyUnusedLocal
                rv *= self.coeff
            if not dot_product:
                dest.sliced_tensor[...] = rv
        #----------------------------------------#
        # Rough automatic factorization...
        else:
            p = lambda x: None  # Initialize p so the code inspector will be happy
            if self.print_factorization:
                file = self.print_factorization if self.print_factorization is not True else sys.stdout
                p = lambda x: print(x, file=file)
                p("Factorization of contraction {} <= {}".format(str(dest), str(self)))
            left = self.tensors[0]
            for i, right in enumerate(self.tensors[1:]):
                contr = EinsumContraction(left, right, coeff=1.0)
                if i < len(self.tensors) - 2:
                    out_idxs = contr.external_indices
                    if self.print_factorization:
                        p('   {{intermediate {}}}["{}"] = {} * {}'.format(
                            i+1,  ",".join(out_idxs),  str(left),  str(right)))
                    out_shape = []
                    for idx in out_idxs:
                        if idx in left.sub_indices:
                            out_shape.append(left.shape[left.sub_indices.index(idx)])
                        else:
                            out_shape.append(right.shape[right.sub_indices.index(idx)])
                    # The new left tensor will be the contraction intermediate
                    if dest._has_only_known_indices():
                        # TODO Test code for this case
                        left = EinsumTensor(out_idxs,
                            known_indices=True,
                            index_range_set=dest._tensor.index_range_set)
                    else:
                        left = EinsumTensor(out_idxs, shape=out_shape)
                    left._tensor.name = "{{intermediate {}}}".format(i+1)
                    contr.contract(dest=left)
                else:
                    if self.print_factorization:
                        p('    {} = {} * {} '.format(
                            str(dest), str(left), str(right)))
                    contr.contract(dest=dest)
            # Return the sliced tensor, just like we do in the sliced part
            rv = dest.sliced_tensor
            if self.coeff != 1.0:
                #noinspection PyUnusedLocal
                rv *= self.coeff
        #========================================#
        return rv


##################################################################################
#                                  EinsumSum                                     #
##################################################################################

class EinsumSum(object):
    """
    """

    #################
    # Inner Classes #
    #################

    class termwise_printing(object):

        defaults = {
            'printing' : False,
            'print_cummulative' : False,
            'result_name' : None,
            'handler' : sys.stdout,
            'names' : None,
            'formatting_keywords' : {}
        }
        kwargs = None

        def __init__(self, **kwargs):
            self.enabled = kwargs.pop('enabled', True)
            self.kwargs = kwargs

        def __enter__(self):
            for kw in self.kwargs:
                if hasattr(EinsumSum, "_termwise_" + kw):
                    setattr(EinsumSum, "_termwise_" + kw, self.kwargs[kw])
                else:
                    raise TypeError("Unknown keyword '{}'".format(kw))
            EinsumSum._termwise_printing = self.enabled

        def __exit__(self, type, value, traceback):
            EinsumSum._termwise_printing = False
            for kw, default in self.defaults.items():
                setattr(EinsumSum, kw, default)

    ##############
    # Attributes #
    ##############

    summands = None
    coeff = None

    ############################
    # Private Class Attributes #
    ############################

    _termwise_printing = False
    _termwise_handler = None
    _termwise_names = None
    _termwise_print_cummulative = None
    _termwise_result_name = None
    _termwise_formatting_keywords = {}

    ##################
    # Initialization #
    ##################

    def __init__(self, left, right=None):
        if isinstance(left, (list, tuple)):
            self.summands = list(left)
        else:
            self.summands = [left]
        if right is not None:
            self.summands.append(right)
        new_summands = []
        # Flatten any inner sums
        for term in self.summands:
            if isinstance(term, EinsumSum):
                new_summands.extend(term.summands)
            else:
                new_summands.append(term)
        self.summands = new_summands

    ###################
    # Special Methods #
    ###################

    def __neg__(self):
        new_terms = []
        for term in self.summands:
            new_terms.append(-term)
        return EinsumSum(new_terms)

    def __add__(self, other):
        if isinstance(other, (EinsumTensor, EinsumContraction)):
            self.summands.append(other)
            return self
        elif isinstance(other, EinsumSum):
            self.summands.extend(other.summands)
            return self
        else: # pragma: no cover
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, (EinsumTensor, EinsumContraction)):
            self.summands.append(-other)
            return self
        elif isinstance(other, EinsumSum):
            self.summands.extend([-o for o in other.summands])
            return self
        else: # pragma: no cover
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Number):
            for term in self.summands:
                term.coeff *= other
            return self
        else: # pragma: no cover
            return NotImplemented

    ###########
    # Methods #
    ###########

    def sum_into(self, dest, accumulate=False):
        # Copy so as to avoid overwriting data that might be on the right hand side.
        dst_slice = copy(dest.sliced_tensor)
        if not accumulate:
            dst_slice[...] = 0.0
        for tidx, term in enumerate(self.summands):
            # Contract if it is a contraction rather than a normal EinsumTensor
            if isinstance(term, EinsumContraction):
                dst_idxs = [idx for idx in dest.sub_indices if idx in term.external_indices]
                dst_shape = tuple(dst_slice.shape[i] for i, idx in enumerate(dest.sub_indices) if idx in term.external_indices)
                to_add = EinsumTensor(dst_idxs, shape=dst_shape)
                to_add._tensor.name = "(term {} of {})".format(tidx+1, dest.name)
                term.contract(dest=to_add)
                # Now continue in the normal fashion, treating the contracted tensor as a regular tensor
                term = to_add
            # Add the contribution from the given term
            if len(term.sub_indices) <= len(dest.sub_indices):
                # the term is not internally contracted
                # TODO: fix the following edge case
                #   This could fail if an internally contracted term is being broadcast back into the
                #   destaination, e.g. T[i,j,k,l] += k[j,c,c].  I can't think of an instance in which this
                #   would be used, but it should be at least recognized as a weakness
                if sanity_checking_enabled and not all(t in dest.sub_indices for t in term.sub_indices):
                    raise IndexError("index mismatch:  can't map '{}' to '{}'".format(
                        term.sub_indices, dest.sub_indices))
                out_axes = [np.newaxis] * len(dst_slice.shape)
                srt_idxs = sorted(term.sub_indices, key=dest.sub_indices.index)
                out_perm = list(term.sub_indices.index(i) for i in srt_idxs)
                srt_term = np.transpose(term.sliced_tensor, axes=out_perm)
                for idx in term.sub_indices:
                    out_axes[dest.sub_indices.index(idx)] = slice(None)
                to_add = term.coeff * srt_term.view(Tensor)[tuple(out_axes)]
                dst_slice += to_add
                #========================================#
                # printing for debugging purposes
                if EinsumSum._termwise_printing:
                    p = lambda x: print(x, file=EinsumSum._termwise_handler)
                    def tban(x):
                        p('+-' + ('-' * len(x)) + '-+')
                        p('| ' + x + ' |')
                        p('+-' + ('-' * len(x)) + '-+')
                    if EinsumSum._termwise_names is not None:
                        name = EinsumSum._termwise_names[tidx]
                    elif EinsumSum._termwise_result_name is not None:
                        name = "Term #{} of {}".format(tidx+1, EinsumSum._termwise_result_name)
                    elif dest._tensor.name is not None:
                        name = "Term #{} of {}".format(tidx+1, dest._tensor.name)
                    else:
                        name = "Term #{}".format(tidx+1)
                    tban(name)
                    if not isinstance(to_add, Tensor):
                        p(Tensor(to_add).formatted_string(**EinsumSum._termwise_formatting_keywords))
                    else:
                        p(to_add.formatted_string(**EinsumSum._termwise_formatting_keywords))
                    p("\n")
                    if EinsumSum._termwise_print_cummulative:
                        if EinsumSum._termwise_result_name is not None:
                            name = EinsumSum._termwise_result_name
                        elif dest._tensor.name is not None:
                            name = dest._tensor.name
                        else:
                            name = "(unlabeled tensor)"
                        tban("Sum of first {} term{} of {}".format(
                            tidx + 1,
                            '' if tidx == 0 else 's',
                            name
                        ))
                        p(dest.formatted_string(**EinsumSum._termwise_formatting_keywords))
            else:
                # Add the internal contraction to dest
                idxs = list(set(dest.sub_indices + term.sub_indices))
                idxmap = {}
                for n, i in enumerate(idxs):
                    idxmap[i] = n
                dst_slice += term.coeff * np.einsum(
                    term.sliced_tensor,
                    [idxmap[i] for i in term.sub_indices],
                    [idxmap[i] for i in dest.sub_indices])
        dest.sliced_tensor[...] = dst_slice
        return dest



#####################
# Dependent Imports #
#####################

from grendel.gmath.tensor import Tensor
from grendel.gmath.einsum_indices import is_subrange, IndexRange

