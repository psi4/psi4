# -*- coding: utf-8 -*-
from grendel import sanity_checking_enabled, type_checking_enabled
from grendel.util.decorators import with_flexible_arguments
from grendel.util.metaclasses import SubscriptableClass
from grendel.util.overloading import listify_args
from grendel.util.strings import classname

__all__ = [
    "IndexRange",
    "DeclareIndexRange",
    "IndexRangeSet",
    "IndexingContext",
    "issubrange",
    "is_subrange"
]

def is_subrange(child_range, parent_range):
    """
    Returns True if `child_range` is `parent_range` or `child_range` is a child of `parent_range`
    (analogous to `isinstance` from the python standard library)
    """
    if child_range is parent_range:
        return True
    parents = []
    spot = child_range
    while spot.parent is not None:
        parents.append(spot.parent)
        spot = spot.parent
    return parent_range in parents
issubrange = is_subrange

def DeclareIndexRange(indices, begin_index_or_size, end_index=None, name=None, **kwargs):
    """ Alias for `IndexRange` constructor that more accurately describes what is actually going on.
    """
    return IndexRange(indices, begin_index_or_size, end_index, name, **kwargs)

class IndexRangeSet(object):
    """
    """

    ##############
    # Attributes #
    ##############

    known_ranges = None

    ##################
    # Initialization #
    ##################

    def __init__(self):
        self.known_ranges = {}

    ###################
    # Special Methods #
    ###################

    def __getitem__(self, item):
        return self.known_ranges[item]

    def __setitem__(self, key, item):
        self.known_ranges[key] = item
# Useful alias when used as a transparent context
IndexingContext = IndexRangeSet

class IndexRange(object):
    """ Allows for the definition of special indices that cover specific ranges and subranges of Tensor axes
    in einstein summations.

    Examples
    --------
    >>> IndexRange.clear_known_ranges()  # Don't do this in your program; this is just for the doctest
    >>> p = IndexRange('p,q,r,s', 5).with_subranges(
    ...   IndexRange('i,j,k,l', 0, 2),
    ...   IndexRange('a,b,c,d', 2, 5)
    ... )
    >>> p
    <IndexRange object covering slice 0:5 represented by indices ['p', 'q', 'r', 's']>
    >>> p.subranges[0]
    <IndexRange object covering slice 0:2 represented by indices ['i', 'j', 'k', 'l']>
    >>> p.subranges[1]
    <IndexRange object covering slice 2:5 represented by indices ['a', 'b', 'c', 'd']>
    >>> prime = IndexRange["m, m',m''", 10, "m range"].with_subranges(
    ...            IndexRange["n,n',n''", ..., 6, "n range"],
    ...            IndexRange["o,o',o''", 6, ...]
    ... )
    >>> prime
    <IndexRange named 'm range' represented by indices ['m', 'm'', 'm''']>
    >>> prime.subranges[0]
    <IndexRange named 'n range' represented by indices ['n', 'n'', 'n''']>
    >>> prime.subranges[1]
    <IndexRange object covering slice 6:10 represented by indices ['o', 'o'', 'o''']>


    """
    __metaclass__ = SubscriptableClass

    ####################
    # Class Attributes #
    ####################

    global_index_range_set = IndexRangeSet()

    ############################
    # Private Class Attributes #
    ############################

    _global_index_context_stack = []

    ##############
    # Attributes #
    ##############

    indices = None
    subranges = None
    index_range_set = None

    name = None
    """ A description of the range covered by the indices.  Useful for debugging. """

    ######################
    # Private Attributes #
    ######################

    _parent = None
    _begin_is_ellipsis = None
    _end_is_ellipsis = None
    _slice = None


    ##################
    # Initialization #
    ##################

    # TODO Error checking (e.g. begin <= end, etc)
    # TODO Ellipsis used to mean "to the end of the parent range" or "to the beginning of the parent range"
    # TODO use the SubscriptableClass __class_getitem__ to allow writing of constructors with (literal) ellipses and slices
    @with_flexible_arguments(
        optional=[
            ('index_range_set', 'set', 'in_set', 'idx_set')
        ]
    )
    def __init__(self, indices,
            begin_index_or_size_or_slices,
            end_index=None,
            name=None,
            parent=None,
            index_range_set=None):
        self.indices, self.sub_indices = EinsumTensor.split_indices(indices, include_sub=True)
        if self.indices != self.sub_indices:
            raise ValueError("can't include indices with sub-indices (such as '{}') in"
                             " IndexRange declaration.".format(
                [self.indices[i] + "_" + self.sub_indices[i] for i in xrange(len(self.indices)) if self.indices[i] != self.sub_indices[i]][0]
            ))
        self.subranges = []
        if parent:
            self.parent = parent
        if index_range_set is None:
            if self.parent is not None:
                self.index_range_set = self.parent.index_range_set
            else:
                self.index_range_set = IndexRange.global_index_range_set
        else:
            self.index_range_set = index_range_set
        self.name = name
        self._begin_is_ellipsis = False
        self._end_is_ellipsis = False
        self._parent = None
        if isinstance(begin_index_or_size_or_slices, int):
            if end_index is not None:
                if end_index is Ellipsis:
                    self._end_is_ellipsis = True
                    self.slice = slice(begin_index_or_size_or_slices, -1)
                else:
                    self.slice = slice(begin_index_or_size_or_slices, end_index)
            else:
                self.slice = slice(0, begin_index_or_size_or_slices)
        elif isinstance(begin_index_or_size_or_slices, slice):
            self.slice = begin_index_or_size_or_slices
        elif begin_index_or_size_or_slices is Ellipsis:
            self._begin_is_ellipsis = True
            if isinstance(end_index, int):
                self.slice = slice(-1, end_index)
            else:
                raise TypeError("Unsupport type for third argument to IndexRange constructor of type {0}".format(classname(end_index)))
        else:
            raise TypeError("Unsupport type for third argument to IndexRange constructor of type {0}".format(classname(begin_index_or_size_or_slices)))
        for idx in self.indices:
            if idx in self.index_range_set.known_ranges:
                raise IndexError(u"Index {0} is already part of index range {1}.".format(idx, self.index_range_set.known_ranges[idx]))
            self.index_range_set[idx] = self


    ##############
    # Properties #
    ##############

    @property
    def size(self):
        return self.end_index - self.begin_index

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, new_parent):
        if self._begin_is_ellipsis:
            self.begin_index = new_parent.begin_index
        if self._end_is_ellipsis:
            self.end_index = new_parent.end_index
        self._parent = new_parent

    @property
    def begin_index(self):
        return self._slice.start

    @begin_index.setter
    def begin_index(self, idx):
        if idx > self.end_index:
            raise IndexError("new begin_index {0} is after end_index {1}".format(idx, self.end_index))
        self.slice = slice(idx, self._slice.stop, self._slice.step)

    @property
    def end_index(self):
        return self._slice.stop

    @end_index.setter
    def end_index(self, idx):
        if sanity_checking_enabled:
            if idx < self.begin_index:
                raise IndexError("new end_index {0} is before begin_index {1}".format(idx, self.begin_index))
        self.slice = slice(self._slice.start, idx, self._slice.step)

    @property
    def slice(self):
        if sanity_checking_enabled:
            if self._end_is_ellipsis and self._slice.stop == -1:
                raise IndexError("Ellipsis used as end_index to denote remainder of parent range for "
                                 "IndexRange {0} that has no parent range.".format(self))
            if self._begin_is_ellipsis and self._slice.start == -1:
                raise IndexError("Ellipsis used as begin_index denote remainder of parent range for "
                                 "IndexRange {0} that has no parent range.".format(self))
        return self._slice

    @slice.setter
    def slice(self, new_slice):
        if type_checking_enabled:
            if not isinstance(new_slice, slice):
                raise TypeError("IndexRange.slice must an instance of 'slice'")
        if sanity_checking_enabled:
            if not new_slice.step is None and not new_slice.step == 1:
                raise NotImplementedError("IndexRange objects for slices with non-unitary steps are not yet implemented.")
            pass

        self._slice = new_slice


    #########################
    # Special Class Methods #
    #########################

    @classmethod
    def __class_getitem__(cls, args):
        """
        When given multiple arguments, this functions as a serogate constructor, allowing Ellipsis literals and
        slice literals to be given.

        """
        items = listify_args(args)
        if len(items) == 3 and isinstance(items[2], basestring):
            name = items.pop()
            return IndexRange(*items, name=name)
        else:
            return IndexRange(*items)

    ###################
    # Special Methods #
    ###################

    def __len__(self):
        return self.size

    def __str__(self):
        if self.name is None:
            return "<IndexRange object covering slice {0.begin_index}:{0.end_index}"\
                   " represented by indices ['{1}']>".format(self, "', '".join(self.indices))
        else:
            return "<IndexRange named '{0.name}'"\
                   " represented by indices ['{1}']>".format(self, "', '".join(self.indices))
    __repr__ = __str__

    #################
    # Class Methods #
    #################

    @classmethod
    def clear_known_ranges(cls):
        """ Clears the known ranges.  Use with care and only if you're sure you know what you're doing!
        """
        cls.global_index_range_set.known_ranges = {}
    reset_ranges = clear_known_ranges

    @classmethod
    def set_global_index_range_set(cls, new_global_set):
        cls._global_index_context_stack.append(cls.global_index_range_set)
        cls.global_index_range_set = new_global_set
    # Aliases illustrating the new way to think about IndexRangeSets as a transparent context
    set_indexing_context = set_global_index_range_set
    begin_indexing_context = set_global_index_range_set

    @classmethod
    def unset_global_index_range_set(cls, context_to_end=None):
        # check to make sure there's a context to take it's place
        if len(cls._global_index_context_stack) == 0:
            raise ValueError("No indexing context to end.")
        # Allow the user to give a context as a safety check
        if context_to_end is not None and cls.global_index_range_set is not context_to_end:
            raise ValueError("Requested context to end is not the current indexing context."
                             "  Perhaps you forgot to end another context created inside of this one?")
        cls.global_index_range_set = cls._global_index_context_stack.pop()
    # Alias illustrating the new way to think about IndexRangeSets as a transparent context
    end_indexing_context = unset_global_index_range_set


    ###########
    # Methods #
    ###########

    def add_subrange(self, subrange):
        """ Add a subrange to the list of subranges for `self`.
        Note that this returns `subrange` (with `subrange.parent` modified) to allow for a 'pass-through' like usage
        (see examples below)

        Examples
        --------
        >>> IndexRange.reset_ranges()
        >>> p = IndexRange('p,q,r,s', 4, name="orbital space")
        >>> # Convenient "pass-through" usage for saving subranges:
        >>> i = p.add_subrange(IndexRange('i,j,k,l', 0, 2, name="occupied space"))
        >>> a = p.add_subrange(IndexRange('a,b,c,d',2,4))
        >>> a
        <IndexRange object covering slice 2:4 represented by indices ['a', 'b', 'c', 'd']>
        >>> p.subranges[0]
        <IndexRange named 'occupied space' represented by indices ['i', 'j', 'k', 'l']>
        >>> p.subranges[1]
        <IndexRange object covering slice 2:4 represented by indices ['a', 'b', 'c', 'd']>
        >>> a.parent
        <IndexRange named 'orbital space' represented by indices ['p', 'q', 'r', 's']>

        """
        if isinstance(subrange, IndexRange):
            if not subrange._begin_is_ellipsis and subrange.begin_index < self.begin_index:
                raise IndexError("Subrange falls outside of parent range:  Subrange start ("
                                    + str(subrange.begin_index) + ") is before parent range start ("
                                    + str(self.begin_index) + ")"
                                )
            elif not subrange._end_is_ellipsis and subrange.end_index > self.end_index:
                raise IndexError("Subrange falls outside of parent range:  Subrange end ("
                                 + str(subrange.begin_index) + ") is after parent end ("
                                 + str(self.begin_index) + ")"
                )
            else:
                self.subranges.append(subrange)
                subrange.parent = self
                # Pass though for easy assigning...
                return subrange
        else:
            raise TypeError("subrange must be an instance of IndexRange")

    def with_subranges(self, *subranges):
        """ Basically the same thing as calling `add_subrange()` multiple times, except returns `self` instead of
        the subrange, allowing for a different syntax (see below)

        Examples
        --------
        >>> IndexRange.reset_ranges()
        >>> orb = DeclareIndexRange('p,q,r,s', 10, name="Orbital space").with_subranges(
        ...           DeclareIndexRange('i,j,k,l', 0, 3, name="Occupied space").with_subranges(
        ...               DeclareIndexRange("i*,j*,k*,l*", 0, 1, name="Core"),
        ...               DeclareIndexRange("i',j',k',l'", 1, 3)
        ...           ),
        ...           DeclareIndexRange('a,b,c,d', 3, 10, name="Virtual space")
        ...       )
        >>> orb
        <IndexRange named 'Orbital space' represented by indices ['p', 'q', 'r', 's']>
        >>> len(orb.subranges)
        2
        >>> len(orb.subranges[0].subranges)
        2

        """
        for subrange in listify_args(subranges):
            self.add_subrange(subrange)
            # TODO raise an error if subrange's index_range_set is set to something other than the same thing as ours
        return self

    def slice_in(self, parent_range):
        """ Gets the slice of `parent_range` represented by self
        (`parent_range` need not be a *direct* parent of self, but it should be a parent.  See `is_subrange()`)

        Examples
        --------
        >>> IndexRange.reset_ranges()
        >>> p = IndexRange('p,q,r,s',4)
        >>> i = p.add_subrange(IndexRange('i,j,k,l',0,2))
        >>> a = p.add_subrange(IndexRange('a,b,c,d',2,4))
        >>> a.slice_in(p)
        slice(2, 4, None)

        """
        start = self.begin_index - parent_range.begin_index
        return slice(start, start + self.size)



#####################
# Dependent Imports #
#####################

from grendel.gmath.einsum import EinsumTensor

