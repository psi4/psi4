""" Miscellaneous tools for iteration and acting on iterables that don't belong somewhere else.
"""
from __future__ import print_function
from collections import Iterable, Iterator
from itertools import izip_longest, chain, tee, cycle, islice, combinations_with_replacement
import operator

from grendel.external.combinatorics import unlabeled_balls_in_unlabeled_boxes, labeled_balls_in_unlabeled_boxes, unlabeled_balls_in_labeled_boxes, labeled_balls_in_labeled_boxes

#--------------------------------------------------------------------------------#
#                                   Iterators                                    #
#--------------------------------------------------------------------------------#
from grendel.util.sentinal_values import ArgumentNotGiven

def stutter(iterable, ntimes):
    """ repeat each element n times before moving on.

    :Examples:

    >>> [s for s in stutter([1, 2, 3, 4], 3)]
    [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
    >>> ''.join(a for a in stutter('hello', 2))
    'hheelllloo'

    """
    return chain(*zip(*tee(iterable, ntimes)))

def flattened(iterable, keep_types=None, debug=False):
    """ Flatten across all levels of iteration.

    :Examples:

    >>> [f for f in flattened([1,2,3,[4,5,[6]]])]
    [1, 2, 3, 4, 5, 6]
    >>> [f for f in flattened([1,2,3,[4,5,[6],[[xrange(1,17)]]]])]
    [1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    >>> [f for f in flattened([range(1,4),1,2,3,[4,5,[6],[[xrange(1,17)]]]])]
    [1, 2, 3, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    >>> [f for f in flattened([1, [3, 4], [["hello", ["world"]]], {"yes": True}])]
    [1, 3, 4, 'h', 'e', 'l', 'l', 'o', 'w', 'o', 'r', 'l', 'd', 'y', 'e', 's']
    >>> [f for f in flattened([1, [3, 4], [["hello", ["world"]]], {"yes": True}], keep_types=str)]
    [1, 3, 4, 'hello', 'world', 'yes']
    >>> [f for f in flattened([1, [3, 4], [["hello", ["world"]]], {"yes": True}], keep_types=(str, dict))]
    [1, 3, 4, 'hello', 'world', {'yes': True}]
    >>> [f for f in flattened([1, [3, 4], [["hello", ["world"]]], {"yes": True}], keep_types=(str, dict, list))]
    [[1, [3, 4], [['hello', ['world']]], {'yes': True}]]

    """
    keep_types = tuple(keep_types if isinstance(keep_types, Iterable) else [keep_types]) if keep_types else ()
    curr = iterable
    if isinstance(curr, keep_types):
        yield curr
        raise StopIteration
    stack = []
    num_to_pop = []
    first_time = True
    while first_time or len(stack) > 0:
        if first_time: first_time = False
        if debug:
            print("curr is {}\nstack is {}\nnum_to_pop is {}\n----------------------------------".format(
                str(curr), str(stack), str(num_to_pop)))
        if isinstance(curr, Iterator):
            if num_to_pop[-1] == 1:
                num_to_pop.append(1)
            stack.append(curr)
            try:
                curr = next(curr)
            except StopIteration:
                stack.pop()
                n = num_to_pop.pop()
                if n == 2:
                    stack.pop()
                try:
                    curr = stack.pop()
                except IndexError:
                    raise StopIteration
        elif isinstance(curr, Iterable):
            if len(stack) >= 2 and curr is stack[-2]:
                # it's an object like a single-character string, whose first element of it's iterator is itself
                yield curr
                stack.pop()
                n = num_to_pop.pop()
                if n == 2:
                    stack.pop()
                try:
                    curr = stack.pop()
                except IndexError:
                    raise StopIteration
            elif len(keep_types) > 0 and isinstance(curr, keep_types):
                yield curr
                curr = stack.pop()
            else:
                stack.append(curr)
                curr = iter(curr)
                num_to_pop.append(2)
        else:
            yield curr
            curr = stack.pop()


def ordered_partitions_iter(sequence, length):
    """
    This iterates over the P^k,m operator from Allen, et al. Mol. Phys. 89 (1996), 1213-1221
    See the explanation of its funtion therein.  This is needed for the arbitrary order B
    tensor formulae.

    :Examples:


    >>> [tuple(''.join(part) for part in parts) for parts in ordered_partitions_iter('ABCD', 2)]
    [('A', 'BCD'), ('AB', 'CD')]
    >>> [tuple(''.join(part) for part in parts) for parts in ordered_partitions_iter('ABCDEF', 3)]
    [('A', 'B', 'CDEF'), ('A', 'BC', 'DEF'), ('AB', 'CD', 'EF')]


    """
    n = len(sequence)
    for partitions in unlabeled_balls_in_unlabeled_boxes(n, [n]*length):
        if 0 in partitions:
            continue
        spart = sorted(partitions)
        yield tuple(partitioned(sequence, spart))

def all_partitions_iter(sequence, length, allow_zero=True):
    """
    Iterate over all possible partitions of a sequence

    :Examples:


    """
    n = len(sequence)
    for partitions in unlabeled_balls_in_labeled_boxes(n, [n]*length):
        if not allow_zero and 0 in partitions:
            continue
        to_yield = tuple(partitioned(sequence, partitions))[:length]
        yield to_yield

def ordered_partitions(sequence, length):
    """
    see `ordered_partitions_iter`
    """
    return [p for p in ordered_partitions_iter(sequence, length)]

def brace_notation_iter(sequence_of_sequences):
    """
    This iterates over the brace notation combinations from Allen, et al. Mol. Phys. 89 (1996), 1213-1221
    See the explanation of its funtion therein.  This is needed for the arbitrary order B
    tensor formulae.

    :Examples:

    >>> [tuple(''.join(part) for part in parts) for parts in brace_notation_iter(['AB', 'C', 'D'])]
    [('AB', 'C', 'D'), ('AC', 'B', 'D'), ('AD', 'B', 'C'), ('BC', 'A', 'D'), ('BD', 'A', 'C'), ('CD', 'A', 'B')]

    """
    seq = [p for p in sequence_of_sequences]
    joined = [item for item in chain(*seq)]
    for indices in labeled_balls_in_unlabeled_boxes(len(joined), [len(s) for s in seq]):
        yield tuple(tuple(joined[idx] for idx in subset) for subset in indices)

def partitioned(iterable, partitions):
    """ Iterate over iterable in chunks of sizes given in partitions.
    If `partitions` runs out before `iterable` does, `partitions` will
    be cycled (i.e. the next partition size will be the first element of `partitions`,
    the following size will be the second, and so on).  If there are not
    enough elements remaining to fill the last partition, the
    last item yielded will be shorter (just as if a slice running over
    the end of the list were taken).

    :Examples:


    >>> [''.join(a) for a in partitioned('ABCDEFG', (1, 3, 2, 1))]
    ['A', 'BCD', 'EF', 'G']
    >>> [''.join(a) for a in partitioned('ABCDEFGHIJ', (1, 3, 2, 1))]
    ['A', 'BCD', 'EF', 'G', 'H', 'IJ']
    >>> [''.join(a) for a in partitioned('ABCDEFGHIJK', (2,))]
    ['AB', 'CD', 'EF', 'GH', 'IJ', 'K']
    >>> [''.join(a) for a in partitioned('ABC', (1, 3, 2, 1))]
    ['A', 'BC']
    >>> [i for i in partitioned('abcdef', (0,3))]
    [(), ('a', 'b', 'c'), (), ('d', 'e', 'f'), ()]
    >>> [i for i in partitioned('abcdef', (3,1,0,3))]
    [('a', 'b', 'c'), ('d',), (), ('e', 'f')]
    >>> [i for i in partitioned('abcdef', (3,1,0))]
    [('a', 'b', 'c'), ('d',), (), ('e', 'f')]

    """
    prev_part, this_part = tee(cummulative_sum(cycle(partitions)), 2)
    start_index, end_index = 0, next(this_part)
    num_yielded = 0
    to_yield = tuple(item for item in islice(iterable, 0, end_index))
    while not (len(to_yield) == 0 and start_index != end_index):
        yield to_yield
        num_yielded += 1
        start_index, end_index = next(prev_part), next(this_part)
        to_yield = tuple(item for item in islice(iterable, start_index, end_index))

def cummulative_sum(iterable, initial_value=0, op=operator.add):
    """ Iterate over the progressive cumulative sum of the items in the iterable.

    :Examples:


    >>> [i for i in cummulative_sum(xrange(6))]
    [0, 1, 3, 6, 10, 15]
    >>> [a for a in cummulative_sum('ABCD', '')]
    ['A', 'AB', 'ABC', 'ABCD']
    >>> [a for a in cummulative_sum(cummulative_sum('ABCD', ''), '')]
    ['A', 'AAB', 'AABABC', 'AABABCABCD']
    >>> [i for i in cummulative_sum(xrange(1,6), 1, op=operator.mul)]
    [1, 2, 6, 24, 120]


    """
    sum = initial_value
    for i in iterable:
        sum = op(sum, i)
        yield sum

#------------------------------------------------------------------------------------#
# Unique permutations iterator                                                       #
# Credit: http://stackoverflow.com/questions/6284396/permutations-with-unique-values #
# (slightly modified to fit our needs)                                               #
#------------------------------------------------------------------------------------#

class unique_element:
    def __init__(self,value,occurrences):
        self.value = value
        self.occurrences = occurrences

def unique_permutations(elements):
    eset=set(elements)
    listunique = [unique_element(i,elements.count(i)) for i in eset]
    u=len(elements)
    return perm_unique_helper(listunique,[0]*u,u-1)

def perm_unique_helper(listunique,result_list,d):
    if d < 0:
        yield tuple(result_list)
    else:
        for i in listunique:
            if i.occurrences > 0:
                result_list[d]=i.value
                i.occurrences-=1
                for g in  perm_unique_helper(listunique,result_list,d-1):
                    yield g
                i.occurrences+=1

#----------------------------------#
# END Unique permutations iterator #
#----------------------------------#

#---------------------------------------#
# Combinatorial (ul-notation) iterators #
#---------------------------------------#
# See Hollman and Schaefer paper (forthcoming)
def I_uu(nballs, nboxes, labels):
    if nballs == 0:
        yield (tuple(),) * nboxes
    else:
        for occs in unlabeled_balls_in_unlabeled_boxes(nballs, [nballs]*nboxes):
            stop_points = (0,) + tuple(cummulative_sum(occs))
            yield tuple(labels[stop_points[i]:stop_points[i+1]] for i in xrange(nboxes))

def I_uubar(nballs, nboxes, labels):
    if nballs == 0:
        yield (tuple(),) * nboxes
    else:
        for occs in unlabeled_balls_in_unlabeled_boxes(nballs, [nballs]*nboxes):
            if 0 in occs:
                continue
            stop_points = (0,) + tuple(cummulative_sum(occs))
            yield tuple(labels[stop_points[i]:stop_points[i+1]] for i in xrange(nboxes))

def I_ul(nballs, nboxes, labels):
    if nballs == 0:
        yield (tuple(),) * nboxes
    else:
        for occs in unlabeled_balls_in_labeled_boxes(nballs, [nballs]*nboxes):
            stop_points = (0,) + tuple(cummulative_sum(occs))
            yield tuple(labels[stop_points[i]:stop_points[i+1]] for i in xrange(nboxes))

def I_ulbar(nballs, nboxes, labels):
    if nballs == 0:
        yield (tuple(),) * nboxes
    else:
        for occs in unlabeled_balls_in_labeled_boxes(nballs, [nballs]*nboxes):
            if 0 in occs:
                continue
            stop_points = (0,) + tuple(cummulative_sum(occs))
            yield tuple(labels[stop_points[i]:stop_points[i+1]] for i in xrange(nboxes))

def I_lu(nballs, nboxes, labels):
    if nballs == 0:
        yield (tuple(),) * nboxes
    else:
        for idxs in labeled_balls_in_unlabeled_boxes(nballs, [nballs]*nboxes):
            yield tuple(tuple(labels[i] for i in subset) for subset in idxs)

def I_lubar(nballs, nboxes, labels):
    if nballs != 0:
        for idxs in labeled_balls_in_unlabeled_boxes(nballs, [nballs]*nboxes):
            if not all(len(i) > 0 for i in idxs):
                continue
            else:
                yield tuple(tuple(labels[i] for i in subset) for subset in idxs)

def I_ll(nballs, nboxes, labels):
    if nballs == 0:
        yield (tuple(),) * nboxes
    else:
        for idxs in labeled_balls_in_labeled_boxes(nballs, [nballs]*nboxes):
            yield tuple(tuple(labels[i] for i in subset) for subset in idxs)

def I_llbar(nballs, nboxes, labels):
    if nballs != 0:
        for idxs in labeled_balls_in_labeled_boxes(nballs, [nballs]*nboxes):
            if not all(len(i) > 0 for i in idxs):
                continue
            else:
                yield tuple(tuple(labels[i] for i in subset) for subset in idxs)


#--------------------------------------------------------------------------------#
#                  Functions that act on iterables                               #
#--------------------------------------------------------------------------------#

def first(iterable, predicate=lambda x: True, default_value=ArgumentNotGiven):
    for item in iterable:
        if predicate(item):
            return item
    if default_value is not ArgumentNotGiven:
        return default_value
    else:
        raise ValueError("item not found")

def last(iterable, predicate=lambda x: True, default_value=ArgumentNotGiven):
    for item in reversed(iterable):
        if predicate(item):
            return item
    if default_value is not ArgumentNotGiven:
        return default_value
    else:
        raise ValueError("item not found")

#--------------------------------------------------------------------------------#
#     Copied from the recipes section of the itertools module documentation      #
#--------------------------------------------------------------------------------#

def grouper(n, iterable, fillvalue=None):
    """group(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"""
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


