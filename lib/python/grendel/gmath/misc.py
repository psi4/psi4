""" Miscellaneous mathematical utility functions
"""
from grendel import type_checking_enabled
from grendel.util.decorators import typechecked, SequenceOf, AnyType
from grendel.util.metaclasses import SubscriptableClass

class delta(object):
    __metaclass__ = SubscriptableClass

    @classmethod
    def __class_getitem__(cls, item):
        if type_checking_enabled and not isinstance(item, tuple):
            raise TypeError("delta must be subscripted with at least two arguments")
        return 1 if all(i == item[0] for i in item[1:]) else 0

    def __new__(cls, *args):
        return cls[args]

def get_permutation(relative_to, permuted):
    """
    Returns the indices of `permuted` in `relative_to`
    """
    relative_to = tuple(relative_to)
    permuted = tuple(permuted)
    srel = set(relative_to)
    if srel != set(permuted):
        raise ValueError("'{}' is not a permutation of '{}'".format(permuted, relative_to))
    if len(srel) != len(relative_to):
        raise ValueError("Permutation of '{}', containing repeats, is ambiguous.".format(relative_to))
    return tuple(relative_to.index(item) for item in permuted)

def permutation_is_even(perm):
    p = 0
    v = [0] * (len(perm) + 1)
    j = len(perm)
    while j > 0:
        j -= 1
        if v[j] != 0:
            p += 1
        else:
            x = perm[j]
            v[x] = 1
            while x != j:
                x = perm[x]
                v[x] = 1
    return (p & 1) == 0

