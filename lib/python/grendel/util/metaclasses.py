""" Miscellaneous metaclasses that don't logically belong somewhere else.
"""

# Nothing in this module thus far should be exposed to the public interface.
from abc import ABCMeta
from functools import wraps
from itertools import groupby

__all__ = []

# Slightly more efficient than using @operable_class, but much less general
class SubscriptableClass(type):
    """ Meta-class that allows a class to define "classmethod" equivalents of `__getitem__` and `__setitem__`
    called `__class_getitem__` and `__class_setitem__`
    """

    def __getitem__(cls, item):
        if '__class_getitem__' in cls.__dict__:
            return cls.__class_getitem__(item)
        else:
            # Otherwise, act like we're not here
            return NotImplemented

    def __setitem__(cls, key, value):
        if '__class_setitem__' in cls.__dict__:
            return cls.__class_setitem__(key, value)
        else:
            # Otherwise, act like we're not here
            return NotImplemented

# TODO @must_call_super decorator
#def semi_abstract_method(group):
#    return with_attributes(__semiabstractgroup__=group)
#
#
#class AbstractMeta(ABCMeta):
#    """ A metaclass with a few customizations to and improvements on `abc.ABCMeta`.
#    You should still use `abc.ABCMeta` unless you specifically need some of the features implemented here,
#    since the built-in library abstract class mechanism is much faster.
#    """
#
#    def __new__(mcs, name, bases, namespace):
#        cls = super(AbstractMeta, mcs).__new__(mcs, name, bases, namespace)
#        cls.__semiabstracts__ = {}
#        for key, grp in groupby(namespace.items(), lambda pair: getattr(pair[1], '__semiabstractgroup__', None)):
#            if key is None: continue
#            names = frozenset([0] for pair in grp)
#            try:
#                cls.__semiabstracts__[key] = names
#            except TypeError:
#                raise TypeError("key '{0}' with unhashable type '{1}' cannot be used to define a group "
#                                "of semiabstract methods.".format(key, type(key).__name__))
#        for base in bases:
#            for key, names in getattr(base, '__semiabstracts__', {}).items():
#                pass
#        cls.__abstract_meta_oldnew = cls.__new__
#        # TODO replace with the wrapper from the decorator module that preserves function signature
#        @wraps(cls.__new__)
#        def new_with_check(cls, *args, **kwargs):
#            for base in cls.__mro__:
#                for key, names in getattr(base, '__semiabstracts__', {}):
#                    if not any(name in cls.__dict__ for name in names):
#                        # Missing from cls...check the parent classes up to base...
#                        for base_to_check in cls.__mro__[1:]:
#                            # TODO Finsih this
#                            # also: use inspect.getclasstree to distinguish siblings from parents
#
#
#
#
#
#
#        return cls

#####################
# Dependent Imports #
#####################

from grendel.util.decorators import with_attributes

