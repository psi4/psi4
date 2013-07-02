"""
Abstract base classes that don't logically belong somewhere else.
"""
import inspect
import types
from abc import ABCMeta, abstractmethod
from collections import Iterable

class NonstringIterable(object):
    """ An `Iterable` that is not an instance of basestring
    """
    __metaclass__ = ABCMeta

    @classmethod
    def __subclasshook__(cls, subclass):
        if issubclass(subclass, basestring):
            return False
        else:
            return issubclass(subclass, Iterable)
NonStringIterable = NonstringIterable


class FunctionLike(object):
    """ Specifies the requirements for a class that grendel can comfortably assume is like a function.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def getargspec(self):
        return NotImplemented

    @abstractmethod
    def __call__(self, *args, **kwargs):
        return NotImplemented

class MethodLike(FunctionLike):
    """ Specifies the minimum requirements for a class that grendel can also assume is bindable like a method.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __get__(self, instance, owner):
        return NotImplemented


class EasilyCopied(object):
    """
    Abstract base class defining a simple copying protocol in which keyword arguments for the
    constructor in `__copy__` are built up hierarchally using the `__copy_kwargs__` method.
    """
    __metaclass__ = ABCMeta

    ###################
    # Special Methods #
    ###################

    def __copy__(self):
        return self.__class__(**self.__copy_kwargs__())

    @abstractmethod
    def __copy_kwargs__(self):
        return dict()

class EasilyDeepCopied(object):
    """
    Abstract base class defining a simple deep copying protocol in which keyword arguments for the
    constructor in `__copy__` are built up hierarchally using the `__copy_kwargs__` method.
    """
    __metaclass__ = ABCMeta

    ###################
    # Special Methods #
    ###################

    def __deepcopy__(self):
        return self.__class__(**self.__deepcopy_kwargs__())

    @abstractmethod
    def __deepcopy_kwargs__(self):
        return dict()

