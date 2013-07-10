"""
Protocol which, if carefully implemented for a class, allows instances of the class to be "frozen", or at least
act as such to some extent or another.

:Examples:

>>> class Foo(Freezable):
...     bar = FreezableAttribute('bar')
...     lst = FreezableListAttribute('lst')
...     once = SetOnceAttribute('once')
...     def __init__(self):
...         self.bar = None
...         self.lst = []
...     def __freeze__(self):
...         print "it's getting cold in here..."
...
>>> f = Foo()
>>> f.bar = 'hello'
>>> f.once = 'something'
>>> f.lst.append('Bart')
>>> f.lst.extend(['Marge', 'Lisa', 'Homer', 'Maggie'])
>>> f.lst
['Bart', 'Marge', 'Lisa', 'Homer', 'Maggie']
>>> f.once = 'something else' # doctest: +ELLIPSIS
Traceback (most recent call last):
    ...
FrozenObjectError: Can't modify attribute 'once' of instance <...> of type Foo once it has been set.
>>> f.bar
'hello'
>>> freeze(f)
it's getting cold in here...
>>> f.bar
'hello'
>>> f.lst
('Bart', 'Marge', 'Lisa', 'Homer', 'Maggie')
>>> f.bar = 'world'  # doctest: +ELLIPSIS
Traceback (most recent call last):
    ...
FrozenObjectError: Can't modify attribute 'bar' of instance <...> of type Foo once the instance has been frozen.
>>> f.lst.append("Ned")
Traceback (most recent call last):
    ...
AttributeError: 'tuple' object has no attribute 'append'


"""
from collections import defaultdict

from grendel import type_checking_enabled
from grendel.util.overloading import caller
from grendel.util.strings import classname

class FrozenObjectError(RuntimeError):
    """ Raised when a frozen class is modified.
    """

def freeze(obj):
    if type_checking_enabled and not isinstance(obj, Freezable):
        raise TypeError("Cannot freeze unfreezable object of type {0}".format(classname(obj)))
    else:
        obj.freeze()

class Freezable(object):
    """ Superclass for classes whose instances can become immutable at some point in time.
    """

    frozen = False

    def freeze(self):
        if hasattr(self, '__freeze__'):
            self.__freeze__()
        for item in self.__class__.__dict__:
            if hasattr(self.__class__.__dict__[item], '__freeze__'):
                self.__class__.__dict__[item].__freeze__(self)
        self.frozen = True

    def fail_if_frozen(self):
        if self.frozen:
            raise FrozenObjectError("Can't call method {0}() on instance {1} of {2} once "
                                                "it has been frozen.".format(caller(), self, classname(self)))

class FreezableAttribute(object):

    _name = None
    _private_name = None
    _default = None
    _setter_func = None

    def __init__(self, name, default=None, doc=None):
        self.__doc__ = doc
        self._name = name
        self._private_name = "_" + self._name
        self._default = default

    def __get__(self, instance, owner):
        if instance is None:
            return self._default
        else:
            return getattr(instance, self._private_name, self._default)

    def __set__(self, instance, value):
        if getattr(instance, 'frozen', False):
            raise FrozenObjectError("Can't modify attribute '{0}' of instance {1} of type {2} once the "
                                    "instance has been frozen.".format(self._name, instance, classname(instance)))
        elif self._setter_func is not None:
            self._setter_func(instance, value)
        else:
            setattr(instance, self._private_name, value)

    def setter(self, f):
        self._setter_func = f
        return self


#TODO typecheck for list in __set__
class FreezableListAttribute(FreezableAttribute):

    def __freeze__(self, instance):
        setattr(instance, self._private_name, tuple(getattr(instance, self._private_name, self._default)))


class SetOnceAttribute(FreezableAttribute):

    suffix = '__SetOnceAttribute_isset'

    required_before_freeze = None

    def __init__(self, name, default=None, doc=None, required_before_freeze=False):
        self.required_before_freeze = required_before_freeze
        super(SetOnceAttribute, self).__init__(name, default, doc)

    def __freeze__(self, instance):
        if self.required_before_freeze:
            if not self.is_set_for_instance(instance):
                raise FrozenObjectError("Attribute '{0}' of instance {1} of type {2} must be set before "
                                        "freezing can occur.".format(self._name, instance, classname(instance)))
        setattr(instance, self._private_name + SetOnceAttribute.suffix, True)

    def __set__(self, instance, value):
        if self.is_set_for_instance(instance):
            raise FrozenObjectError("Can't modify attribute '{0}' of instance {1} of type {2} once "
                                    "it has been set.".format(self._name, instance, classname(instance)))
        else:
            super(SetOnceAttribute, self).__set__(instance, value)
            if value is not None and value is not NotImplemented:
                setattr(instance, self._private_name + SetOnceAttribute.suffix, True)

    def is_set_for_instance(self, instance):
        return getattr(instance, self._private_name + SetOnceAttribute.suffix, False)


