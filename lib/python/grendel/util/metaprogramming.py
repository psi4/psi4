""" Various utility functions that don't logically belong anywhere else for generating parts of classes or functions
"""

# Nothing yet that should be exposed to the public interface
from abc import ABCMeta, abstractmethod
from collections import Iterable
from grendel import type_checking_enabled
import types

__all__ = []

# TODO ClassAttribute

# TODO Reorganize this.  Multiple "action" attributes should be constructable (with little extra work) simply by nesting.
class SpecialAttribute(object):
    """ Abstract superclass for special attribute descriptor
    """
    __metaclass__ = ABCMeta

    private_name = None
    default = None

    def __init__(self, name, doc=None, private_name=None, default=None):
        """ Default constructor for special attributes.
        The attribute is stored as `private_name` in the parent class, which defaults to `'_' + name`.
        If `default` is given, it acts as the class value of the attribute as well as the default value
        when `private_name` is not defined for an instance.
        """
        self.__name__ = name
        self.private_name = private_name
        if self.private_name is None is None:
            self.private_name = '_' + name
        self.__doc__ = doc
        self.default = default

    @abstractmethod
    def __attr_init__(self, *args):
        return NotImplemented

    @abstractmethod
    def __set_check__(self, instance, value):
        """ Called before setting an attribute.  The attribute set only proceeds if __set_check__ returns True.
        This allows attributes with multiple classes to be created.
        """
        return NotImplemented

    def __get__(self, instance, owner):
        if instance is None:
            return self.default
        else:
            return getattr(instance, self.private_name, self.default)

    def __set__(self, instance, value):
        if self.__set_check__(instance, value):
            setattr(instance, self.private_name, value)

class ReadOnlyAttribute(SpecialAttribute):
    """ A read-only attribute.
    See `SpecialAttribute` for initialization variables standard to all `SpecialAttribute` subclasses.

    :Attributes:

    error_text : str
        The text for the error raised when setting of the attribute is attempted.  Defaults to "cannot set read-only attribute"

    """

    error_text = "cannot set read-only attribute"

    def __init__(self, name, doc=None, private_name=None, default=None, error_text=None):
        super(ReadOnlyAttribute, self).__init__(name, doc, private_name, default)
        self.__attr_init__(error_text)

    def __attr_init__(self, error_text):
        if error_text:
            self.error_text = error_text

    def __set_check__(self, instance, value):
        raise AttributeError("cannot set read-only attribute")


# TODO migrate documentation
# TODO Make this simpler...there's no reason something should be both read-only and type_checked...
class TypeCheckedAttribute(SpecialAttribute):
    """An Attribute that raises a TypeError if it is set to some `new_val` that returns
    `False` for `is_instance(new_val, types)`.  If `strong` is True, the property created will
    run the typecheck even if `grendel.type_checking_enabled` is `False` (don't do this unless
    you have a very good reason to).  If `type_error_text` is given, it is used for the text of the TypeError to be
    raised when an incorrect type is given.  Aliased as `type_checked`.
    See `read_only_attribute` for explaination of other arguments.

    :Examples:

    >>> class SuperSafe(object):
    ...     name = TypeCheckedAttribute('name', str)
    ...     number = TypeCheckedAttribute('number', (int, float, long), default=0.0)
    ...     def __init__(self):
    ...         self._name = 'hello world'
    ...
    >>> foo = SuperSafe()
    >>> foo.name
    'hello world'
    >>> foo.name = 'goodbye'
    >>> foo.name
    'goodbye'
    >>> foo.number
    0.0
    >>> foo.number = 25
    >>> foo.number
    25
    >>> foo.name = 75
    Traceback (most recent call last):
    ...
    TypeError: Attribute 'name' in class 'SuperSafe' must be an instance of 'str' (got 'int')
    >>> foo.number = 'bar'
    Traceback (most recent call last):
    ...
    TypeError: Attribute 'number' in class 'SuperSafe' must be an instance of one of the following: ('int', 'float', 'long'); (got 'str')

    """

    type_error_text = None
    strong = False

    def __init__(self, name, types, doc=None, private_name=None, default=None, strong=False, type_error_text=None):
        super(TypeCheckedAttribute, self).__init__(name, doc, private_name, default)
        self.__attr_init__(types, strong, type_error_text)

    def __attr_init__(self, types, strong, type_error_text):
        self.types = tuple(listify_args(types))
        if type_error_text is None:
            self.type_error_text = "Attribute '{name}'".format(name=self.__name__)
            self.type_error_text += " in class '{cls}' "
            if isinstance(types, Iterable):
                second_part = "must be an instance of one of the following: ('{types}');".format(types="', '".join(map(classname, types)))
            else:
                second_part = "must be an instance of '{types}'".format(types=classname(types))
            self.type_error_text += second_part + " (got '{got}')"
        else:
            self.type_error_text = type_error_text
        self.strong = strong

    def __set_check__(self, instance, value):
        if type_checking_enabled:
            if not isinstance(value, self.types):
                raise TypeError(self.type_error_text.format(cls=classname(type(instance)), got=classname(type(value))))
        return True


#--------------------------------------------------------------------------------#
#                               Importable Class                                 #
#--------------------------------------------------------------------------------#

class ImportableClass(types.ModuleType):

    def __import__(self, globals=None, locals=None, fromlist=None, level=None):
        return self


########################
# Aliases and Wrappers #
########################

def read_only_attribute(*args, **kwargs):
    return ReadOnlyAttribute(*args, **kwargs)



#####################
# Dependent Imports #
#####################

from grendel.util.strings import classname
from grendel.util.overloading import listify_args

