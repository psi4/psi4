""" Miscellaneous metaclasses that don't logically belong somewhere else.
"""

# Nothing in this module thus far should be exposed to the public interface.

__all__ = []

#region | Immutable |

class Immutable(type):
    """
    Metaclass for classes intended to be immutable.  Immutable classes can only
    have attributes that are themselves immutable (or built-in types known to be
    immutable, like str, int, float, complex, tuple, etc.).  Attributes can only
    be set between the beginning and the end of the __init__ method (i.e. methods
    called by __init__ can also make changes, but only if they are called from
    __init__).  Methods intended to be called from __init__ that mutate the class
    should have the @mutator decorator (though forgetting to do so will only
    lead to a more confusing error message now, it may be used later).
    All base classes of an Immutable class must themselves be immutable
    """

    IMMUTABLE_BUILTINS = (
        str,
        int,
        long,
        float,
        complex,
        bool,
        slice
    )

    IMMUTABLE_SEQUENCES = (
        frozenset,
        tuple
    )

    #TODO support for lazily evaluated attributes

    def __new__(mcs, name, bases, cls_dict):
        # Make sure the bases are immutable
        for base in bases:
            if not isinstance(base, Immutable) and base is not object:
                raise TypeError("Immutable classes can only inherit from other Immutable classes"
                                " or 'object';  '{0}' inherits from '{1}', which is not Immutable".format(
                    name, base.__name__
                ))
        #----------------------------------------#
        # Patch __init__ to activate the immutability guard after it finishes
        if "__init__" in cls_dict:
            init_func = cls_dict['__init__']
        else:
            init_func = None
        def __init__(self, *args, **kwargs):
            object.__setattr__(self, "__allow_mutations__", True)
            if init_func is not None:
                init_func(self, *args, **kwargs)
            if self.__class__.__name__ == name:
                # only stop allowing mutations if we are in our
                #   own __init__
                self.__allow_mutations__ = False
        if init_func is not None:
            __init__.__doc__ = init_func.__doc__
            __init__.__dict__.update(init_func.__dict__)
        cls_dict['__init__'] = __init__
        #----------------------------------------#
        # Immutability checker
        def is_immutable(obj):
            if isinstance(obj, Immutable.IMMUTABLE_BUILTINS) or isinstance(type(obj), Immutable):
                return True
            elif isinstance(obj, Immutable.IMMUTABLE_SEQUENCES):
                return all(is_immutable(item) for item in obj)
            else:
                return False
        #----------------------------------------#
        # Make __setattr__ raise an error
        if "__setattr__" in cls_dict:
            raise TypeError("Immutable classes cannot define '__setattr__'")
        def __setattr__(self, key, value):
            if self.__allow_mutations__:
                # This means we're coming from __init__ or some function that __init__ calls
                # But we still need to check if the value is immutable
                if not is_immutable(value):
                    raise AttributeError("Attributes of Immutable objects must be themselves Immutable.")
                object.__setattr__(self, key, value)
            else:
                raise TypeError("Can't change attribute '{0}' of object of immutable type '{1}'".format(
                    key, name
                ))
        cls_dict['__setattr__'] = __setattr__
        #----------------------------------------#
        # Make __delattr__ raise an error
        if "__delattr__" in cls_dict:
            raise TypeError("Immutable classes cannot define '__delattr__'")
        def __delattr__(self, item):
            if self.__allow_mutations__:
                object.__delattr__(self, item)
            else:
                raise TypeError("Can't delete attribute '{0}' of object of immutable type '{1}'".format(
                    item, name
                ))
        cls_dict['__delattr__'] = __delattr__
        #----------------------------------------#
        return type.__new__(mcs, name, bases, cls_dict)

def mutator(func):
    """
    Does nothing so far
    """
    func.__is_mutator__ = True
    return func

#TODO StrictlyImmutable submetaclass that checks hash before and after every method call?
#endregion

#--------------------------------------------------------------------------------#

#region | SubscriptableClass                                                             |

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

#endregion

#--------------------------------------------------------------------------------#

#region | Commented out code                                                             |

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
#endregion

#####################
# Dependent Imports #
#####################

from grendel.util.decorators import with_attributes

