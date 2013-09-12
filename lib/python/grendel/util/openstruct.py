from grendel.util.strings import classname

class DescriptorCallingOpenStruct(type):
    """ Allows the assignment of arbitrary attributes to an instance without predeclaring their existance.
    This particular open struct treats descriptors as though they were part of the parent class dictionary
    (see http://docs.python.org/reference/datamodel.html#implementing-descriptors).
    (Note: most of this is trivially available for all instances in Python, but things like iteration and
    the `CaseInsensativeOpenStruct` subclass are not; this is why this is here).

    :Examples:

    >>> hello = DescriptorCallingOpenStruct()
    >>> hello.world = 'foo'
    >>> hello.moon = 'bar'
    >>> print hello.world + hello.moon
    foobar
    >>> # Iteration: (notice that order is not guarenteed)
    >>> for attr, val in hello:
    ...     print attr + ' => ' + val
    moon => bar
    world => foo
    >>> # It inherits from type, but it can't act as a metaclass...
    >>> class CannotBe(object):
    ...     __metaclass__ = DescriptorCallingOpenStruct
    ...     pass
    ...
    Traceback (most recent call last):
        ...
    TypeError: 'DescriptorCallingOpenStruct' is not a real metaclass.  It merely inherits from 'type' to make descriptors work properly.


    """

    ##################
    # Initialization #
    ##################

    def __new__(mcs, *args, **kwargs):
        ret_val = type.__new__(mcs, classname(mcs), (), {})
        mcs.__init__(ret_val, *args, **kwargs)
        return ret_val

    def __init__(cls, *args, **dic):
        if not len(args) == 0:
            if len(args) == 3 and all(isinstance(a, t) for a, t in zip(args, [basestring, tuple, dict])):
                raise TypeError("'{0}' is not a real metaclass.  It merely inherits from 'type' to "
                                "make descriptors work properly.".format(classname(cls)))
            else:
                raise TypeError("{1} constructor takes exactly 1 argument, besides keyword arguments"
                                " ({0} given)".format(len(args) + 1, classname(cls)))
        for item in dic:
            setattr(cls, item, dic[item])

    ###################
    # Special Methods #
    ###################

    def __getattr__(cls, item):
        if item in cls.__dict__:
            ret_val = cls.__dict__[item]
            # respect descriptors
            if hasattr(ret_val, '__get__'):
                return ret_val.__get__(None, cls)
            else:
                return ret_val
        else:
            raise AttributeError(item)

    def __setattr__(cls, key, value):
        if key in cls.__dict__:
            # respect descriptors
            if hasattr(cls.__dict__[key], '__set__'):
                cls.__dict__[key].__set__(cls, value)
            else:
                super(DescriptorCallingOpenStruct, cls).__setattr__(key, value)
        else:
            super(DescriptorCallingOpenStruct, cls).__setattr__(key, value)
        return value

    def __iter__(cls):
        for key in cls.__dict__:
            if key not in cls.__class__.__dict__ and key not in type.__dict__ and key not in ['__weakref__']:
                yield key, cls.__dict__[key]


class CaseInsensativeOpenStruct(DescriptorCallingOpenStruct):
    """ Acts like an open struct, but attribute access is case-insensative.

    :Examples:

    >>> hello = CaseInsensativeOpenStruct()
    >>> hello.World = "foo"
    >>> hello.mOOn = "bar"
    >>> print hello.woRLD + hello.MOON
    foobar
    >>> # Iteration yields the original names as keys along with the values
    >>> for attr, val in hello:
    ...     print attr + " => " + val
    ...
    mOOn => bar
    World => foo

    """

    ##############
    # Attributes #
    ##############

    __original_names__ = None


    ##################
    # Initialization #
    ##################

    def __init__(cls, *args, **dic):
        super(CaseInsensativeOpenStruct, cls).__init__(*args, **(dict(zip(map(str.lower, dic.keys()), dic.values()))))
        cls.__original_names__ = {}

    ###################
    # Special Methods #
    ###################

    def __getattr__(cls, item):
        if item.lower() in cls.__dict__:
            return super(CaseInsensativeOpenStruct, cls).__getattr__(item.lower())
        else:
            raise AttributeError(item)

    def __setattr__(cls, key, value):
        if key.lower() in cls.__dict__:
            super(CaseInsensativeOpenStruct, cls).__setattr__(key.lower(), value)
        else:
            super(CaseInsensativeOpenStruct, cls).__setattr__(key.lower(), value)
            cls.__original_names__[key.lower()] = key
        return value

    def __iter__(cls):
        for key in cls.__dict__:
            if key not in cls.__class__.__dict__ and key not in type.__dict__ and key not in ['__weakref__']:
                yield cls.__original_names__[key], cls.__dict__[key]

#####################
# Dependent Imports #
#####################

