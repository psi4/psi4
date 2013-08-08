import inspect
from grendel.util.aliasing import bindablepartial

def operable_class(*ops):
    """ Decorator to make class work with certain 'magic' double-underscored methods
    (specified as arguments *without* the proceeding and trailing double-underscores).
     This only works with new-style classes, as far as I can tell.


    :Examples:

    >>> @operable_class('or', 'and')
    ... class Awesome(object):
    ...     def __or__(self, other):
    ...         print "instance or"
    ...     @classmethod
    ...     def __class_or__(cls, other):
    ...         print "class or"
    ...     @classmethod
    ...     def __class_and__(cls, other):
    ...         print "class and"
    ...
    >>> Awesome|Awesome
    class or
    >>> a, b = Awesome(), Awesome()
    >>> a|b
    instance or
    >>> Awesome|b
    class or
    >>> a|Awesome
    instance or
    >>> Awesome & Awesome
    class and
    >>> a & b
    Traceback (most recent call last):
        ...
    TypeError: unsupported operand type(s) for &: 'Awesome' and 'Awesome'

    *Subclassing:*

    >>> class MoreAwesome(Awesome):
    ...     @classmethod
    ...     def __class_and__(cls, other):
    ...         print "more awesome and"
    ...
    >>> MoreAwesome & Awesome
    more awesome and
    >>> Awesome & MoreAwesome
    class and
    >>> @operable_class('xor')
    ... class AwesomeXor(MoreAwesome):
    ...     @classmethod
    ...     def __class_xor__(cls, other):
    ...         print "class xor"
    ...
    >>> AwesomeXor ^ Awesome
    class xor
    >>> Awesome ^ AwesomeXor
    Traceback (most recent call last):
        ...
    TypeError: unsupported operand type(s) for ^: 'type' and 'type'

    *Metaclasses:*

    >>> class AwesomeMeta(type):
    ...     def __new__(mcs, name, bases, dict):
    ...         dict['foo'] = 'metacls was here'
    ...         return type.__new__(mcs, name, bases, dict)
    ...
    >>> @operable_class('add')
    ... class TheNewAwesome(object):
    ...     __metaclass__ = AwesomeMeta
    ...     @classmethod
    ...     def __class_add__(cls, other):
    ...         print "class add with foo = " + cls.foo
    ...     def __add__(self, other):
    ...         print "instance add"
    ...
    >>> TheNewAwesome + Awesome
    class add with foo = metacls was here
    >>> @operable_class('sub')
    ... class TheNewerAwesome(TheNewAwesome):
    ...     @classmethod
    ...     def __class_sub__(cls, other):
    ...         print "class sub with foo = " + cls.foo
    ...
    >>> TheNewerAwesome - TheNewAwesome
    class sub with foo = metacls was here
    >>> TheNewerAwesome + TheNewAwesome
    class add with foo = metacls was here
    >>> TheNewerAwesome() + a
    instance add

    :Raises:

    TypeError :  If the method `__class_<op>__` is not defined as a classmethod.  For instance,

        >>> @operable_class('mul')
        ... class Fail(object):
        ...     def __class_mul__(cls, obj):
        ...         print 'class mul'
        ...
        Traceback (most recent call last):
            ...
        TypeError: class magic method '__class_mul__' must be defined as a classmethod

    TypeError : If the thing to be decorated is not a class.  For example,

        >>> @operable_class('mul')
        ... def fail(obj): pass
        ...
        Traceback (most recent call last):
            ...
        TypeError: operable_class is a decorator for classes, not functions

    """
    # A __new__ for the automatically-created metaclass of the new class
    def _meta_new(mcs, cls_or_name, bases=None, dct=None):
        if inspect.isclass(cls_or_name):
            name = cls_or_name.__name__
            bases = inspect.getclasstree((cls_or_name,))[-1][0][1]
            # Respect the class's metaclass
            return type(mcs).__new__(mcs, name, bases, dict(cls_or_name.__dict__))
        else:
            # This allows subclasses of cls to be constructed normally...
            return type(mcs).__new__(mcs, cls_or_name, bases, dct)
    # Now create the decorator
    def _decorate(cls):
        # Make sure we're decorating a class
        if not inspect.isclass(cls):
            raise TypeError("operable_class is a decorator for classes, not functions")
        # Strategy: make an individual metaclass for the class being created so we
        # can make 'class versions' of the requested operators
        # Also, we make the new metaclass inherit from the class's metaclass
        cls_meta = type.__new__(type, '_' + cls.__name__ + '___auto__meta_', (type(cls),), {'__new__': _meta_new})
        # And disguise the metaclass as it's parent...
        cls_meta.__name__ = type(cls).__name__
        # Now iterate over the operation names passed in
        for op in ops:
            cmeth_name = '__class_' + op + '__'
            # Only define the metaclass version if the class has a '__class_<op>__' function
            if hasattr(cls, cmeth_name):
                # Make sure the attribute is a classmethod
                attr = getattr(cls, cmeth_name)
                if getattr(attr, "im_self", None) is None:
                    raise TypeError("class magic method '{name}' must be defined as a classmethod".format(name=cmeth_name))
                # Make the metaclass wrapper
                def _cls_version(cls, *args, **kwargs):
                    attr = getattr(cls, kwargs.pop('_cls_version_method_call'))
                    return attr(*args, **kwargs)
                # and make it look like the class version...
                _cls_version.__name__ = '__' + op + '__'
                _cls_version.__doc__ = getattr(cls, cmeth_name).__doc__
                # now "freeze" the method name so we can redefine the wrapper for the next op
                part = bindablepartial(_cls_version, _cls_version_method_call=cmeth_name)
                # and set the attribute of the metaclass
                setattr(cls_meta, '__' + op + '__', part)
        return cls_meta(cls)
    return _decorate




