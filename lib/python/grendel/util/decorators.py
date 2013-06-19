""" Custom decorators for use in PyGrendel.
These are not really for the end-user.  If you are a developer, these will often save you a lot of time.
Note that if it makes sense for a decorator to be put in another module, it will be there rather than here.
"""
from __future__ import print_function
#TODO Add the decorator library (http://pypi.python.org/pypi/decorator), (from source) and make decorators preserve function signatures
from abc import ABCMeta
from collections import Iterable, Sequence, MutableSequence, Set, Mapping
from functools import wraps

#from new import instancemethod
from types import MethodType as instancemethod
import inspect
from inspect import getargspec
import sys
from grendel import caching_enabled, type_checking_enabled
import types
from grendel.util.abstract_bases import NonstringIterable, MethodLike
from grendel.util.aliasing import function_alias, bindablepartial
from grendel.util.exceptions import ProgrammerNeedsMoreCoffeeError
from grendel.util.strings import indented, make_safe_identifier, classname

# Nothing here needs to be imported when the end-user says 'from grendel import *'
__all__ = []

# Useful reference: http://wiki.python.org/moin/PythonDecoratorLibrary

#--------------------------------------------------------------------------------#
#                             Caching decorators                                 #
#--------------------------------------------------------------------------------#

#region Caching

# TODO @language-hack Take an option (to __init__, which will need to be changed to args/kwargs form) to cache only certain arg sets (if function takes arguments)
# TODO @language-hack Take an option (to __init__, which will need to be changed to args/kwargs form) that is a list of cached methods that depend on self so that calling reset on this calls reset on the other method
# TODO @language-hack option for automatic state checking and resetting


class CachedMethod(object):
    """
    Descriptor that caches return values of methods by their arguments.
    For now, if the method call depends on some mutable attribute of the
    class, do not use this.  Perhaps later a state-checking mechanism
    will be implemented that allows checking to see if dependent mutable
    attributes have changed.
    Note:  For now, functions with variable-length arguments and arbitrary
    keyword arguments are not supported
    """

    def __init__(self, func):
        self.func = func
        self.__doc__ = self.func.__doc__
        self.__name__ = self.func.__name__
        self.argspec = getargspec(self.func)

    def __get__(self, instance, owner):
        if instance is None:
            # Unbound method version
            return self
        else:
            # Bound method version
            if not hasattr(instance, CachedMethod.cache_name(self.__name__)):
                # Circumvent immutability checking, just for the cache
                if isinstance(owner, Immutable):
                    setattr_owner = object
                else:
                    setattr_owner = owner
                setattr_owner.__setattr__(instance, CachedMethod.cache_name(self.__name__), dict())
            return BoundCachedMethod(instance, self.func, owner)

    @staticmethod
    def cache_name(func_name):
        return "__" + func_name + "_method_call_cache__"

class BoundCachedMethod(object):

    def __init__(self, im_self, im_func, im_class):
        self.im_self = self.__self__ = im_self
        self.im_func = self.__func__ = im_func
        self.im_class = im_class
        self.__doc__ = self.__func__.__doc__
        self.__name__ = self.__func__.__name__
        self.argspec = getargspec(self.__func__)
        self.cache = getattr(self.im_self, CachedMethod.cache_name(self.__name__))

    def __call__(self, *args, **kwargs):
        args_needed = self.argspec.args[1:]
        ArgNotFilled = object()
        arg_key = [ArgNotFilled] * len(args_needed)
        defaults = self.argspec.defaults
        if defaults is None: defaults = []
        for iarg, arg in enumerate(args):
            try:
                arg_key[iarg] = arg
            except IndexError:
                if len(defaults) == 0:
                    raise TypeError("{0}() takes exactly {1} argument{2} ({3} given)".format(
                        self.__name__, len(args_needed)+1, 's' if len(args_needed) != 0 else '',
                        len(args) + len(kwargs) + 1
                    ))
                else:
                    raise TypeError("{0}() takes at most {1} argument{2} ({3} given)".format(
                        self.__name__, len(args_needed)+1, 's' if len(args_needed)-len(defaults) != 0 else '',
                        len(args) + len(kwargs) + 1
                    ))
        for key, value in kwargs.items():
            try:
                idx = args_needed.index(key)
            except ValueError:
                raise TypeError("{0}() got an unexpected keyword argument '{1}'".format(
                    self.__name__, key
                ))
            if arg_key[idx] is not ArgNotFilled:
                raise TypeError("{0}() got multiple values for keyword argument '{1}'".format(
                    self.__name__, key
                ))
            arg_key[idx] = value
        for idef, default in enumerate(defaults[::-1]):
            if arg_key[-idef-1] is ArgNotFilled:
                arg_key[-idef-1] = default
        if ArgNotFilled in arg_key:
            if len(defaults) == 0:
                raise TypeError("{0}() takes exactly {1} argument{2} ({3} given)".format(
                    self.__name__, len(args_needed)+1, 's' if len(args_needed) != 0 else '',
                    len(args) + len(kwargs) + 1
                ))
            else:
                raise TypeError("{0}() takes at least {1} argument{2} ({3} given)".format(
                    self.__name__, len(args_needed)-len(defaults)+1, 's' if len(args_needed)-len(defaults) != 0 else '',
                    len(args) + len(kwargs) + 1
                ))
        arg_key = tuple(arg_key)
        if arg_key in self.cache:
            return self.cache[arg_key]
        else:
            rv = self.__func__(self.im_self, *arg_key)
            self.cache[arg_key] = rv
            return rv

    def purge_cache(self):
        self.cache.clear()



def cached_method(func):
    """
    DEPRECATED
    Decorator for functions/methods with cached return values.
    For now, this can only be used.  The decorator creates a method that stores the return value of the method call
    in a variable named _[func_name] (where [func_name] is the name of the function.

    Examples
    --------

    >>> class Factorializer(object):
    ...
    ...    n = None
    ...
    ...    def __init__(self, n):
    ...        self.n = n
    ...
    ...    @cached_method
    ...    def compute(self):
    ...        if self.n == 0: return 1
    ...        ret_val = 1
    ...        for i in xrange(2, self.n+1):
    ...            ret_val *= i
    ...        return ret_val
    ...
    >>> f = Factorializer(3000)
    >>> # Takes some time...
    >>> var = f.compute()
    >>> # Should take no time...
    >>> for i in xrange(100):
    ...     var = f.compute()
    ...

    """
    def decorated(self):
        if not hasattr(self, "_" + func.__name__):
            self.__dict__["_" + func.__name__] = None
        if caching_enabled:
            if self.__dict__["_" + func.__name__] is not None:
                return self.__dict__["_" + func.__name__]
            self.__dict__["_" + func.__name__] = func(self)
            return self.__dict__["_" + func.__name__]
        else:
            return func(self)

    decorated.__name__ = func.__name__
    decorated.__doc__ = func.__doc__
    decorated.__dict__.update(func.__dict__)

    return decorated

# Simpler than CachedProperty, but does essentially the same thing...
class LazyProperty(object):

    def __init__(self, func):
        self._func = func
        self.__name__ = func.__name__
        self.__doc__ = func.__doc__

    def __get__(self, obj, klass=None):
        if obj is None: return None
        result = obj.__dict__[self.__name__] = self._func(obj)
        return result

class CachedProperty(object):
    """ Decorator for properties with cached return values.
    The decorator creates a property that stores the return value of the method call
    in a variable named _[func_name] (where [func_name] is the name of the function).

    Examples
    --------

    >>> class Factorializer(object):
    ...
    ...    n = None
    ...
    ...    def __init__(self, n):
    ...        self.n = n
    ...
    ...    @CachedProperty
    ...    def thirty_plus_n_fact(self):
    ...        ret_val = 1
    ...        for i in xrange(2, self.n + 31):
    ...            ret_val *= i
    ...        return ret_val
    ...
    >>> f = Factorializer(3000)
    >>> # Takes some time...
    >>> var = f.thirty_plus_n_fact
    >>> # Should take no time...
    >>> for i in xrange(100):
    ...     var = f.thirty_plus_n_fact
    ...
    >>> # Be careful... Remember to invalidate the cache when the appropriate variables are updated...
    >>> f.n = 2500
    >>> # Takes no time but gives the wrong answer:
    >>> var = f.thirty_plus_n_fact
    >>> # Invalidate the cache:
    >>> f._thirty_plus_n_fact = None
    >>> # This should take some time again, but now it will get the right answer
    >>> var = f.thirty_plus_n_fact



    """

    func = None
    name = None
    cached_var_name = None

    def __init__(self, func):
        """
        """
        self.func = func
        self.name = func.__name__
        self.cached_var_name = "_" + self.name
        self.__doc__ = func.__doc__
        self.__dict__.update(func.__dict__)

    def __get__(self, obj, obj_class=None):
        """
        """
        # First handle the case of call from class rather than instance
        if obj is None:
            return None
        # Now for instance call...
        if not hasattr(obj, self.cached_var_name):
            setattr(obj, self.cached_var_name, None)
        if caching_enabled:
            attr = getattr(obj, self.cached_var_name, None)
            if attr is not None:
                return attr
            else:
                ret_val = self.func(obj)
                setattr(obj, self.cached_var_name, ret_val)
                return ret_val
        else:
            return self.func(obj)


#endregion

#--------------------------------------------------------------------------------#
#                         Delegator defining decorators                          #
#--------------------------------------------------------------------------------#

#region | Delegator defining decorators |

#TODO explain delegation in documentation, make docstrings of delegates reference this explanation
#TODO check to make sure the user is not nesting delegation, which will break the decorator
def forwardable(cls):
    """ Class decorator that allows for the definition of delegators (in conjunction with `def_delegators`),
    just like the Forwardable module in Ruby.  *Note that this pattern does not work for nested classes*

    Examples
    --------
    >>> class NameContainer(object):
    ...     name = ""
    ...     def __init__(self, name):
    ...         self.name = name
    ...
    >>> @forwardable
    ... class AppendOnlyList(object):
    ...     list = []
    ...     def_delegators('list', 'append')
    ...     # Attributes can also be delegated:
    ...     name_container = None
    ...     def_delegators('name_container', 'name')
    ...     def __init__(self, name):
    ...         self.list = []
    ...         self.name_container = NameContainer(name)
    ...     # Double-underscore methods can't be delegated.  Instead, do something like this:
    ...     def __repr__(self):
    ...         return self.list.__repr__()
    ...
    >>> foo = AppendOnlyList("foobar")
    >>> foo.append('a')
    >>> foo.append(2)
    >>> foo.append('hello')
    >>> foo
    ['a', 2, 'hello']
    >>> foo.name
    'foobar'

    """
    if 'delegates' not in forwardable.__dict__:
        return cls
    for attr in forwardable.__dict__['delegates']:
        symbols = forwardable.__dict__['delegates'][attr]
        for symb in symbols:
            if hasattr(cls, symb):
                raise ValueError("Tried to define delegator of '{0}' to '{1}', but '{0}' is already defined in" \
                                 " class {2}.".format(symb, attr, classname(cls)))
            elif (symb[0], symb[1], symb[-2], symb[-1]) == ('_', '_', '_', '_'):
                raise ValueError("You are not allowed to delegate double-underscored methods.  Do these by hand, "
                                 "keeping in mind how Python actually accesses the particular double-underscored "
                                 "method you're trying to create.")
            else:
                class _wrapped_descriptor(object):
                    attr_name = None
                    symb_name = None
                    def __init__(self, in_attr, in_symb):
                        self.attr_name = in_attr
                        self.symb_name = in_symb
                    def __get__(self, obj, type=None):
                        return getattr(getattr(obj, self.attr_name), self.symb_name)
                _wrapped = _wrapped_descriptor(attr, symb)
                _wrapped.__name__ = attr
                _wrapped.__doc__ = """Delegated to `self.{0}`.""".format(attr)
                setattr(cls, symb, _wrapped)
    # Now clear the delegates for the next use.
    forwardable.__dict__['delegates'] = {}
    return cls

def def_delegators(attribute, *symbols_to_delegate):
    """ Define delegators for certain symbols in a class.
    """
    if 'delegates' not in forwardable.__dict__:
        forwardable.__dict__['delegates'] = {}
    delegate_dict = forwardable.__dict__['delegates']
    if attribute not in delegate_dict:
        delegate_dict[attribute] = []
    delegate_dict[attribute].extend(symbols_to_delegate)

#endregion

#--------------------------------------------------------------------------------#
#                       Assigning attributes to functions                        #
#--------------------------------------------------------------------------------#

#region | Assigning attributes to functions |

def with_attributes(**kwargs):
    """ Decorator allowing the assignment of attributes to a function.
    """
    def attribute_it(f):
        f.__dict__.update(kwargs)
        return f
    return attribute_it

#endregion

#--------------------------------------------------------------------------------#
#                             Flexible arguments                                 #
#--------------------------------------------------------------------------------#

#region | Flexible arguments |

def with_flexible_arguments(required=None, optional=None, what_to_call_it=None):
    """ Allows the creation of functions with case insensative, alternately named keyword arguments.

    Examples
    --------
    >>> @with_flexible_arguments(
    ...     required=[
    ...         ('name', 'greet', 'name_to_greet'),
    ...         ('from_where', 'source')
    ...     ],
    ...     optional=[
    ...         ('greeting', 'hello_phrase', 'salutation'),
    ...         ('message',)
    ...     ]
    ... )
    ... def say_hello(name, from_where, greeting='Hello', message=''):
    ...     print(greeting + ', ' + name + ', from ' + from_where + '!' + message)
    ...
    >>> say_hello('moon', 'Earth')
    Hello, moon, from Earth!
    >>> say_hello('world', source='the moon')
    Hello, world, from the moon!
    >>> say_hello(source='France', name_to_greet='visitor', salutation='Bonjour')
    Bonjour, visitor, from France!
    >>> say_hello('earthlings', 'outer space', 'Greetings', message='  We come in peace!')
    Greetings, earthlings, from outer space!  We come in peace!

    """
    required = required or []
    def wrap_it(f):
        if hasattr(f, 'getargspec'):
            fargspec = f.getargspec()
        else:
            fargspec = inspect.getargspec(f)
        fargs = fargspec.args
        allow_kwargs=False
        allow_args=False
        if fargspec.keywords is not None:
            allow_kwargs=True
        if fargspec.varargs is not None:
            allow_args=True
        ndefaults = len(fargspec.defaults) if fargspec.defaults else 0
        req_given = fargs[:-ndefaults] if len(fargs) > 0 else []
        opt_given = fargs[-ndefaults:] if len(fargs) > 0 and ndefaults > 0 else []
        # TODO Warning (for sanity check on) when fargs member is not in required/optional?
        # First, check for the presence of required arguments:
        for req_name in req_given:
            if optional is not None:
                bad = [o for o in optional if o[0] == req_name]
                if len(bad) != 0:
                    raise ValueError("Keyword '{0}' listed as optional, but given in argument"
                                     " specification for {1} without default value.".format(
                        bad[0],
                        what_to_call_it or (f.__name__ + '()'),
                    ))
            if not any(r[0] == req_name for r in required):
                required.append((req_name,))
        for opt_name in opt_given:
            bad = [r for r in required if r[0] == opt_name]
            if len(bad) != 0:
                print(fargspec)
                raise ValueError("Keyword '{0}' listed as required, but given in argument"
                                 " specification for {1} with default value '{2}'.".format(
                    bad[0],
                    what_to_call_it or (f.__name__ + '()'),
                    fargspec.defaults[opt_given.index(opt_name)]
                ))
            if not any(o[0] == opt_name for o in optional):
                optional.append((opt_name,))
        @wraps(f)
        def _wrapped(*args, **kwargs):
            new_kwargs = {}
            found_as_args = 0
            # First check for required arguments
            for req in required:
                found_arg = False
                # see if it can be found as a regular arg first...
                if req[0] in fargs and fargs.index(req[0]) < len(args):
                    found_as_args += 1
                    found_arg = True
                for possible in req:
                    for kwarg in kwargs.keys():
                        if possible.lower() == kwarg.lower():
                            val = kwargs.pop(kwarg)
                            if req[0] in new_kwargs or found_arg:
                                raise TypeError("{0} got multiple values for keyword argument '{1}'{2}.".format(
                                    what_to_call_it or (f.__name__ + '()'),
                                    req[0],
                                    " (which is the same thing as {0})".format(possible) if possible != req[0] else ''
                                ))
                            else:
                                new_kwargs[req[0]] = val
                if not found_arg and not req[0] in new_kwargs:
                    # The second part is in case it was found as an arg rather than a kwarg
                    raise TypeError("{2} missing required keyword argument '{0}'{1}.  "
                                    "Check your spelling and try again.".format(
                                        req[0],
                                        " or one of its alternate keywords: ('{0}')".format("', '".join(req[1:])) if len(req) > 1 else '',
                                        what_to_call_it or (f.__name__ + '()')
                                    ))
            #--------------------------------------------------------------------------------#
            # Now check for optional arguments and substitute alternate names for standard ones:
            if optional is not None:
                for opt in optional:
                    found_arg = found = False
                    # see if it can be found as a regular arg first...
                    if opt[0] in fargs and fargs.index(opt[0]) < len(args):
                        found_as_args += 1
                        found_arg = True
                    for possible in opt:
                        for kwarg in kwargs.keys():
                            if possible.lower() == kwarg.lower():
                                val = kwargs.pop(kwarg)
                                if opt[0] in new_kwargs or found_arg:
                                    raise TypeError("{0} got multiple values for keyword argument '{1}'{2}.".format(
                                        what_to_call_it or (f.__name__ + '()'),
                                        opt[0],
                                        " (which is the same thing as {0})".format(possible) if possible != opt[0] else ''
                                    ))
                                else:
                                    new_kwargs[opt[0]] = val
            #--------------------------------------------------------------------------------#
            # Now check the args and raise errors if necessary
            if allow_kwargs:
                new_kwargs.update(kwargs)
            elif len(kwargs) != 0:
                raise TypeError("{0} got unexpected keyword argument '{1}'".format(
                    what_to_call_it or (f.__name__ + '()'),
                    kwargs.keys()[0]
                ))
            num_not_included = len([a for a in fargs if not any(p[0] == a for p in (required or []) + (optional or []))])
            if len(args) != found_as_args + num_not_included and not allow_args:
                if len(fargs) == 0:
                    raise TypeError("{name} takes no arguments ({num} given)".format(
                        name = what_to_call_it or (f.__name__ + '()'),
                        num=len(args)
                    ))
                else:
                    num_given = len([a for a in kwargs if a in fargs]) + len(args)
                    raise TypeError("{name} takes exactly {num} argument{plural} ({given} given)".format(
                        name=what_to_call_it or (f.__name__ + '()'),
                        num=len(fargs),
                        plural='s' if len(fargs) != 1 else '',
                        given=num_given
                    ))
            #--------------------------------------------------------------------------------#
            return f(*args, **new_kwargs)
        return _wrapped
    return wrap_it

#endregion

#--------------------------------------------------------------------------------#
#                                 Aliases                                        #
#--------------------------------------------------------------------------------#

#region | Aliases |

# DOESN'T WORK
#class MethodAlias(object):
#    """ A method that is an alias for another method.
#    """
#
#    func = None
#
#    def __init__(self, func):
#        self.__doc__ = "Alias for `" + func.__name__ + "`"
#        self.func = func
#
#    # TODO Modify stacktrace
#    def __call__(self, *args, **kwargs):
#        self.func(*args, **kwargs)


# TODO get something like this working...
#def alias_for(func, argmap=None):
#    """ Returns an alias for the method `f`
#
#    Examples
#    --------
#
#    >>> class A(object):
#    ...     def foo(self, arg):
#    ...         \""" Says Hello \"""
#    ...         print("Hello, " + arg)
#    ...     # the contents of bar are ignored
#    ...     @alias_for(foo)
#    ...     def bar(self): pass
#    ...
#    >>> A.foo("world")
#    Hello, world
#    >>> A.bar("moon")
#    Hello, moon
#    >>> A.foo.__doc__
#    'Says Hello'
#    >>> A.bar.__doc__
#    'Alias for `foo`'
#
#    """
#    def _decorator(f):
#        # TODO Raise an error unless the contents of f are "pass"
#        if argmap is None:
#            @wraps(f)
#            def wrapped(*args, **kwargs):
#                return func(*args, **kwargs)
#            wrapped.__doc__ = "Alias for `" + func.__name__ + "`"
#        else:
#            # TODO Finish this
#            raise NotImplementedError
#            #fargs = inspect.getargspec(f)
#            #funcargs = inspect.getargspec(func)
#            #rev_argmap = {}
#            #for key in argmap:
#            #    rev_argmap =
#            #@wraps(f)
#            #def wrapped(*args, **kwargs):
#            #    # Build the kwargs to pass to func
#            #    pass_kwargs = {}
#            #    for n, arg in enumerate(args):
#            #        pass_kwargs[argmap]
#        return wrapped
#    return _decorator

#endregion

#--------------------------------------------------------------------------------#
#                               Delayed Decoration                               #
#--------------------------------------------------------------------------------#

#region | Delayed Decoration |

# Source: http://code.activestate.com/recipes/577993/

class DelayedDecorator(object):
    """Wrapper that delays decorating a function until it is invoked.

    This class allows a decorator to be used with both ordinary functions and
    methods of classes. It wraps the function passed to it with the decorator
    passed to it, but with some special handling:

      - If the wrapped function is an ordinary function, it will be decorated
        the first time it is called.

      - If the wrapped function is a method of a class, it will be decorated
        separately the first time it is called on each instance of the class.
        It will also be decorated separately the first time it is called as
        an unbound method of the class itself (though this use case should
        be rare).
    """

    def __init__(self, deco, func):
        # The base decorated function (which may be modified, see below)
        self._func = func
        # The decorator that will be applied
        self._deco = deco
        # FiniteDifferenceVariable to monitor calling as an ordinary function
        self.__decofunc = None
        # FiniteDifferenceVariable to monitor calling as an unbound method
        self.__clsfunc = None

    def _decorated(self, cls=None, instance=None):
        """Return the decorated function.

        This method is for internal use only; it can be implemented by
        subclasses to modify the actual decorated function before it is
        returned. The ``cls`` and ``instance`` parameters are supplied so
        this method can tell how it was invoked. If it is not overridden,
        the base function stored when this class was instantiated will
        be decorated by the decorator passed when this class was instantiated,
        and then returned.

        Note that factoring out this method, in addition to allowing
        subclasses to modify the decorated function, ensures that the
        right thing is done automatically when the decorated function
        itself is a higher-order function (e.g., a generator function).
        Since this method is called every time the decorated function
        is accessed, a new instance of whatever it returns will be
        created (e.g., a new generator will be realized), which is
        exactly the expected semantics.
        """
        return self._deco(self._func)

    def __call__(self, *args, **kwargs):
        """Direct function call syntax support.

        This makes an instance of this class work just like the underlying
        decorated function when called directly as an ordinary function.
        An internal reference to the decorated function is stored so that
        future direct calls will get the stored function.
        """
        if not self.__decofunc:
            self.__decofunc = self._decorated()
        return self.__decofunc(*args, **kwargs)

    def __get__(self, instance, cls):
        """Descriptor protocol support.

        This makes an instance of this class function correctly when it
        is used to decorate a method on a user-defined class. If called
        as a bound method, we store the decorated function in the instance
        dictionary, so we will not be called again for that instance. If
        called as an unbound method, we store a reference to the decorated
        function internally and use it on future unbound method calls.
        """
        if instance:
            deco = instancemethod(self._decorated(cls, instance), instance, cls)
            # This prevents us from being called again for this instance
            setattr(instance, self._func.__name__, deco)
        elif cls:
            if not self.__clsfunc:
                self.__clsfunc = instancemethod(self._decorated(cls), None, cls)
            deco = self.__clsfunc
        else:
            raise ValueError("Must supply instance or class to descriptor.")
        return deco

#endregion

#--------------------------------------------------------------------------------#
#                            Method Decorators                                   #
#--------------------------------------------------------------------------------#

#region | Method Decorators |
see_abstract_doc = with_attributes(__needs_abstract_doc__=True)
see_abstract_doc.__doc__ = """
Add to the documentation of the method (or create it, if it is empty) a reference
to the documentation of the first abstractmethod in the containing class's __mro__.

Examples
--------
>>> from abc import ABCMeta, abstractmethod
>>> class A(object):
...     __metaclass__ = ABCMeta
...     @abstractmethod
...     def foo(self):
...         ''' Useful description of what foo() should do in subclasses '''
...         return NotImplemented
...
>>> @uses_method_decorators
... class B(A):
...     @see_abstract_doc
...     def foo(self):
...         print("hello world")
...
>>> b = B()
>>> b.foo()
hello world
>>> print(B.foo.__doc__)
Concrete implementation of abstract method `A.foo()`

"""

# TODO doc, test
# deprecated.  Use aliased_function instead
#def aliased_method(*names):
#    """
#    Alias a method as multiple names in a class.
#
#
#    """
#    return with_attributes(__aliases__=listify_args(names))
#
#method_aliased_as = function_alias('method_aliased_as', aliased_method)

# TODO property alias as a simple subclass of property with one extra attribute (then implement in uses_method_decorators)


def uses_method_decorators(cls):
    """
    Some decorators for methods such as `alias` and `see_abstract_doc`
    can't be complete without access to the containing class object.
    This decorator goes through and finishes the job.  To see what
    this does for individual 'method decorators', see the documentation
    for the given decorator.
    """
    if not inspect.isclass(cls):
        raise SyntaxError("the '@uses_method_decorators' decorator is only for decorating classes")
    for name, val in cls.__dict__.items():
        #--------------------------------------------------------------------------------#
        if hasattr(val, '__needs_abstract_doc__'):
            #--------------------------------#
            # Completion of see_abstract_doc #
            #--------------------------------#
            abs_meths = []
            for base in cls.__mro__[1:]:
                if name in base.__dict__:
                    meth = getattr(base, name, None)
                    if getattr(meth, '__isabstractmethod__', False):
                        abs_meths.append((meth, base))
            if len(abs_meths) == 0:
                raise SyntaxError("see_abstract_doc decorator cannot be completed:"
                                  "  '{0}' is not an implementation of an abstract method in"
                                  " any of the base classes of '{1}'".format(name, cls.__name__))
            elif len(abs_meths) == 1:
                doc_add = "Concrete implementation of abstract method `{0}.{1}()`".format(abs_meths[0][1].__name__, abs_meths[0][0].__name__)
            else:
                doc_add = "\nConcrete implementation of abstract methods (`{0}` and `{1}`)".format(
                    '`, `'.join("{0}.{1}()".format(base.__name__, meth.__name__) for meth, base in abs_meths[:-1]),
                    "{0}.{1}()".format(abs_meths[-1][1].__name__, abs_meths[-1][0].__name__)
                )
            if val.__doc__:
                val.__doc__ += "\n" + doc_add
            else:
                val.__doc__ = doc_add
        #--------------------------------------------------------------------------------#
        if hasattr(val, '__aliases__'):
            # DEPRECATED
            #------------------------------#
            # Completion of aliased_method #
            #------------------------------#
            # TODO Signature-preserving wrapping
            for alias in val.__aliases__:
                if not isinstance(alias, str):
                    raise TypeError('aliases must be strings')
                if alias in cls.__dict__:
                    raise NameError("cannot alias '{0}' as '{1}' since the name '{1}' is already defined in '{2}'".format(
                        val.__name__,
                        alias,
                        cls.__name__
                    ))
                def _aliased_func(_method_decorator_alias_for, *args, **kwargs):
                    return _method_decorator_alias_for(*args, **kwargs)
                aliaspart = bindablepartial(_aliased_func, _method_decorator_alias_for=val)
                aliaspart.__name__ = alias
                aliaspart.__doc__ = "Alias for `{0}.{1}()`".format(
                    cls.__name__,
                    val.__name__
                )
                setattr(cls, alias, aliaspart)
        #--------------------------------------------------------------------------------#
    return cls

#endregion

#--------------------------------------------------------------------------------#
#                               Type Checking                                    #
#--------------------------------------------------------------------------------#

#region | Type Checking |
# helper class, not a decorator
#TODO @language-hack has_attribute-like "type" recognition (doable already with lambdas, but it would be nice to give a readable output)
#TODO @language-hack SubclassOf (or TypeIsSubclassOf?) (or something more intuitive that will not be confusing (!))
#TODO @language-hack DictOf
class argtypespec(object):
    """ Contains the type specifications for a set of arguments to a method.
    Utility class for use with `overloaded`.  This used to be an inner class, but
    I think I might use it for type checking in the future.
    """

    ##############
    # Attributes #
    ##############

    func = None
    argspec = None
    argtypes = None

    ##################
    # Initialization #
    ##################

    def __init__(self, func, typespec_dict):
        self.func = func
        self.argspec = inspect.getargspec(func)
        if len(typespec_dict) != len(self.argspec.args):
            raise ValueError("dimension mismatch ({} != {})".format(len(typespec_dict), len(argspec.args)))
            # any other type checking and sanity checking should be done at this point
        self.argtypes = {}
        self.argtypes = {}
        for arg, ty in typespec_dict.items():
            if isinstance(ty, NonstringIterable):
                self.argtypes[arg] = ty
            else:
                self.argtypes[arg] = (ty,)

    ##############
    # Properties #
    ##############

    @property
    def required(self):
        if self.argspec.defaults:
            return self.argspec.args[:-len(self.defaults)]
        else:
            return self.argspec.args

    @property
    def optional(self):
        if self.argspec.defaults:
            return self.argspec.args[-len(self.defaults):]
        else:
            return []

    @property
    def defaults(self):
        return self.argspec.defaults or tuple([])

    @property
    def signature(self):
        # TODO multiline output if necessary
        ret_val = self.func.__name__ + "("
        for name in self.required:
            types = self.argtypes[name]
            ret_val += name
            if types == (None,) or types is AnyType or AnyType in types:
                ret_val += ": <any type>"
            elif len(types) == 1:
                ret_val += ": " + self._arg_type_str(name, types[0])
            elif len(types) > 1:
                ret_val += ": ({})".format(', '.join(self._arg_type_str(name, tyty) for tyty in types))
            if name != self.required[-1]:
                ret_val += ', '
        for name, default in zip(self.optional, self.defaults):
            types = self.argtypes[name]
            ret_val += ', '
            ret_val += name
            if types == (None,) or types is AnyType or AnyType in types:
                ret_val += ": <any type>"
            elif len(types) == 1:
                ret_val += ": " + self._arg_type_str(name, types[0])
            elif len(types) > 1:
                ret_val += ": ({})".format(', '.join(self._arg_type_str(name, tyty) for tyty in types))
            ret_val += '=' + repr(default)
        if self.argspec.varargs is not None:
            if len(self.required) + len(self.optional) > 0:
                ret_val += ', '
            ret_val += "*" + self.argspec.varargs
        if self.argspec.keywords is not None:
            if len(self.required) + len(self.optional) > 0 or self.argspec.varargs is not None:
                ret_val += ', '
            ret_val += "**" + self.argspec.keywords
        ret_val += ')'
        return ret_val

    ###################
    # Special Methods #
    ###################

    @staticmethod
    def arg_is_instance_of(arg, ty, globals=None):
        # If it's a string, evaluate it in the context of the function's globals (with no locals)
        # and use the result as the type instead
        from_string = False
        if isinstance(ty, basestring):
            from_string = True
            ty = eval(ty, globals, {})
        # Now check the possibilities...
        if isinstance(ty, type):
            if isinstance(arg, ty):
                # It's the correct type
                return True
        elif isinstance(ty, IterableOf):
            if not ty.arg_type_okay(arg):
                return False
            else:
                return all(argtypespec.arg_is_instance_of(a, ty.types, globals) for a in arg)
        elif callable(ty):
            if ty(arg):
                # It responded True to the call, so assume it's the right type
                return True
        elif ty is None:
            if arg is None and not from_string:
                # We'll let None stand in for NoneType, since NoneType isn't all that understandable
                # unless it's the result of a string evaluation...
                return True
        return False

    @staticmethod
    def get_call_signature(_get_call_signature_name, *args, **kwargs):
        iters_to_expand=(list, tuple, set,)
        ret_val = _get_call_signature_name + "("
        argtypes = []
        for a in args:
            if isinstance(a, iters_to_expand):
                argtypes.append(argtypespec.iterable_of_call_signature(a))
            else:
                argtypes.append(type(a).__name__)
        ret_val += ', '.join(argtypes)

        if len(kwargs) > 0:
            ret_val += ', '
            kwargexprs = []
            for k, v in kwargs.items():
                if isinstance(v, iters_to_expand):
                    kwargexprs.append('{}={}'.format(
                        str(k),
                        argtypespec.iterable_of_call_signature(v)
                    ))
                else:
                    kwargexprs.append('{}={}'.format(str(k), type(v).__name__))
            ret_val += ', '.join(kwargexprs)

        ret_val += ')'
        return ret_val

    @staticmethod
    def iterable_of_call_signature(iterable, iters_to_expand=(list, tuple, set,)):
        typenamelist = set()
        for item in iterable:
            # TODO @fix @edge-case Make sure to avoid infinite recursion if iter(item).next() returns something that is the same as item
            if isinstance(item, iters_to_expand):
                typenamelist.add("<Iterable of {}>".format(argtypespec.iterable_of_call_signature(item)))
            else:
                typenamelist.add(type(item).__name__)
        if len(typenamelist) == 1:
            return '<{} of {}>'.format(type(iterable).__name__, typenamelist.pop())
        else:
            return '<{} of ({})>'.format(type(iterable).__name__, ', '.join(typenamelist))


    ###########
    # Methods #
    ###########

    def get_arg_list(self, *args, **kwargs):
        names = self.argspec.args
        specified = names[:len(args)]
        required = self.required
        unspecified = [r for r in required if r not in specified]
        #========================================#
        arg_dict = {}
        args = list(args)
        for name in specified if len(specified) < len(required) else required:
            if name in kwargs:
                # the argument is specified multiple times...not a valid call
                return None, "got multiple values for keyword argument '{}'".format(name)
            arg_dict[name] = args.pop(0)
        for name in unspecified:
            if name not in kwargs:
                # missing required argument...not a valid call
                return None, "missing required argument '{}'".format(name)
            arg_dict[name] = kwargs.pop(name)
        #========================================#
        # optional arguments
        optional = self.optional
        defaults_dict = dict(zip(optional, self.defaults))
        for name in optional:
            if len(args) > 0:
                if name in kwargs:
                    # the argument is specified multiple times...not a valid call
                    return None, "argument '{}' is specified multiple times".format(name)
                arg_dict[name] = args.pop(0)
            else:
                arg_dict[name] = kwargs.pop(name, defaults_dict[name])
        #========================================#
        # Check the types
        if len(args) > 0 and self.argspec.varargs is None:
            return None, None
        for name in names:
            types = self.argtypes[name]
            arg = arg_dict[name]
            if not (types == (None,) or types is AnyType or AnyType in types):
                type_found = False
                for ty in types:
                    if argtypespec.arg_is_instance_of(arg, ty, self.func.func_globals):
                        type_found = True
                        break
                if not type_found:
                    return None, "argument '{}' is not of the correct type".format(name)
        #========================================#
        # build the argument list and return
        ret_val = [arg_dict[name] for name in names]
        if self.argspec.varargs is None:
            if len(args) > 0:
                # remaining args, but no varargs variable in function spec
                return None, "additional arguments given, but no variable arguments allowed"
            if self.argspec.keywords is None:
                if len(kwargs) > 0:
                    # remaining kwargs, but no keyword arguments variable in function spec
                    return None, "additional keyword arguments given, but no keyword arguments allowed"
                return ret_val, {}
            else:
                return ret_val, kwargs
        elif self.argspec.keywords is None:
            if len(kwargs) > 0:
                # remaining kwargs, but no keyword arguments variable in function spec
                return None, "additional keyword arguments given, but no keyword arguments allowed"
            return ret_val + args, {}
        else:
            return ret_val + args, kwargs


    ###################
    # Private Methods #
    ###################

    def _arg_type_str(self, arg, tyty):
        if isinstance(tyty, type):
            return tyty.__name__
        elif callable(tyty):
            return '<something for which {0}({1}) returns True>'.format(tyty.__name__, arg)
        elif isinstance(tyty, basestring):
            return tyty
        elif isinstance(tyty, IterableOf):
            if isinstance(tyty, types.GeneratorType):
                raise NotImplementedError('generator type-checking not implemented.  Just make it into a list for now')
            return '<{} of {}>'.format(
                tyty.collection_type.__name__,
                '({})'.format(
                    self._arg_type_str('<item>', t) for t in tyty.types) \
                        if isinstance(tyty.types, NonstringIterable) else self._arg_type_str('<item>', tyty.types)
            )
        elif tyty is AnyType:
            return '<any type>'
        elif tyty is None:
            return 'NoneType'
        else: # pragma: no cover
            raise ProgrammerNeedsMoreCoffeeError("Type specification passed to argtypespec" \
                                                 " was '{}' of type '{}', which doesn't" \
                                                 " make sense.".format(
                tyty,
                type(tyty).__name__
            ))

AnyType = object()

#TODO @URGENT serious bug: generators passed in don't get reset
class IterableOf(object):
    collection_type = Iterable
    types = None

    def __init__(self, types):
        self.types = types

    def arg_type_okay(self, arg):
        return isinstance(arg, self.collection_type)

class SequenceOf(IterableOf):
    collection_type = Sequence

class MutableSequenceOf(IterableOf):
    collection_type = MutableSequence

class SetOf(IterableOf):
    collection_type = Set

class MappingOf(IterableOf):
    collection_type = Mapping

class DictOf(IterableOf):
    collection_type = dict

class ListOf(IterableOf):
    collection_type = list

class TupleOf(IterableOf):
    collection_type = tuple

if type_checking_enabled:
    # TODO documentation and testing
    class typechecked_function(MethodLike):
        """
        Type checking for functions and methods.

        """

        ######################
        # Private Attributes #
        ######################

        _typespec = None

        ##################
        # Initialization #
        ##################

        def __new__(cls, *typespecs_or_func, **kwargs):
            self = object.__new__(cls)
            # TODO check to make sure type specifications are valid
            if sys.version_info > (3, 0)\
                    and len(kwargs) == 0\
                    and len(typespecs_or_func) == 1\
                    and inspect.isfunction(typespecs_or_func[0])\
                    and not inspect.isbuiltin(typespecs_or_func[0]):
                func = typespecs_or_func[0]
                if hasattr(func, 'func_annotations'):
                    # Python 3 version with func_annotations:
                    argspec = inspect.getargspec(func)
                    typespec_dict = {}
                    for arg in argspec:
                        typespec_dict[arg] = func.func_annotations[arg]
                    self._typespec = argtypespec(func, typespec_dict)
                    self.__name__ = func.__name__
                    self.__doc__ = func.__doc__
                    self.__dict__.update(func.__dict__)
                    return self
                else:
                    raise TypeError("when not using Python 3 function annotations, "
                                    "{}.overload_with must\nbe called with "
                                    "type specifications as arguments".format(self.__name__))
            else:
                def _get_func_decorator(f):
                    argspec = inspect.getargspec(f)
                    typespec_dict = {}
                    # remember zip always truncates at the length of the shortest iterator
                    for spec, argname in zip(typespecs_or_func, argspec.args):
                        typespec_dict[argname] = spec
                    for argname, spec in kwargs.items():
                        if argname in typespec_dict:
                            raise TypeError("multiple type specifications given for"
                                            " argument {} to typechecked function {}()".format(argname, f.__name__))
                        if argname not in argspec.args:
                            raise TypeError("type specifiation for unknown argument {} given to typechecked"
                                            " function {}()".format(argname, f.__name__))
                        typespec_dict[argname] = spec
                        # fill in None for anything missing, indicating any type will do
                    for argname in argspec.args:
                        if argname not in typespec_dict:
                            typespec_dict[argname] = None
                    self._typespec = argtypespec(f, typespec_dict)
                    self.__name__ = f.__name__
                    self.__doc__ = f.__doc__
                    self.__dict__.update(f.__dict__)
                    return self
                return _get_func_decorator

        ###################
        # Special Methods #
        ###################

        def __get__(self, instance, owner):
            # Mimic the behavior of the built-in function type
            return types.MethodType(self, instance, owner)

        def __call__(self, *args, **kwargs):
            fargs, fkwargs = self._typespec.get_arg_list(*args, **kwargs)
            if fargs is not None:
                return self._typespec.func(*fargs, **fkwargs)
            else:
                raise TypeError("invalid function signature {0}.\n  Function call must have" \
                                " the form:\n{1}\n    which is not a valid signature because: {2}".format(
                    argtypespec.get_call_signature(self.__name__, *args, **kwargs),
                    indented(self._typespec.signature),
                    fkwargs
                ))

        ###########
        # Methods #
        ###########

        def getargspec(self):
            """ Pass-through to function to conform with the `FunctionLike` protocol.
            """
            return inspect.getargspec(self._typespec.func)
else:
    # Typechecking is off, try to streamline things as much as possible
    def typechecked_function(*args, **kwargs):
        def _passthrough(f):
            return f
        return _passthrough
# aliases
typechecked = typechecked_function
typechecked_method = typechecked_function
typecheck = typechecked_function

#endregion

#####################
# Dependent Imports #
#####################

from grendel.util.metaclasses import Immutable

###########
# Aliases #
###########

#cached_method = function_alias('cached_method', CachedMethod)


