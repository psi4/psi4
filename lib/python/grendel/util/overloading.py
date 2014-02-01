""" Utilities for making it easier to have functions be callable in multiple forms.
"""
from collections import Iterable
import sys
import inspect
import re
import types
from grendel.util.abstract_bases import MethodLike, NonstringIterable
from grendel.util.aliasing import function_alias, bindablepartial
from grendel.util.decorators import with_flexible_arguments
from grendel.util.strings import indented

__all__ = [
    "caller",
    "listify_args",
    "get_multikwarg", "get_kwarg",
    "pop_multikwarg", "pop_kwarg",
]

#TODO move full documentation to a separate overloading.rst file

class OverloadedFunctionCallError(TypeError):
    """ Raised when something goes wrong with an overloaded function.
    Right now, this is mostly just a placeholder for the kind of thing that could be raised by the
    default implementation of an overloaded function (i.e. the contents of the argument to the overloaded constructor)
    if the function is not supposed to support a default implementation.  This functionality is not yet implemented.
    However, for future-proofness, all "default implementations" should raise an `OverloadedFunctionCallError` for
    now.
    """

def listify_args(*args, **kwargs):
    """ A simple way to accept arguments that can be an Iterable or a bunch of arguments that should be treated as a single Iterable.
    (This method is pretty trivial, but I find myself doing it a whole lot.)  The reason you can't just use `list()`
    for this purpose is that if an Iterable is given as the only argument, `list(args)` will return a list of Iterables
    rather than a single Iterable.

    :Examples:

    >>> def test(*args):
    ...     return listify_args(*args)
    ...
    >>> test(1,2,3)
    [1, 2, 3]
    >>> test([1,2,3])
    [1, 2, 3]
    >>> test([1,2,3], 4)
    [[1, 2, 3], 4]
    >>> test("1,2,3", 4)
    ['1,2,3', 4]
    >>> test("1,2,3")
    ['1,2,3']

    """
    ignore = kwargs.pop('ignore', basestring)
    if len(args) == 1 and isinstance(args[0], Iterable) and not isinstance(args[0], ignore):
        return list(args[0])
    else:
        return list(args)


def tuplify_args(*args, **kwargs):
    return tuple(listify_args(*args, **kwargs))

def get_multikwarg(kwarg_dict, *args):
    """ Utility function for getting the value of a keyword argument that can be named multiple things.
    This function also checks to make sure no more than one of the possible names is specified in the
    keyword dictionary `kwarg_dict`.  If the None of the keywords given are found in `kwarg_dict`, `None`
    is returned.

    """
    names = listify_args(*args)
    ret_val = None
    found = None
    for name in names:
        if name in kwarg_dict:
            if found:
                raise TypeError("Multiple keywords meaning the same thing given in call to {0}():"
                                " You may specify '{1}' or '{2}', but not both.".format(caller(), found, name))
            else:
                ret_val = kwarg_dict[name]
                found = name
    return ret_val
get_kwarg = get_multikwarg

# Deprecated...ish
def pop_multikwarg(kwarg_dict, *args):
    """ Utility function for getting the value of a keyword argument that can be named multiple things and popping
    that value off of the kwarg dictionary.  If the None of the keywords given are found in `kwarg_dict`, `None` is
     returned.

     :See Also:

     get_multikwarg

    """
    names = listify_args(*args)
    ret_val = None
    found = None
    for name in names:
        if name in kwarg_dict:
            if found:
                raise TypeError("Multiple keywords meaning the same thing given in call to {0}():"
                                " You may specify '{1}' or '{2}', but not both.".format(caller(), found, name))
            else:
                ret_val = kwarg_dict.pop(name)
                found = name
    return ret_val
pop_kwarg = pop_multikwarg

# TODO move this to a more logical place?
def caller():
    """ Get the name of the calling function as a `str`

    :Returns:

    out : str
        The name of the function that called the function whose context `caller()` is called from.  (Hopefully the
        examples makes this clearer.)

    :Examples:


    >>> def foo():
    ...     print caller()
    ...
    >>> def bar():
    ...     foo()
    ...
    >>> bar()
    bar
    >>> def foobar():
    ...     bar()
    ...
    >>> foobar()
    bar

    """
    frame = inspect.currentframe()
    try:
        return inspect.getouterframes(frame)[2][3]
    finally:
        del frame


#TODO @language-hack common-to-all optional arguments
#TODO @documentation modify (or replace?) numpydoc style to accomodate overloaded functions cleanly
#TODO @CRITICAL calls to __init__ in overloaded function go to subclass init if it is overloaded
class overloaded(MethodLike):
    """ Overload a function the Pythonic (at least, Python 2) way.

    Types can be specified using the `overloaded` instance decorator method `overload_with` (which is aliased as
    `overload`, `submethod`, and `subfunction`).

    :Examples:

    >>> @overloaded
    ... def test():
    ...     raise TypeError
    ...
    >>> @test.overload_with(int, int)
    ... def test(a, b=5):
    ...     return a+b
    ...
    >>> @test.overload_with(str)
    ... def test(string, *args, **kwargs):
    ...     return string + (' ' if len(args) else '') + ', '.join(list(args) + kwargs.keys())
    ...
    >>> @test.overload_with(list, func=callable)
    ... def test(lst, func, *args):
    ...     for item in lst:
    ...         if func(item, *args):
    ...             return item
    ...
    >>> test(1, 3)
    4
    >>> test(17)
    22
    >>> test("Hello")
    'Hello'
    >>> test("Hello", "World")
    'Hello World'
    >>> test("Hello", "Earth", "Mars", "Jupiter")
    'Hello Earth, Mars, Jupiter'
    >>> test([1, 2, 3], lambda x: x**2 % 2 == 0)
    2
    >>> test([1, 2, 3], lambda x, y: x**2 % y == 0, 1)
    1
    >>> test(3.14159, 3.234)
    Traceback (most recent call last):
        ...
    TypeError: invalid function signature test(float, float).
      Available signatures are:
        test(a: int, b: int=5)
            Not valid signature because: argument 'a' is not of the correct type
        test(string: str, *args, **kwargs)
            Not valid signature because: argument 'string' is not of the correct type
        test(lst: list, func: <something for which callable(func) returns True>, *args)
            Not valid signature because: argument 'lst' is not of the correct type
    >>> # giving multiple values for a keyword argument means there is no match.  For instance:
    >>> test('Hello', 'World', string='this is wrong')
    Traceback (most recent call last):
        ...
    TypeError: invalid function signature test(str, str, string=str).
      Available signatures are:
        test(a: int, b: int=5)
            Not valid signature because: argument 'a' is not of the correct type
        test(string: str, *args, **kwargs)
            Not valid signature because: got multiple values for keyword argument 'string'
        test(lst: list, func: <something for which callable(func) returns True>, *args)
            Not valid signature because: argument 'lst' is not of the correct type
    >>> # Methods in classes
    >>> class Foo(object):
    ...     @overloaded
    ...     def bar(self, *args, **kwargs): pass
    ...
    ...     @bar.overload_with(a=int, b=int)
    ...     def bar(self, a, b):
    ...         if not hasattr(self, 'total'):
    ...             self.total = 0
    ...         self.total += a + b
    ...         return a + b
    ...
    >>> foo = Foo()
    >>> foo.bar(1,2)
    3
    >>> foo.bar(4,5)
    9
    >>> foo.total
    12
    >>> Foo.bar(foo, 7, 8)
    15
    >>> foo.total
    27

    """

    ######################
    # Private Attributes #
    ######################

    _orig_func = None
    _versions = None

    ##################
    # Initialization #
    ##################

    def __init__(self, f):
        self._orig_func = f
        self.__name__ = f.__name__
        self.__doc__ = f.__doc__
        self.__dict__.update(f.__dict__)
        self._versions = []

    ###################
    # Special Methods #
    ###################

    def __get__(self, instance, owner):
        # Mimic the behavior of the built-in function type
        return types.MethodType(self, instance, owner)

    def __call__(self, *args, **kwargs):
        fail_reasons = []
        for typespec in self._versions:
            fargs, fkwargs_or_fail_reason = typespec.get_arg_list(*args, **kwargs)
            if fargs is not None:
                return typespec.func(*fargs, **fkwargs_or_fail_reason)
            else:
                fail_reasons.append(fkwargs_or_fail_reason or '<unknown reason for failed match>')
        # Should we fall back on the original function's code here?  If so, we should catch and reraise any type errors
        raise TypeError("invalid function signature {0}.\n  Available signatures are:\n{1}".format(
            argtypespec.get_call_signature(self.__name__, *args, **kwargs),
            indented('\n'.join(line + "\n" + indented("Not valid signature because: " + reason)
                        for line, reason in zip(self._build_allowed_string().splitlines(), fail_reasons)))

        ))

    ###########
    # Methods #
    ###########

    def getargspec(self):
        """ Pass-through to original function to conform with the `FunctionLike` protocol
        """
        return inspect.getargspec(self._orig_func)

    def overload_with(self, *typespecs_or_func, **kwargs):
        """ Add a overloaded version of the decorated function to the available options for a given overloaded function.
        Types can be specified several ways.  The simplest is to give the types as plain old arguments to
        the `overload_with` decorator.  The order of the arguments given corresponds to the order of the arguments in
        the function definition.  For example:

            >>> @overloaded
            ... def foo():
            ...     pass
            ...
            >>> @foo.overload_with(int, str)
            ... def foo(num, word):
            ...     return word * num
            ...
            >>> foo(3, 'hello')
            'hellohellohello'
            >>> foo('hello', 3)  # doctest: +ELLIPSIS
            Traceback (most recent call last):
                ...
            TypeError: invalid function signature ...

        If fewer types are specified than arguments, it is assumed that each of the remaining arguments can be any type.
        (specifying `None` for the type has the same effect).  Keyword arguments may also be given to `overload_with`,
        where the keyword must correspond to the argument name in the function definition.  Here is an example of this
        specification technique (which may be mixed with the previous technique), using the overloaded function `foo`
        from above:

            >>> @foo.overload_with(int, int, c=int, string=str)
            ... def foo(a, b, string, c=1):
            ...     return string*a*b*c
            ...
            >>> foo(2, 3, 'a')
            'aaaaaa'

        Finally, types may be specified by way of Python 3 function annotations (this is the prefered technique
        if the code is only to be used in Python 3) by using the `overload_with` decorator without arguments.

        The type of a given argument need not be an instance of `type`.  It can also be any callable that can
        be called with one argument, `None`, an `IterableOf` instance, a string, or a tuple of `type`
        any of these.  The way these work is as follows:

        * If the type is a callable that does not return `False` (or anything with a truth value
          of False) when called with the argument value, the argument is considered to be of the correct type.

        * If one of the types in a sequence is `None` and the argument is `None`, the argument is considered
          to have the correct type. [#f2]_

        * If the type is an instance if `IterableOf` (or one of its subtypes), the argument is considered to be
          the correct type if and only if all of the items in the container are considered the correct type when
          the criteria described here are applied with the `types` attribute of the `IterableOf` instance

        * If the type is a string, the string is evaluated in the context of the overloaded function's globals, and
          the return value of this evaluation is used as the type (which may, in turn, be any of the options described
          here except for a string).  This is particularly useful when the function to be overloaded is a
          method that takes as an argument an instance of that method's parent class.  Since the parent class
          will not be defined at the time of the function definition, using the parent class's name as an identifier
          will not work.

        * Finally, if a `tuple` (or any other `Iterable`) of type specifications is given, the argument is considered
          to have the correct type if any of the type specifications in the list lead to the correct type.

        For example:

            >>> def str_reversed(var):
            ...     return ''.join(reversed(var)) if isinstance(var, str) else NotImplemented
            ...
            >>> @foo.overload_with(
            ...     bar=lambda x: x == str_reversed(x),
            ...     n=(int, lambda x: isinstance(x, float) and round(x) == x)
            ... )
            ... def foo(bar, n):
            ...     return bar * int(n) + ' is a palendrome!'
            ...
            >>> @foo.overload_with(baz=str)
            ... def foo(baz, n=3):
            ...     return baz * int(n) + ' is not a palendrome.'
            ...
            >>> foo('racecar', 2.0)
            'racecarracecar is a palendrome!'
            >>> foo('asdf')
            'asdfasdfasdf is not a palendrome.'

        The above example also illustrates the order of precidence for overloaded functions:  if a call is
        allowed to multiple versions of the function, the first to be defined is used.  (Notice that the first
        call in the above example matches both overloaded versions of `foo` defined in the example).

        ..warning ::
            To allow for the use Python 3 function annotations, specifying a single argument to `overload_with` that
            is a callable will assume that the Python 3 type specification method is being used, as described above.
            This could cause problems if what is actually intended is that there be only one type-constrained
            argument which returns `True` for an arbitrary callable.  To get around this, simply specify the type
            callable for the single argument using the keyword form or as a single item tuple.  If the callable
            you want to use for this specification is a builtin such as `callable` or `int`, this shouldn't be
            an issue.

        .. rubric:: Footnotes

        .. [#f2] Note that if the type is `None` and it is not in a sequence, or a length one tuple `(None,)` is
                 given, this is considered a match for all types.  This is mostly to allow arguments to be left out
                 of the specification (in particular, the first argument of most methods, `self`)

        """
        # TODO check to make sure type specifications are valid
        if sys.version_info > (3, 0) \
                and len(kwargs) == 0 \
                and len(typespecs_or_func) == 1 \
                and inspect.isfunction(typespecs_or_func[0]) \
                and not inspect.isbuiltin(typespecs_or_func[0]):
            func = typespecs_or_func[0]
            if hasattr(func, 'func_annotations'):
                # Python 3 version with func_annotations:
                argspec = inspect.getargspec(func)
                # for now, ignore varargs and keywords annotations
                annotations = func.func_annotations
                typespec_dict = {}
                for arg in argspec:
                    typespec_dict[arg] = func.func_annotations[arg]
                self._versions.append(argtypespec(func, typespec_dict))
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
                        raise TypeError("multiple type specifications given for "
                                        "argument {} to overloaded function {}()".format(argname, self.__name__))
                    if argname not in argspec.args:
                        raise TypeError("type specifiation for unknown argument {} given to overloaded "
                                        "function {}()".format(argname, self.__name__))
                    typespec_dict[argname] = spec
                # fill in None for anything missing, indicating any type will do
                for argname in argspec.args:
                    if argname not in typespec_dict:
                        typespec_dict[argname] = None
                self._versions.append(argtypespec(f, typespec_dict))
                return self
            return _get_func_decorator
    submethod = function_alias('submethod', overload_with)
    subfunction = function_alias('subfuntion', overload_with)
    overload = function_alias('overload', overload_with)

    ###################
    # Private Methods #
    ###################

    def _build_allowed_string(self):
        return '\n'.join(typespec.signature for typespec in self._versions)

#--------------------------------------------------------------------------------#
#                       Partially Constructed Classes                            #
#--------------------------------------------------------------------------------#

class PartiallyConstructed(type):
    """
    Allow classes to be "partially constructed" by wrapping the class's `__init__`
     and `__new__` with a method that sets given attributes or appends certain keyword
     arguments to the initialzation call.
    """

    ####################
    # Class Attributes #
    ####################

    counter = 0

    ##################
    # Initialization #
    ##################

    def __new__(mcs, class_to_construct, *args, **kwargs):

        attr_set = {}

        if '_PartiallyConstructed' in class_to_construct.__name__:
            if len(args) != 0:
                raise ValueError("only keyword arguments may be given when partially constructing"
                                 " a partially constructed class.  Use the 'with_args' method"
                                 " of the existing partially constructed class if you really need"
                                 " to do this.")
            args += class_to_construct.__partial_args__
            kwargs.update(class_to_construct.__partial_kwargs__)
            attr_set.update(class_to_construct.__partial_attributes__)
            class_to_construct = class_to_construct.__partial_parent__

        def _partial_new__(cls, *args, **kwargs):
            kwargs.update(cls.__partial_kwargs__)
            ret_val = cls.__partial_parent__.__new__(cls, *(cls.__partial_args__ + args), **kwargs)
            return ret_val

        def _partial_init__(self, *args, **kwargs):
            kwargs.update(self.__partial_kwargs__)
            self.__partial_parent__.__init__(self, *(self.__partial_args__ + args), **kwargs)
            for attr, val in self.__partial_attributes__.items():
                setattr(self, attr, val)

        namespace = {
            '__new__': _partial_new__,
            '__init__': _partial_init__,
            '__partial_parent__': class_to_construct,
            'with_attributes': classmethod(mcs.__dict__['with_attributes']),
            'with_args': classmethod(mcs.__dict__['with_args']),
            'with_kwargs': classmethod(mcs.__dict__['with_kwargs']),
            '__partial_args__' : args,
            '__partial_kwargs__' : kwargs,
            '__partial_attributes__' : attr_set
        }

        cls = type.__new__(
            type(class_to_construct),
            # Name mangling--who says you don't need it in Python? :-)
            '_PartiallyConstructed__' + class_to_construct.__name__ + '_' + str(PartiallyConstructed.counter),
            (class_to_construct,),
            namespace
        )
        # Remember variables starting with a double-underscore are effectively private in that they are prefixed
        # with the class name and thus can't be 'accidentally' used by superclasses or subclasses
        cls.__attrset = attr_set
        PartiallyConstructed.counter += 1
        return cls

    def __init__(cls, cls_to_construct, *args, **kwargs):
        # Nothing to do here but call the superclass init correctly
        type.__init__(cls, cls.__name__, (cls_to_construct,), cls.__dict__)

    #################
    # Class Methods #
    #################

    def with_attributes(cls, **kwargs):
        """ Append some attributes to be set after the initialization is performed.
        Note that these override anything done in the class's __init__, so don't do
        anything stupid.
        """
        cls.__partial_attributes__.update(kwargs)

    def with_kwargs(cls, **kwargs):
        cls.__partial_kwargs__.update(kwargs)

    def with_args(cls, *args):
        cls.__partial_args__ += args

#####################
# Dependent Imports #
#####################

from grendel.util.decorators import argtypespec

