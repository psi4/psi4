""" Exceptions for use in Grendel.
Note that this module intentionally violates the one-class-per-file rule, since exceptions are fairly trivial classes.
"""
from collections import Iterable

# Nothing in this module thus far should be exposed to the public interface.
__all__ = []

class ProgrammerNeedsMoreCoffeeError(Exception):
    """Raised when something doesn't even make sense.
    **The user should never see this exception!!!**
    In fact, this should not even get thrown in tests (and hopefully, the lines that it is thrown from should
    eventually be removed from the code).  To exclude the branches that throw this error from the coverage report,
    use the syntax::
        # pragma: no cover

    In other words, this is a placeholder for future development work that should only get thrown if something
    is done that completely doesn't make sense *but* can only be reached if the developer is modifying/messing with
    the innards of Grendel
    """

    def __init__(self, msg=''): # pragma: no cover
        message = "The programmers need more coffee, since the following non-sensical problem happened:\n"
        message += "   " + msg
        message += "\nIf you are seeing this message as a user, please find a programmer and wack them on the head with a book "
        message += "(then buy them a Cappucchino and tell them to get to work)."
        super(ProgrammerNeedsMoreCoffeeError, self).__init__(message)

class ConflictingKwargs(TypeError):
    """ Error to raise when two different keyword arguments are specified that are not allowed to be simultaneously specified.
    """

    def __init__(self, keyword1, keyword2, func=None):
        if caller is None:
            callfunc = caller()
        else:
            callfunc = func
        message = ("Function {0}() called with conflicting keyword arguments:"
                   "  It can be called with the keyword {1} or {2}, but not both.").format(callfunc, keyword1, keyword2)
        super(ConflictingKwargs, self).__init__(message)

class ArgTypeError(TypeError):
    """ Convenience error to raise when an argument is not of the right types.
    """

    def __new__(cls, argname, arg, *allowed_types, **kwargs):
        allowed = listify_args(allowed_types)
        my_caller = get_kwarg(kwargs, 'caller') or (caller() + '()')
        if len(allowed) == 1:
            return TypeError('Argument {0} to {1} must be {2} (got {3})'.format(
                argname,
                my_caller,
                classname(allowed[0]),
                type(arg).__name__
            ))
        else:
            return TypeError("Argument {0} to {1} must be one of ({2}) (got {3})".format(
                argname,
                my_caller,
                ', '.join(classname(a) for a in allowed),
                type(arg).__name__
            ))

class ChemistryError(ValueError):
    """ Used for when something doesn't make sense chemically
    """
    pass

class SuperCallMissingError(RuntimeError):
    """ Raised when a subclass fails to call super() version of a method that requires this behavior.
    """
    pass

##################
# Helper Methods #
##################

def min_arg_count(min, given):
    """ Shorthand for raising a TypeError telling the user that fewer than the minimum number of arguments was given.
    Uses the standard python form of the expression.
    """
    return TypeError(caller() + "() takes at least " + str(min)
                     + " argument" + ("" if min == 1 else "s")
                     + "(" + str(given) + " given)")

def max_arg_count(max, given):
    """ Shorthand for raising a TypeError telling the user that more than the maximum number of arguments was given.
    Uses the standard python form of the expression.
    """
    return TypeError(caller() + "() takes at most " + str(max)
                     + " argument" + ("" if max == 1 else "s")
                     + "(" + str(given) + " given)")

def exact_arg_count(num, given):
    """ Shorthand for raising a TypeError telling the user that more than the maximum number of arguments was given.
    Uses the standard python form of the expression.
    """
    return TypeError(caller() + "() takes exactly " + str(num)
                     + " argument" + ("" if num == 1 else "s")
                     + "(" + str(given) + " given)")

def raises_error(callable_obj, *args, **kwargs):
    """ Returns True if a call of `callable` with `args` raises an error, and False if not.
    Note that the callable will get called, so don't do anything that takes a long time or changes things in a way
    you don't want them to be changed.
    If the optional keyword argument 'error' is given as either an Exception subclass or an Iterable, `raises_error`
    returns True if the call of callable raises the error given or one of the errors given in the list.  (The keyword
    argument can also be named 'errors').  Any other keyword arguments are passed through to `callable`.

    *Technical note:  As per the Python manual, only exceptions that subclass from Exception (and not BaseException
    directly) will be recognized.  According to the python users manual, you should never implement a user exception
    that subclasses from BaseException directly.*

    :Examples:


    >>> raises_error(int, "5")
    False
    >>> raises_error(int, "0x5ab7", 0)
    False
    >>> raises_error(lambda x: int(x), "5")
    False
    >>> raises_error(int, "abc")
    True
    >>> raises_error(int, "abc", error = ValueError)
    True
    >>> raises_error(int, "abc", error = BufferError)
    False
    >>> raises_error(int, "abc", error = (BufferError, EnvironmentError))
    False
    >>> raises_error(int, "abc", error = [BufferError, ValueError])
    True
    >>> raises_error(float, "abc")
    True
    >>> # Use raises_error to see what errors raises_error raises
    ... # Raise a type error if the first argument is not a callable (i.e. 17(25) doesn't make sense to Python)
    ... raises_error(raises_error, 17, 25, error = TypeError)
    True
    >>> # Raise a TypeError if the error keyword argument is not a subclass of Exception
    ... raises_error(lambda x: raises_error(int, x, error = "hello world"), 25, error = TypeError)
    True
    >>> raises_error(lambda x: raises_error(int, x, error = BaseException), 25, error = TypeError)
    True
    """

    if "error" in kwargs and "errors" in kwargs:
        raise TypeError("Come on now.  You're just trying to mess things up."
                        "  Please don't give raises_error keyword arguments of 'error' and 'errors' at the same time.")

    errors = kwargs.pop("error") if "error" in kwargs else kwargs.pop("errors") if "errors" in kwargs else None

    if not hasattr(callable_obj, '__call__'):
        raise TypeError(repr(callable_obj) + " is not callable.")

    if errors is None:
        try:
            callable_obj(*args, **kwargs)
        except Exception:
            return True

        return False

    elif isinstance(errors, Iterable):
        for error in errors:
            if not issubclass(error, Exception): raise TypeError(repr(error) + " is not a subclass of Exception or Iterable.")
            try:
                callable_obj(*args, **kwargs)
                return False
            except error:
                return True
            except Exception:
                continue

        return False

    elif isinstance(errors, type) and issubclass(errors, Exception):
        try:
            callable_obj(*args, **kwargs)
            return False
        except errors:
            return True
        except Exception:
            return False

    else:
        raise TypeError(repr(errors) + " is not a subclass of Exception or Iterable.")

def class_property_type(attr, cls, types, got):
    type_error_text = "Property '{0}'".format(attr)
    type_error_text += " in class '{0}' "
    if isinstance(types, Iterable):
        second_part = "must be an instance of one of the following: ('{0}');".format("', '".join(map(classname, types)))
    else:
        second_part = "must be an instance of '{0}'".format(classname(types))
    type_error_text += second_part + " (got '{1}')"
    type_error_text.format(classname(cls), classname(got))

#####################
# Dependent Imports #
#####################

from grendel.util.overloading import caller, listify_args, get_kwarg
from grendel.util.strings import classname

