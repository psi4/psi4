from functools import wraps
import re

# TODO set up a check to prevent chained aliases
def function_alias(alias, function):
    """
    Alias a function as something else.  Set the documentation of the new method to
    "Alias for <function>()" and add "Aliased as <alias>()" to the documentation of
    `function`.  Not inteded to be used as a decorator.

    Examples
    --------
    >>> def foo(a, b, c):
    ...     '''Print a test message.'''
    ...     print "testing {0}, {1}, {2}".format(a, b, c)
    ...
    >>> foo(1, 2, 3)
    testing 1, 2, 3
    >>> bar = function_alias('bar', foo)
    >>> bar(4, 5, 6)
    testing 4, 5, 6
    >>> print bar.__doc__
    Alias for `foo()`
    >>> print foo.__doc__
    Print a test message.
    <BLANKLINE>
    Aliased as `bar()`
    >>> baz = function_alias('baz', foo)
    >>> baz(7, 8, 9)
    testing 7, 8, 9
    >>> print foo.__doc__
    Print a test message.
    <BLANKLINE>
    Aliased as `baz()`, `bar()`

    *Inside a class*

    >>> class FooBar(object):
    ...     def testfunc(self, *args):
    ...         '''Print out the args.'''
    ...         print "testing " + ', '.join(str(a) for a in args)
    ...     another_name = function_alias('another_name', testfunc)
    ...     test_func = function_alias('test_func', testfunc)
    ...
    >>> f = FooBar()
    >>> f.test_func(1, 2, 3, 4)
    testing 1, 2, 3, 4
    >>> print FooBar.test_func.__doc__
    Alias for `testfunc()`
    >>> print FooBar.testfunc.__doc__
    Print out the args.
    <BLANKLINE>
    Aliased as `test_func()`, `another_name()`


    """
    if not isinstance(alias, str):
        raise TypeError('aliases must be strings')
    # TODO Preserve signature
    def _aliased_func(*args, **kwargs):
        _aliased_function_alias_for = kwargs.pop('_aliased_function_alias_for')
        return _aliased_function_alias_for(*args, **kwargs)
    aliaspart = bindablepartial(_aliased_func, _aliased_function_alias_for=function)
    aliaspart.__name__ = alias
    aliaspart.__doc__ = "Alias for `{0}()`".format(
        function.__name__
    )
    if function.__doc__:
        if re.search(r'Aliased as `.+?`', function.__doc__):
            function.__doc__ = re.sub(r'Aliased as `(.+?)`', r'Aliased as `{0}()`, `\1`'.format(alias), function.__doc__)
        else:
            function.__doc__ += '\n\nAliased as `{0}()`'.format(alias)
    else:
        function.__doc__ = 'Aliased as `{0}()`'.format(alias)
    return aliaspart

def bindablepartial(func, *args, **keywords):
    """ Since the default python partial is not bindable, we need to make our own...
    """
    @wraps(func)
    def newfunc(*fargs, **fkeywords):
        newkeywords = keywords.copy()
        newkeywords.update(fkeywords)
        return func(*(args + fargs), **newkeywords)
    return newfunc

# TODO @language-hack this would probably be better
#class bindablepartial(MethodLike)