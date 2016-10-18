"""cloud_sptheme.utils -- internal helper utilities"""
#=============================================================================
# imports
#=============================================================================
# core
from functools import update_wrapper
import logging; log = logging.getLogger(__name__)
import sys
# site
# pkg
# local
__all__ = [
    # py2/3 compat
    'PY2', 'PY3', 'u', 'ru',

    # monkeypatch helpers
    "patchapplier",
    "monkeypatch",
]

#=============================================================================
# internal py2/3 compat helpers
#=============================================================================
PY2 = sys.version_info < (3,0)
PY3 = not PY2

# FIXME: these aren't very rigorous / correct, but they work for current purposes.
if PY2:
    def u(s):
        return s.decode("unicode_escape")
    def ru(s):
        return s.decode("ascii")
else:
    def u(s):
        return s
    ru = u

#=============================================================================
# monkeypatch helpers
#=============================================================================
def patchapplier(func):
    """
    function decorator to help functions that apply a monkeypatch.
    makes them only run once.
    """
    def wrapper():
        if wrapper.patched:
            return False
        func()
        wrapper.patched = True
        logging.getLogger(func.__module__).debug("%s: patch applied", func.__name__)
        return True
    wrapper.patched = False
    update_wrapper(wrapper, func)
    return wrapper

def monkeypatch(target, name=None):
    """
    helper to monkeypatch another object.
    the decorated function is wrapped around the existing function in
    :samp:`target.{name}`, and used to replace it.

    **name** defaults to the name of the function being decorated.

    the original value is passed in as the first positional argument to the function.
    """
    def builder(func):
        attr = name or func.__name__
        wrapped = getattr(target, attr)
        def wrapper(*args, **kwds):
            return func(wrapped, *args, **kwds)
        update_wrapper(wrapper, wrapped)
        wrapper.__wrapped__ = wrapped # not set by older update_wrapper() versions
        setattr(target, attr, wrapper)
        return func # return orig func so we can use it again
    return builder


#=============================================================================
# eof
#=============================================================================
