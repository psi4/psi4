""" Miscellaneous classes that conform to the Python descriptor protocol and that don't logically belong somewhere else.
"""

# Nothing in this module thus far should be exposed to the public interface.
__all__ = []

class RaiseOnAccessDescriptor(object):
    """ A descriptor that raises an error when it gets accessed.
    Useful for optional features that require external packages but are not required for core functionality.
    A `RaiseOnAccessDescriptor` can be stored in place of some variable that requires an optional package.

    :Examples:

    >>> class AwesomeThingDoer(object):
    ...     try:
    ...         from awesome_module_that_no_one_has import thing
    ...         awesome_thing = thing()
    ...     except ImportError:
    ...         awesome_thing = RaiseOnAccessDescriptor(
    ...             ImportError, "You need a special package to continue!")
    ...
    >>> # The code can now initialize, but if someone tries to access awesome thing...
    >>> a = AwesomeThingDoer()
    >>> a.awesome_thing
    Traceback (most recent call last):
        ...
    ImportError: You need a special package to continue!

    """

    to_be_raised = None
    msg = None

    def __init__(self, exc, msg=None):
        self.to_be_raised = exc
        self.msg = msg

    def __get__(self, obj, objtype=None):
        raise self.to_be_raised(self.msg)
