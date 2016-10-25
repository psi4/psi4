"""
Miscellaneous enhancements to help autodoc along.
"""

import inspect
import sys
import types

from sphinx.ext.autodoc import AttributeDocumenter, ModuleDocumenter
from sphinx.util.inspect import isdescriptor

if sys.version_info[0] == 3:
    class_types = (type,)
else:
    class_types = (type, types.ClassType)


MethodDescriptorType = type(type.__subclasses__)


# See
# https://github.com/astropy/astropy-helpers/issues/116#issuecomment-71254836
# for further background on this.
def type_object_attrgetter(obj, attr, *defargs):
    """
    This implements an improved attrgetter for type objects (i.e. classes)
    that can handle class attributes that are implemented as properties on
    a metaclass.

    Normally `getattr` on a class with a `property` (say, "foo"), would return
    the `property` object itself.  However, if the class has a metaclass which
    *also* defines a `property` named "foo", ``getattr(cls, 'foo')`` will find
    the "foo" property on the metaclass and resolve it.  For the purposes of
    autodoc we just want to document the "foo" property defined on the class,
    not on the metaclass.

    For example::

        >>> class Meta(type):
        ...     @property
        ...     def foo(cls):
        ...         return 'foo'
        ...
        >>> class MyClass(metaclass=Meta):
        ...     @property
        ...     def foo(self):
        ...         \"\"\"Docstring for MyClass.foo property.\"\"\"
        ...         return 'myfoo'
        ...
        >>> getattr(MyClass, 'foo')
        'foo'
        >>> type_object_attrgetter(MyClass, 'foo')
        <property at 0x...>
        >>> type_object_attrgetter(MyClass, 'foo').__doc__
        'Docstring for MyClass.foo property.'

    The last line of the example shows the desired behavior for the purposes
    of autodoc.
    """

    for base in obj.__mro__:
        if attr in base.__dict__:
            if isinstance(base.__dict__[attr], property):
                # Note, this should only be used for properties--for any other
                # type of descriptor (classmethod, for example) this can mess
                # up existing expectations of what getattr(cls, ...) returns
                return base.__dict__[attr]
            break

    return getattr(obj, attr, *defargs)


# Provided to work around a bug in Sphinx
# See https://github.com/sphinx-doc/sphinx/pull/1843
class AttributeDocumenter(AttributeDocumenter):
    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        non_attr_types = cls.method_types + class_types + \
            (MethodDescriptorType,)
        isdatadesc = isdescriptor(member) and not \
            isinstance(member, non_attr_types) and not \
            type(member).__name__ == "instancemethod"
        # That last condition addresses an obscure case of C-defined
        # methods using a deprecated type in Python 3, that is not otherwise
        # exported anywhere by Python
        return isdatadesc or (not isinstance(parent, ModuleDocumenter) and
                              not inspect.isroutine(member) and
                              not isinstance(member, class_types))


def setup(app):
    # Must have the autodoc extension set up first so we can override it
    app.setup_extension('sphinx.ext.autodoc')
    # Need to import this too since it re-registers all the documenter types
    # =_=
    import sphinx.ext.autosummary.generate

    app.add_autodoc_attrgetter(type, type_object_attrgetter)

    if sphinx.version_info < (1,4,2):
        # this is a really ugly hack to supress a warning that sphinx 1.4
        # generates when overriding an existing directive (which is *desired*
        # behavior here).  As of sphinx v1.4.2, this has been fixed:
        # https://github.com/sphinx-doc/sphinx/issues/2451
        # But we leave it in for 1.4.0/1.4.1 .  But if the "needs_sphinx" is
        # eventually updated to >= 1.4.2, this should be removed entirely (in
        # favor of the line in the "else" clause)
        _oldwarn = app._warning
        _oldwarncount = app._warncount
        try:
            try:
                # *this* is in a try/finally because we don't want to force six as
                # a real dependency.  In sphinx 1.4, six is a prerequisite, so
                # there's no issue. But in older sphinxes this may not be true...
                # but the inderlying warning is absent anyway so we let it slide.
                from six import StringIO
                app._warning = StringIO()
            except ImportError:
                pass
            app.add_autodocumenter(AttributeDocumenter)
        finally:
            app._warning = _oldwarn
            app._warncount = _oldwarncount
    else:
        app.add_autodocumenter(AttributeDocumenter)

