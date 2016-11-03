import inspect
import sys
import re
import os
from warnings import warn

from sphinx.ext.autosummary.generate import find_autosummary_in_docstring


def find_mod_objs(modname, onlylocals=False):
    """ Returns all the public attributes of a module referenced by name.

    .. note::
        The returned list *not* include subpackages or modules of
        `modname`,nor does it include private attributes (those that
        beginwith '_' or are not in `__all__`).

    Parameters
    ----------
    modname : str
        The name of the module to search.
    onlylocals : bool
        If True, only attributes that are either members of `modname` OR one of
        its modules or subpackages will be included.

    Returns
    -------
    localnames : list of str
        A list of the names of the attributes as they are named in the
        module `modname` .
    fqnames : list of str
        A list of the full qualified names of the attributes (e.g.,
        ``astropy.utils.misc.find_mod_objs``). For attributes that are
        simple variables, this is based on the local name, but for
        functions or classes it can be different if they are actually
        defined elsewhere and just referenced in `modname`.
    objs : list of objects
        A list of the actual attributes themselves (in the same order as
        the other arguments)

    """

    __import__(modname)
    mod = sys.modules[modname]

    if hasattr(mod, '__all__'):
        pkgitems = [(k, mod.__dict__[k]) for k in mod.__all__]
    else:
        pkgitems = [(k, mod.__dict__[k]) for k in dir(mod) if k[0] != '_']

    # filter out modules and pull the names and objs out
    ismodule = inspect.ismodule
    localnames = [k for k, v in pkgitems if not ismodule(v)]
    objs = [v for k, v in pkgitems if not ismodule(v)]

    # fully qualified names can be determined from the object's module
    fqnames = []
    for obj, lnm in zip(objs, localnames):
        if hasattr(obj, '__module__') and hasattr(obj, '__name__'):
            fqnames.append(obj.__module__ + '.' + obj.__name__)
        else:
            fqnames.append(modname + '.' + lnm)

    if onlylocals:
        valids = [fqn.startswith(modname) for fqn in fqnames]
        localnames = [e for i, e in enumerate(localnames) if valids[i]]
        fqnames = [e for i, e in enumerate(fqnames) if valids[i]]
        objs = [e for i, e in enumerate(objs) if valids[i]]

    return localnames, fqnames, objs


def find_autosummary_in_lines_for_automodsumm(lines, module=None, filename=None):
    """Find out what items appear in autosummary:: directives in the
    given lines.
    Returns a list of (name, toctree, template, inherited_members)
    where *name* is a name
    of an object and *toctree* the :toctree: path of the corresponding
    autosummary directive (relative to the root of the file name),
    *template* the value of the :template: option, and *inherited_members*
    is the value of the :inherited-members: option.
    *toctree*, *template*, and *inherited_members* are ``None`` if the
    directive does not have the corresponding options set.

    .. note::

       This is a slightly modified version of
       ``sphinx.ext.autosummary.generate.find_autosummary_in_lines``
       which recognizes the ``inherited-members`` option.
    """
    autosummary_re = re.compile(r'^(\s*)\.\.\s+autosummary::\s*')
    automodule_re = re.compile(
        r'^\s*\.\.\s+automodule::\s*([A-Za-z0-9_.]+)\s*$')
    module_re = re.compile(
        r'^\s*\.\.\s+(current)?module::\s*([a-zA-Z0-9_.]+)\s*$')
    autosummary_item_re = re.compile(r'^\s+(~?[_a-zA-Z][a-zA-Z0-9_.]*)\s*.*?')
    toctree_arg_re = re.compile(r'^\s+:toctree:\s*(.*?)\s*$')
    template_arg_re = re.compile(r'^\s+:template:\s*(.*?)\s*$')
    inherited_members_arg_re = re.compile(r'^\s+:inherited-members:\s*$')
    no_inherited_members_arg_re = re.compile(r'^\s+:no-inherited-members:\s*$')

    documented = []

    toctree = None
    template = None
    inherited_members = None
    current_module = module
    in_autosummary = False
    base_indent = ""

    for line in lines:
        if in_autosummary:
            m = toctree_arg_re.match(line)
            if m:
                toctree = m.group(1)
                if filename:
                    toctree = os.path.join(os.path.dirname(filename),
                                           toctree)
                continue

            m = template_arg_re.match(line)
            if m:
                template = m.group(1).strip()
                continue

            m = inherited_members_arg_re.match(line)
            if m:
                inherited_members = True
                continue

            m = no_inherited_members_arg_re.match(line)
            if m:
                inherited_members = False
                continue

            if line.strip().startswith(':'):
                warn(line)
                continue  # skip options

            m = autosummary_item_re.match(line)
            if m:
                name = m.group(1).strip()
                if name.startswith('~'):
                    name = name[1:]
                if current_module and \
                   not name.startswith(current_module + '.'):
                    name = "%s.%s" % (current_module, name)
                documented.append((name, toctree, template,
                                   inherited_members))
                continue

            if not line.strip() or line.startswith(base_indent + " "):
                continue

            in_autosummary = False

        m = autosummary_re.match(line)
        if m:
            in_autosummary = True
            base_indent = m.group(1)
            toctree = None
            template = None
            continue

        m = automodule_re.search(line)
        if m:
            current_module = m.group(1).strip()
            # recurse into the automodule docstring
            documented.extend(find_autosummary_in_docstring(
                current_module, filename=filename))
            continue

        m = module_re.match(line)
        if m:
            current_module = m.group(2)
            continue

    return documented

