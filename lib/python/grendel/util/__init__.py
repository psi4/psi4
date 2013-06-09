from copy import copy
import sys

# util TODOs:
# TODO Logging and useful logging messages throughout the code

__all__ = [
    'SystemInfo',
]

from grendel.util.openstruct import CaseInsensativeOpenStruct
# Struct for holding system information, such as paths and such
SystemInfo = CaseInsensativeOpenStruct()

#--------------------------------------------------------#
# 'raw lists':  allows for the use of slice literals and #
# Ellipsis literals in listsellipsis and slices in lists #
#--------------------------------------------------------#

__submodules__ = [
    'misc',
    'overloading',
    'decorators',
    'descriptors',
    'metaclasses',
    'strings',
    'exceptions',
    'units',
    'web_getter',
    'iteration',
    'parsing',
    'sentinal_values',
    'context_managers'
]

for name in __submodules__:
    __import__(__name__ + "." + name)
    m = sys.modules[__name__ + "." + name]
    globals()[name] = m
    if hasattr(m, '__all__'):
        attrlist = copy(m.__all__)
    else:
        attrlist = list(filter(lambda x: x[0]!='_', dir(m)))
    for attr in attrlist:
        globals()[attr] = getattr(m, attr)
    if hasattr(m, '__not_parent_all__'):
        for item in m.__not_parent_all__:
            attrlist.remove(item)
    __all__.extend(attrlist)

__all__.extend(__submodules__)


