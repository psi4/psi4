from copy import copy
import sys

__submodules__ = [
    'cartesian_representation',
    'representation',
    'internal_representation',
    'normal_representation'
]

__all__ = []


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


