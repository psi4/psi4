from copy import copy
import sys

__submodules__ = [
    "coordinate",
    "cartesian_coordinate",
    "internal_coordinate",
    "simple_internal_coordinate",
    "normal_coordinate",
    "function_coordinate",
    "internal_cartesians",
    "bond_angle",
    "bond_length",
    "torsion"
]

__all__ = [
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

#import coordinate; from coordinate import *
#__all__.extend(coordinate.__all__)
#
#import cartesian_coordinate; from cartesian_coordinate import *
#__all__.extend(cartesian_coordinate.__all__)
#
#import internal_coordinate; from internal_coordinate import *
#__all__.extend(internal_coordinate.__all__)
#
#import symmetry_internal_coordinate; from symmetry_internal_coordinate import *
#__all__.extend(symmetry_internal_coordinate.__all__)
#
#import simple_internal_coordinate; from simple_internal_coordinate import *
#__all__.extend(simple_internal_coordinate.__all__)
#
#import normal_coordinate; from normal_coordinate import *
#__all__.extend(normal_coordinate.__all__)
#
#import bond_angle; from bond_angle import *
#__all__.extend(bond_angle.__all__)
#
#import bond_length; from bond_length import *
#__all__.extend(bond_length.__all__)
#
#import torsion; from torsion import *
#__all__.extend(torsion.__all__)

