import sys
from copy import copy

__submodules__ = [
    'computation',
    'computation_details',
    'result_getter',
    'input_generator',
    'output_parser',
    'queue',
    'local_parallel_queue',
    'result_getter',
    'runner',
    'legacy_xml'
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

#TODO ComputationCache of some sort?
#  A ComputationCache would be an SQL or similar database that would allow grendel to try and retrieve already-run calculations
#  automatically before running a computation.  The cache could have a cache_level variable that could take on the
#  following values:
#       0:  don't retrieve any cached computations
#       1:  only if input file is an exact match
#       2:  retrieve compatible computations from same program
#       3:  retrieve compatible computations from any program
#       4:  retrieve compatible or better computations from same program
#       5:  retrieve compatible or better computations from any program

