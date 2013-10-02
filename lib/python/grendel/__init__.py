from __future__ import print_function
import os
import atexit
import sys
import traceback

if sys.version_info < (2, 7):
    raise SystemError("grendel needs at least Python version 2.7")
from copy import copy


#--------------------------------------------------------------------------------#
# Helper function
def get_env_bool(env_var_name, default=False):
    var = os.getenv(env_var_name, default)
    if var in [0, '0', False, 'False', 'FALSE', 'false', 'off', 'OFF', 'Off', 'NO', 'no', 'No']:
        return False
    elif var in [1, '1', True, 'True', 'TRUE', 'true', 'on', 'ON', 'On', 'YES', 'yes', 'Yes']:
        return True
    else:
        raise OSError("invalid value for environment variable {}: '{}'.\n"
                      "Please give something like 'yes' or 'no', 'true' or 'false',"
                      " '1' or '0'".format(env_var_name, var))
#========================================#
# Package scope variables from the environment...
type_checking_enabled = get_env_bool('GRENDEL_TYPE_CHECK', True)
sanity_checking_enabled = get_env_bool('GRENDEL_SANITY_CHECK', True)
caching_enabled = not get_env_bool('GRENDEL_NO_CACHE', False)
show_warnings = get_env_bool('GRENDEL_SHOW_WARNINGS', False)
#--------------------------------------------------------------------------------#

# Uhh...technically in this case this could(/should?) be called __subpackages__
#TODO get this as a list of files/folders in the directory?
__submodules__ = [
    "util",
    "gmath",
    "differentiation",
    "coordinates",
    "representations",
    "chemistry",
    "symmetry",
    "interface",
]

__all__ = [
]

# TODO Propagate this pattern throughout grendel packages
for name in __submodules__:
    __import__("grendel." + name)
    m = sys.modules["grendel." + name]
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

#--------------------------------------------------------------------------------#

# TODO save/pickle expensive runtime information here in case of exit on error and provide some convenient/reasonable interface for loading said data
#@atexit.register
#def grendel_exit():
#    pass


