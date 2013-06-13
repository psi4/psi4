
import sys; sys.path.append('.so INSTALL DIRECTORY GOES HERE')
from psi4 import *
from molutil import *
from driver import *


class _options(object):

    def __getattr__(self, obj):
        return get_global_option(obj)

    def __setattr__(self, obj, val):
        set_global_option(obj, val)

options = _options()
