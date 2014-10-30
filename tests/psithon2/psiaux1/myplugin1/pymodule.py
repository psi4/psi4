import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
from p4const import *
from p4util import *
from p4xcpt import *



def run_myplugin1(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    myplugin1 can be called via :py:func:`~driver.energy`.

    >>> energy('myplugin1')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.set_global_option('BASIS', 'sto-3g')
    psi4.set_local_option('MYPLUGIN1', 'PRINT', 1)
    energy('scf', **kwargs)
    returnvalue = psi4.plugin('myplugin1.so')

    return returnvalue


def exampleFN(name, **kwargs):
    psi4.set_variable('CURRENT ENERGY', -74.94550962)
    # Your Python code goes here
    pass


# Integration with driver routines
procedures['energy']['myplugin1'] = exampleFN

