import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from p4xcpt import *


def run_plugin_mp2(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    mollerplesset2 can be called via :py:func:`~driver.energy`.

    >>> energy('mollerplesset2')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.set_local_option('MOLLERPLESSET2', 'PRINT', 1)
    scf_wfn = scf_helper(lowername)
    returnvalue = psi4.plugin('mollerplesset2.so')

    return returnvalue


# Integration with driver routines
procedures['energy']['mollerplesset2'] = run_plugin_mp2


def exampleFN():
    # Your Python code goes here
    pass
