import PsiMod
import re
import os
import input
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
from text import *
from procutil import *

def run_plugin_mp2sort(name, **kwargs):
    """Function encoding sequence of PSI module and plugin calls so that
    plugin_mp2 can be called via :py:func:`driver.energy`.

    >>> energy('plugin_mp2sort')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    energy('scf', **kwargs)

    # If the scf type is DF, then the AO integrals were never generated
    if (PsiMod.get_global_option('scf_type') == 'DF' or PsiMod.get_local_option('scf','scf_type') == 'DF'):
       mints = PsiMod.MintsHelper()
       mints.integrals()

    PsiMod.set_global_option('WFN', 'MP2')
    PsiMod.plugin("plugin_mp2sort.so")
    PsiMod.ccsort()
    PsiMod.mp2()

    return PsiMod.get_variable("CURRENT ENERGY")

# Integration with driver routines
procedures['energy']['plugin_mp2sort'] = run_plugin_mp2sort

def exampleFN():
    # Your Python code goes here
    pass
