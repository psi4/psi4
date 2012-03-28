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

def run_plugin_gpu_ccsd(name, **kwargs):
    """Function encoding sequence of PSI module and plugin calls so that
    plugin_gpu_ccsd can be called via :py:func:`driver.energy`.

    >>> energy('gpu-ccsd(t)')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # triples?
    if (lowername == 'gpu-ccsd'):
        PsiMod.set_global_option('compute_triples', False)
    if (lowername == 'gpu-ccsd(t)'):
        PsiMod.set_global_option('compute_triples', True)

    PsiMod.set_global_option('WFN', 'CCSD')
    energy('scf', **kwargs)
    PsiMod.transqt2()
    PsiMod.plugin("plugin_gpu_ccsd.so")

    return PsiMod.get_variable("CURRENT ENERGY")

# Integration with driver routines
procedures['energy']['gpu-ccsd'] = run_plugin_gpu_ccsd
procedures['energy']['gpu-ccsd(t)'] = run_plugin_gpu_ccsd

def exampleFN():
    # Your Python code goes here
    pass
