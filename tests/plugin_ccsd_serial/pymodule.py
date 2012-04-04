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

def run_plugin_ccsd_serial(name, **kwargs):
    """Function encoding sequence of PSI module and plugin calls so that
    plugin_ccsd_serial can be called via :py:func:`driver.energy`.

    >>> energy('_ccsd(t)')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # triples?
    if (lowername == '_ccsd'):
        PsiMod.set_global_option('compute_triples', False)
    if (lowername == '_ccsd(t)'):
        PsiMod.set_global_option('compute_triples', True)

    PsiMod.set_global_option('WFN', 'CCSD')
    energy('scf', **kwargs)
    PsiMod.transqt2()
    PsiMod.plugin("plugin_ccsd_serial.so")

    return PsiMod.get_variable("CURRENT ENERGY")

# Integration with driver routines
procedures['energy']['_ccsd'] = run_plugin_ccsd_serial
procedures['energy']['_ccsd(t)'] = run_plugin_ccsd_serial

def exampleFN():
    # Your Python code goes here
    pass
