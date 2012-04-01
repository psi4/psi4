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

def run_plugin_df_ccsd(name, **kwargs):
    """Function encoding sequence of PSI module and plugin calls so that
    plugin_df_ccsd can be called via :py:func:`driver.energy`.

    >>> energy('df-ccsd(t)')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # override symmetry:
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    molecule.reset_point_group('c1')
    molecule.fix_orientation(1)
    molecule.update_geometry()

    # triples?
    if (lowername == 'df-ccsd'):
        PsiMod.set_global_option('compute_triples', False)
    if (lowername == 'df-ccsd(t)'):
        PsiMod.set_global_option('compute_triples', True)

    energy('scf', **kwargs)
    PsiMod.plugin("plugin_df_ccsd.so")

    return PsiMod.get_variable("CURRENT ENERGY")

# Integration with driver routines
procedures['energy']['df-ccsd'] = run_plugin_df_ccsd
procedures['energy']['df-ccsd(t)'] = run_plugin_df_ccsd

def exampleFN():
    # Your Python code goes here
    pass
