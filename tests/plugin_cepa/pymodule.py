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

def run_plugin_cepa(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    plugin_cepa can be called via :py:func:`driver.energy`.

    >>> energy('cepa(1)')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # throw an exception for open-shells
    if (PsiMod.get_global_option('reference') != 'RHF' ):
       PsiMod.print_out("\n")
       PsiMod.print_out("Error: %s requires \"reference rhf\".\n" % lowername )
       PsiMod.print_out("\n")
       sys.exit(1)

    # what type of cepa?
    if (lowername == 'cepa(0)'):
        PsiMod.set_global_option('cepa_level', 'cepa0')
    if (lowername == 'cepa(1)'):
        PsiMod.set_global_option('cepa_level', 'cepa1')
    if (lowername == 'cepa(2)'):
        PsiMod.set_global_option('cepa_level', 'cepa2')
    if (lowername == 'cepa(3)'):
        PsiMod.set_global_option('cepa_level', 'cepa3')
    if (lowername == 'sdci'):
        PsiMod.set_global_option('cepa_level', 'cisd')
    if (lowername == 'acpf'):
        PsiMod.set_global_option('cepa_level', 'acpf')
    if (lowername == 'aqcc'):
        PsiMod.set_global_option('cepa_level', 'aqcc')

    PsiMod.set_global_option('WFN', 'CCSD')
    energy('scf', **kwargs)

    # If the scf type is DF, then the AO integrals were never generated
    if (PsiMod.get_global_option('scf_type') == 'DF' or PsiMod.get_local_option('scf','scf_type') == 'DF'):
       mints = PsiMod.MintsHelper()
       mints.integrals()


    PsiMod.transqt2()
    PsiMod.plugin("plugin_cepa.so")

    return PsiMod.get_variable("CURRENT ENERGY")

# Integration with driver routines
procedures['energy']['sdci'] = run_plugin_cepa
procedures['energy']['acpf'] = run_plugin_cepa
procedures['energy']['aqcc'] = run_plugin_cepa
procedures['energy']['cepa(0)'] = run_plugin_cepa
procedures['energy']['cepa(1)'] = run_plugin_cepa
procedures['energy']['cepa(2)'] = run_plugin_cepa
procedures['energy']['cepa(3)'] = run_plugin_cepa

def exampleFN():
    # Your Python code goes here
    pass
