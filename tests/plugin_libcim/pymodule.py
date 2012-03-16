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

#def run_plugin_libcim(name, **kwargs):
#    """Function encoding sequence of PSI module an plugin calls so that
#    plugin_libcim can be called via :py:func:`driver.energy`.
#
#    >>> energy('plugin_libcim')
#
#    """
#    lowername = name.lower()
#    kwargs = kwargs_lower(kwargs)
#
#    # Your plugin's PsiMod run sequence goes here
#    PsiMod.set_global_option('BASIS', 'sto-3g')
#    PsiMod.set_local_option('PLUGIN_LIBCIM', 'PRINT', 1)
#    energy('scf', **kwargs)
#    returnvalue = PsiMod.plugin('plugin_libcim.so')
#
#    return returnvalue


def run_plugin_libcim(name, **kwargs):
    """Function encoding sequence of PSI module and plugin calls so that
    Eugene's plugin_libcim can be called via :py:func:`driver.energy`.

    >>> energy('cim-ccsd(t)')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # override symmetry:
    molecule = PsiMod.get_active_molecule()
    molecule.update_geometry()
    molecule.reset_point_group('c1')
    molecule.fix_orientation(1)
    molecule.update_geometry()

    # what type of cim?
    if (lowername == 'cim-ccsd'):
        PsiMod.set_global_option('compute_triples', False)
    if (lowername == 'cim-ccsd(t)'):
        PsiMod.set_global_option('compute_triples', True)

    # some options are not correct when i load two plugins ... set them here
    PsiMod.set_global_option('r_convergence', 1e-7)
    PsiMod.set_global_option('maxiter', 100)

    energy('scf', **kwargs)
    PsiMod.set_global_option('cim_initialize', True)
    PsiMod.plugin('plugin_libcim.so')

    cluster_ccsd = []
    cluster_ccsdt = []
    cluster_t = []
    built_ccsd = 0.0
    built_t = 0.0
    built_ccsdt = 0.0
    built_energy = 0.0
    escf = PsiMod.get_variable('SCF TOTAL ENERGY')
    PsiMod.set_global_option('cim_initialize', False)
    cim_n = 0
    while cim_n < PsiMod.get_variable('CIM CLUSTERS'):
        # run plugin_libcim:
        PsiMod.set_global_option('CIM_CLUSTER_NUM', cim_n)
        PsiMod.plugin('plugin_libcim.so')

        # accumulate correlation energies
        cluster_ccsd.append(PsiMod.get_variable('CURRENT CLUSTER CCSD CORRELATION ENERGY'))
        built_ccsd += PsiMod.get_variable('CURRENT CLUSTER CCSD CORRELATION ENERGY')
        built_energy = built_ccsd
        if (lowername == 'cim-ccsd(t)'):
            cluster_ccsdt.append(PsiMod.get_variable('CURRENT CLUSTER CCSD(T) CORRELATION ENERGY'))
            cluster_t.append(PsiMod.get_variable('CURRENT CLUSTER (T) CORRELATION ENERGY'))
            built_ccsdt += PsiMod.get_variable('CURRENT CLUSTER CCSD(T) CORRELATION ENERGY')
            built_t += PsiMod.get_variable('CURRENT CLUSTER (T) CORRELATION ENERGY')
            built_energy = built_ccsdt

        cim_n += 1

    PsiMod.set_variable('CURRENT ENERGY', built_energy + escf)
    PsiMod.set_variable('CIM-CCSD CORRELATION ENERGY', built_ccsd)
    PsiMod.set_variable('CIM-CCSD TOTAL ENERGY', built_ccsd + escf)
    if (lowername == 'cim-ccsd(t)'):
        PsiMod.set_variable('CIM-CCSD(T) CORRELATION ENERGY', built_ccsdt)
        PsiMod.set_variable('CIM-CCSD(T) TOTAL ENERGY', built_ccsdt + escf)
        PsiMod.set_variable('CIM-(T) CORRELATION ENERGY', built_t)

    PsiMod.print_out('\n')
    PsiMod.print_out('        CIM-CCSD correlation energy    %20.12lf\n' % PsiMod.get_variable('CIM-CCSD CORRELATION ENERGY'))
    PsiMod.print_out('      * CIM-CCSD total energy          %20.12lf\n' % PsiMod.get_variable('CIM-CCSD TOTAL ENERGY'))
    PsiMod.print_out('\n')
    if (lowername == 'cim-ccsd(t)'):
       PsiMod.print_out('        CIM-(T) correlation energy     %20.12lf\n' % PsiMod.get_variable('CIM-(T) CORRELATION ENERGY'))
       PsiMod.print_out('        CIM-CCSD(T) correlation energy %20.12lf\n' % PsiMod.get_variable('CIM-CCSD(T) CORRELATION ENERGY'))
       PsiMod.print_out('      * CIM-CCSD(T) total energy       %20.12lf\n' % PsiMod.get_variable('CIM-CCSD(T) TOTAL ENERGY'))
       PsiMod.print_out('\n')

    return built_energy + escf

# Integration with driver routines
procedures['energy']['cim-ccsd(t)'] = run_plugin_libcim
procedures['energy']['cim-ccsd'] = run_plugin_libcim

def exampleFN():
    # Your Python code goes here
    pass
