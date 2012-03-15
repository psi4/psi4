"""Module with functions that call upon those in modules
:py:mod:`proc`, :py:mod:`driver`, and :py:mod:`wrappers`.

Place in this file quickly defined procedures such as
   - aliases for complex methods
   - simple modifications to existing methods

"""
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

# Python procedures like these can be run directly from the input file or integrated
#   with the energy(), etc. routines by means of lines like those at the end of this file.


def sherrillgroup_gold_standard(name='mp2', **kwargs):
    """Function to call the quantum chemical method known as 'Gold Standard'
    in the Sherrill group. Uses :py:func:`wrappers.complete_basis_set` to evaluateo
    the following expression. Two-point extrapolation of the correlation energy
    performed according to :py:func:`wrappers.corl_xtpl_helgaker_2`.

    .. math:: E_{total}^{\\text{Au\_std}} = E_{total,\; \\text{SCF}}^{\\text{aug-cc-pVQZ}} \; + E_{corl,\; \\text{MP2}}^{\\text{aug-cc-pV[TQ]Z}} \; + \delta_{\\text{MP2}}^{\\text{CCSD(T)}}\\big\\vert_{\\text{aug-cc-pVTZ}}

    >>> energy('sherrillgroup_gold_standard')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    if not ('func_cbs' in kwargs):
        kwargs['func_cbs'] = energy

    if not ('scf_basis' in kwargs):
        kwargs['scf_basis'] = 'aug-cc-pVQZ'
    if not ('scf_scheme' in kwargs):
        kwargs['scf_scheme'] = highest_1

    if not ('corl_wfn' in kwargs):
        kwargs['corl_wfn'] = 'mp2'
    if not ('corl_basis' in kwargs):
        kwargs['corl_basis'] = 'aug-cc-pV[TQ]Z'
    if not ('corl_scheme' in kwargs):
        kwargs['corl_scheme'] = corl_xtpl_helgaker_2

    if not ('delta_wfn' in kwargs):
        kwargs['delta_wfn'] = 'ccsd(t)'
    if not ('delta_wfn_lesser' in kwargs):
        kwargs['delta_wfn_lesser'] = 'mp2'
    if not ('delta_basis' in kwargs):
        kwargs['delta_basis'] = 'aug-cc-pVTZ'
    if not ('delta_scheme' in kwargs):
        kwargs['delta_scheme'] = highest_1

    return cbs(name, **kwargs)


def run_mp2_5(name, **kwargs):
    """Function that computes MP2.5 energy from results of a DETCI
    MP3 calculation.

    >>> energy('mp2.5')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    # Run detci calculation and collect conventional quantities
    energy('mp3', **kwargs)
    e_scf = PsiMod.get_variable('SCF TOTAL ENERGY')
    ce_mp2 = PsiMod.get_variable('MP2 CORRELATION ENERGY')
    ce_mp3 = PsiMod.get_variable('MP3 CORRELATION ENERGY')
    e_mp2 = e_scf + ce_mp2
    e_mp3 = e_scf + ce_mp3

    # Compute quantities particular to MP2.5
    ce_mp25 = 0.5 * (ce_mp2 + ce_mp3)
    e_mp25 = e_scf + ce_mp25
    PsiMod.set_variable('MP2.5 CORRELATION ENERGY', ce_mp25)
    PsiMod.set_variable('MP2.5 TOTAL ENERGY', e_mp25)
    PsiMod.set_variable('CURRENT CORRELATION ENERGY', ce_mp25)
    PsiMod.set_variable('CURRENT ENERGY', e_mp25)

    # build string of title banner and print results
    banners = ''
    banners += """PsiMod.print_out('\\n')\n"""
    banners += """banner(' MP2.5 ')\n"""
    banners += """PsiMod.print_out('\\n')\n\n"""
    exec banners

    tables = ''
    tables += """  SCF total energy:                        %16.8f\n""" % (e_scf)
    tables += """  MP2 total energy:                        %16.8f\n""" % (e_mp2)
    tables += """  MP2.5 total energy:                      %16.8f\n""" % (e_mp25)
    tables += """  MP3 total energy:                        %16.8f\n\n""" % (e_mp3)
    tables += """  MP2 correlation energy:                  %16.8f\n""" % (ce_mp2)
    tables += """  MP2.5 correlation energy:                %16.8f\n""" % (ce_mp25)
    tables += """  MP3 correlation energy:                  %16.8f\n""" % (ce_mp3)
    PsiMod.print_out(tables)

    return e_mp25


def run_plugin_ccsd_serial(name, **kwargs):
    """Function encoding sequence of PSI module and plugin calls so that
    Eugene DePrince's ccsd_serial plugin can be called via :py:func:`driver.energy`.

    >>> energy('plugin_ccsd_serial')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    plugfile = PsiMod.Process.environment["PSIDATADIR"] + "/../tests/plugin_ccsd_serial/plugin_ccsd_serial.so"
    PsiMod.plugin_load("%s" % (plugfile))
    PsiMod.set_global_option('WFN', 'CCSD')
    run_scf("scf", **kwargs)
    PsiMod.transqt2()
    PsiMod.plugin("plugin_ccsd_serial.so")

    return PsiMod.get_variable("CURRENT ENERGY")


# A direct translation of a plugin input file into a function call. Function calls are the only
#     way to call plugins in sow/reap mode for db(), opt(), etc. This isn't best practices
#     (see run_plugin_serial_ccsd) but is an example of what to do for a more complicated
#     procedure where different options are set for different qc steps.
def run_plugin_omega(name, **kwargs):
    """Function encoding sequence of PSI module and plugin calls, as well
    as typical options, to access Rob Parrish's omega plugin.

    >>> energy('plugin_omega')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    plugfile = PsiMod.Process.environment["PSIDATADIR"] + "/../tests/plugin_omega/plugin_omega.so"
    PsiMod.plugin_load("%s" % (plugfile))

    PsiMod.set_global_option('BASIS', 'AUG-CC-PVDZ')
    PsiMod.set_global_option('DF_BASIS_SCF', 'AUG-CC-PVDZ-RI')
    PsiMod.set_global_option('REFERENCE', 'UHF')
    PsiMod.set_global_option('SCF_TYPE', 'DF')
    energy('scf', **kwargs)

    PsiMod.set_global_option('dft_functional', 'wB97')
    PsiMod.set_global_option('dft_order_spherical', 25)
    PsiMod.set_global_option('dft_num_radial', 35)
    PsiMod.set_global_option('omega_procedure', 'ip')
    PsiMod.set_global_option('maxiter', 50)
    PsiMod.set_global_option('d_convergence', 5)
    PsiMod.set_global_option('e_convergence', 7)
    PsiMod.plugin("plugin_omega.so")

    return PsiMod.get_variable('SCF TOTAL ENERGY')

def run_cim(name, **kwargs):
   """Eugene's CIM driven by Python"""

   lowername = name.lower()
   kwargs = kwargs_lower(kwargs)

   plugfile1 = PsiMod.Process.environment["PSIDATADIR"] + "/../tests/plugin_localcc/plugin_localcc.so"
   PsiMod.plugin_load("%s" % (plugfile1))

   plugfile = PsiMod.Process.environment["PSIDATADIR"] + "/../tests/plugin_libcim/plugin_libcim.so"
   PsiMod.plugin_load("%s" % (plugfile))

   # override symmetry:
   molecule = PsiMod.get_active_molecule()
   molecule.update_geometry()
   molecule.reset_point_group('c1')
   molecule.fix_orientation(1)
   molecule.update_geometry()

   # what type of cim?
   if (name.lower() == 'ccsd'):
      PsiMod.set_global_option('compute_triples', False)
   if (name.lower() == 'ccsd(t)'):
      PsiMod.set_global_option('compute_triples', True)

   # some options aren't set right...
   PsiMod.set_global_option('r_convergence', 1e-7)
   PsiMod.set_global_option('maxiter', 100)

   energy('scf', **kwargs)
   PsiMod.set_global_option('cim_initialize', True)
   PsiMod.plugin("plugin_libcim.so")

   cluster_ccsd  = []
   cluster_ccsdt = []
   cluster_t     = []
   built_ccsd   = 0.0
   built_t      = 0.0
   built_ccsdt  = 0.0
   built_energy = 0.0
   escf = PsiMod.get_variable('SCF TOTAL ENERGY')
   PsiMod.set_global_option('cim_initialize', False)
   cim_n = 0
   while cim_n < PsiMod.get_variable('CIM CLUSTERS'):
       # run plugin_libcim:
       PsiMod.set_global_option('CIM_CLUSTER_NUM', cim_n)
       PsiMod.plugin("plugin_libcim.so")

       # accumulate correlation energies
       cluster_ccsd.append(PsiMod.get_variable('CURRENT CLUSTER CCSD CORRELATION ENERGY'))
       built_ccsd  += PsiMod.get_variable('CURRENT CLUSTER CCSD CORRELATION ENERGY')
       built_energy = built_ccsd
       if (name.lower()=='ccsd(t)'):
          cluster_ccsdt.append(PsiMod.get_variable('CURRENT CLUSTER CCSD(T) CORRELATION ENERGY'))
          cluster_t.append(PsiMod.get_variable('CURRENT CLUSTER (T) CORRELATION ENERGY'))
          built_ccsdt += PsiMod.get_variable('CURRENT CLUSTER CCSD(T) CORRELATION ENERGY')
          built_t     += PsiMod.get_variable('CURRENT CLUSTER (T) CORRELATION ENERGY')
          built_energy = built_ccsdt

       cim_n += 1

   PsiMod.set_variable('CURRENT ENERGY', built_energy+escf)
   PsiMod.set_variable('TOTAL CIM-CCSD CORRELATION ENERGY', built_ccsd)
   if (name.lower()=='ccsd(t)'):
      PsiMod.set_variable('TOTAL CIM-CCSD(T) CORRELATION ENERGY', built_ccsdt)
      PsiMod.set_variable('TOTAL CIM-(T) CORRELATION ENERGY', built_t)

   return built_energy+escf


# Integration with driver routines
procedures['energy']['mp2.5'] = run_mp2_5
procedures['energy']['sherrillgroup_gold_standard'] = sherrillgroup_gold_standard
procedures['energy']['plugin_ccsd_serial'] = run_plugin_ccsd_serial
procedures['energy']['plugin_omega'] = run_plugin_omega
