import PsiMod
import re
import os
import input
import math
import warnings
from driver import *
from wrappers import *
from molecule import *
from text import *

# Place in this file quickly defined procedures such as
#   (1) aliases for complex methods
#   (2) simple modifications to existing methods
#
# Python procedures like these can be run directly from the input file or integrated
#   with the energy(), etc. routines by means of lines like those at the end of this file.

def sherrillgroup_gold_standard(name='mp2', **kwargs):

    if not (kwargs.has_key('func_cbs')):  kwargs['func_cbs']                 = energy

    if not (kwargs.has_key('scf_basis')):  kwargs['scf_basis']               = 'aug-cc-pVQZ'
    if not (kwargs.has_key('scf_scheme')):  kwargs['scf_scheme']             = highest_1

    if not (kwargs.has_key('corl_wfn')):  kwargs['corl_wfn']                 = 'mp2'
    if not (kwargs.has_key('corl_basis')):  kwargs['corl_basis']             = 'aug-cc-pV[TQ]Z'
    if not (kwargs.has_key('corl_scheme')):  kwargs['corl_scheme']           = corl_xtpl_helgaker_2

    if not (kwargs.has_key('delta_wfn')):  kwargs['delta_wfn']               = 'ccsd(t)'
    if not (kwargs.has_key('delta_wfn_lesser')):  kwargs['delta_wfn_lesser'] = 'mp2'
    if not (kwargs.has_key('delta_basis')):  kwargs['delta_basis']           = 'aug-cc-pVTZ'
    if not (kwargs.has_key('delta_scheme')):  kwargs['delta_scheme']         = highest_1

    return cbs(name, **kwargs)

def run_mp2_5(name, **kwargs):

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

    tables  = ''
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
    plugfile = PsiMod.Process.environment["PSIDATADIR"] + "/../tests/plugin_ccsd_serial/plugin_ccsd_serial.so"
    PsiMod.plugin_load("%s" % (plugfile))
    run_scf("scf",**kwargs)
    PsiMod.transqt2()
    PsiMod.plugin("plugin_ccsd_serial.so")

    return PsiMod.get_variable("CURRENT ENERGY")
    
# A direct translation of a plugin input file into a function call. Function calls are the only
#     way to call plugins in sow/reap mode for db(), opt(), etc. This isn't best practices
#     (see run_plugin_serial_ccsd) but is an example of what to do for a more complicated
#     procedure where different options are set for different qc steps.
def run_plugin_omega(name, **kwargs):

    plugfile = PsiMod.Process.environment["PSIDATADIR"] + "/../tests/plugin_omega/plugin_omega.so"
    PsiMod.plugin_load("%s" % (plugfile))

    PsiMod.set_global_option('BASIS', 'AUG-CC-PVDZ')
    PsiMod.set_global_option('RI_BASIS_SCF', 'AUG-CC-PVDZ-RI')
    PsiMod.set_global_option('REFERENCE', 'UHF')
    PsiMod.set_global_option('SCF_TYPE', 'DF')
    energy('scf', **kwargs)

    PsiMod.set_global_option('dft_functional', 'wB97')
    PsiMod.set_global_option('dft_order_spherical', 25) 
    PsiMod.set_global_option('dft_n_radial', 35) 
    PsiMod.set_global_option('omega_procedure', 'ip')
    PsiMod.set_global_option('maxiter', 50)
    PsiMod.set_global_option('d_converge', 5)
    PsiMod.set_global_option('e_converge', 7)
    PsiMod.plugin("plugin_omega.so")
    
    return PsiMod.get_variable('SCF TOTAL ENERGY')


# Integration with driver routines
procedures['energy']['mp2.5'] = run_mp2_5
procedures['energy']['sherrillgroup_gold_standard'] = sherrillgroup_gold_standard
procedures['energy']['plugin_ccsd_serial'] = run_plugin_ccsd_serial
procedures['energy']['plugin_omega'] = run_plugin_omega

