#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#@END LICENSE
#

"""Module with functions that call upon those in modules
:py:mod:`proc`, :py:mod:`driver`, and :py:mod:`wrappers`.

Place in this file quickly defined procedures such as
   - aliases for complex methods
   - simple modifications to existing methods

"""
import PsiMod
import re
import os
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
    r"""Function to call the quantum chemical method known as 'Gold Standard'
    in the Sherrill group. Uses :py:func:`~wrappers.complete_basis_set` to evaluate
    the following expression. Two-point extrapolation of the correlation energy
    performed according to :py:func:`~wrappers.corl_xtpl_helgaker_2`.

    .. math:: E_{total}^{\text{Au\_std}} = E_{total,\; \text{SCF}}^{\text{aug-cc-pVQZ}} \; + E_{corl,\; \text{MP2}}^{\text{aug-cc-pV[TQ]Z}} \; + \delta_{\text{MP2}}^{\text{CCSD(T)}}\big\vert_{\text{aug-cc-pVTZ}}

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


def allen_focal_point(name='mp2', **kwargs):
    r"""Function to call Wes Allen-style Focal
    Point Analysis. Insert typical reference here.  Uses
    :py:func:`~wrappers.complete_basis_set` to evaluate the following
    expression. SCF employs a three-point extrapolation according
    to :py:func:`~wrappers.scf_xtpl_helgaker_3`. MP2, CCSD, and
    CCSD(T) employ two-point extrapolation performed according to
    :py:func:`~wrappers.corl_xtpl_helgaker_2`.  CCSDT and CCSDT(Q)
    are plain deltas.

    .. math:: E_{total}^{\text{FPA}} = E_{total,\; \text{SCF}}^{\text{cc-pV[Q56]Z}} \; + E_{corl,\; \text{MP2}}^{\text{cc-pV[56]Z}} \; + \delta_{\text{MP2}}^{\text{CCSD)}}\big\vert_{\text{cc-pV[56]Z}} \; + \delta_{\text{CCSD}}^{\text{CCSD(T)}}\big\vert_{\text{cc-pV[56]Z}} \; + \delta_{\text{CCSD(T)}}^{\text{CCSDT}}\big\vert_{\text{cc-pVTZ}} \; + \delta_{\text{CCSDT}}^{\text{CCSDT(Q)}}\big\vert_{\text{cc-pVDZ}}

    >>> energy('allen_focal_point')

    """
    lowername = name.lower()
    kwargs = kwargs_lower(kwargs)

    if not ('func_cbs' in kwargs):
        kwargs['func_cbs'] = energy

    # SCF
    if not ('scf_basis' in kwargs):
        kwargs['scf_basis'] = 'cc-pV[Q56]Z'
    if not ('scf_scheme' in kwargs):
        kwargs['scf_scheme'] = scf_xtpl_helgaker_3

    # delta MP2 - SCF
    if not ('corl_wfn' in kwargs):
        kwargs['corl_wfn'] = 'mp2'
    if not ('corl_basis' in kwargs):
        kwargs['corl_basis'] = 'cc-pV[56]Z'
    if not ('corl_scheme' in kwargs):
        kwargs['corl_scheme'] = corl_xtpl_helgaker_2

    # delta CCSD - MP2
    if not ('delta_wfn' in kwargs):
        kwargs['delta_wfn'] = 'mrccsd'
    if not ('delta_wfn_lesser' in kwargs):
        kwargs['delta_wfn_lesser'] = 'mp2'
    if not ('delta_basis' in kwargs):
        kwargs['delta_basis'] = 'cc-pV[56]Z'
    if not ('delta_scheme' in kwargs):
        kwargs['delta_scheme'] = corl_xtpl_helgaker_2

    # delta CCSD(T) - CCSD
    if not ('delta2_wfn' in kwargs):
        kwargs['delta2_wfn'] = 'mrccsd(t)'
    if not ('delta2_wfn_lesser' in kwargs):
        kwargs['delta2_wfn_lesser'] = 'mrccsd'
    if not ('delta2_basis' in kwargs):
        kwargs['delta2_basis'] = 'cc-pV[56]Z'
    if not ('delta2_scheme' in kwargs):
        kwargs['delta2_scheme'] = corl_xtpl_helgaker_2

    # delta CCSDT - CCSD(T)
    if not ('delta3_wfn' in kwargs):
        kwargs['delta3_wfn'] = 'mrccsdt'
    if not ('delta3_wfn_lesser' in kwargs):
        kwargs['delta3_wfn_lesser'] = 'mrccsd(t)'
    if not ('delta3_basis' in kwargs):
        kwargs['delta3_basis'] = 'cc-pVTZ'
    if not ('delta3_scheme' in kwargs):
        kwargs['delta3_scheme'] = highest_1

    # delta CCSDT(Q) - CCSDT
    if not ('delta4_wfn' in kwargs):
        kwargs['delta4_wfn'] = 'mrccsdt(q)'
    if not ('delta4_wfn_lesser' in kwargs):
        kwargs['delta4_wfn_lesser'] = 'mrccsdt'
    if not ('delta4_basis' in kwargs):
        kwargs['delta4_basis'] = 'cc-pVDZ'
    if not ('delta4_scheme' in kwargs):
        kwargs['delta4_scheme'] = highest_1

    return cbs(name, **kwargs)


#def run_mp2_5(name, **kwargs):
#    r"""Function that computes MP2.5 energy from results of a FNOCC
#    MP3 calculation.
#
#    .. math:: E_{total}^{\text{MP2.5}} = E_{total,\; \text{SCF}} \; + E_{corl,\; \text{MP2}} + E_{corl, \; \text{MP3}}
#
#    :PSI variables: 
#
#    .. hlist:: 
#       :columns: 1 
#     
#       * :psivar:`MP2.5 TOTAL ENERGY <MP2.5TOTALENERGY>` 
#       * :psivar:`MP2.5 CORRELATION ENERGY <MP2.5CORRELATIONENERGY>` 
#
#    >>> energy('mp2.5')
#
#    """
#    lowername = name.lower()
#    kwargs = kwargs_lower(kwargs)
#
#    # Run detci calculation and collect conventional quantities
#    energy('mp3', **kwargs)
#    e_scf = PsiMod.get_variable('SCF TOTAL ENERGY')
#    ce_mp2 = PsiMod.get_variable('MP2 CORRELATION ENERGY')
#    ce_mp3 = PsiMod.get_variable('MP3 CORRELATION ENERGY')
#    e_mp2 = e_scf + ce_mp2
#    e_mp3 = e_scf + ce_mp3
#
#    # Compute quantities particular to MP2.5
#    ce_mp25 = 0.5 * (ce_mp2 + ce_mp3)
#    e_mp25 = e_scf + ce_mp25
#    PsiMod.set_variable('MP2.5 CORRELATION ENERGY', ce_mp25)
#    PsiMod.set_variable('MP2.5 TOTAL ENERGY', e_mp25)
#    PsiMod.set_variable('CURRENT CORRELATION ENERGY', ce_mp25)
#    PsiMod.set_variable('CURRENT ENERGY', e_mp25)
#
#    # build string of title banner and print results
#    banners = ''
#    banners += """PsiMod.print_out('\\n')\n"""
#    banners += """banner(' MP2.5 ')\n"""
#    banners += """PsiMod.print_out('\\n')\n\n"""
#    exec(banners)
#
#    tables = ''
#    tables += """  SCF total energy:                        %16.8f\n""" % (e_scf)
#    tables += """  MP2 total energy:                        %16.8f\n""" % (e_mp2)
#    tables += """  MP2.5 total energy:                      %16.8f\n""" % (e_mp25)
#    tables += """  MP3 total energy:                        %16.8f\n\n""" % (e_mp3)
#    tables += """  MP2 correlation energy:                  %16.8f\n""" % (ce_mp2)
#    tables += """  MP2.5 correlation energy:                %16.8f\n""" % (ce_mp25)
#    tables += """  MP3 correlation energy:                  %16.8f\n""" % (ce_mp3)
#    PsiMod.print_out(tables)
#
#    return e_mp25


# A direct translation of a plugin input file into a function call. Function calls are the only
#     way to call plugins in sow/reap mode for db(), opt(), etc. This isn't best practices
#     but is an example of what to do for a more complicated procedure where different options 
#     are set for different qc steps.
#def run_plugin_omega(name, **kwargs):
#    r"""Function encoding sequence of PSI module and plugin calls, as well
#    as typical options, to access Rob Parrish's omega plugin.
#
#    >>> energy('plugin_omega')
#
#    """
#    lowername = name.lower()
#    kwargs = kwargs_lower(kwargs)
#
#    plugfile = PsiMod.Process.environment["PSIDATADIR"] + "/../tests/plugin_omega/plugin_omega.so"
#    PsiMod.plugin_load("%s" % (plugfile))
#
#    PsiMod.set_global_option('BASIS', 'AUG-CC-PVDZ')
#    PsiMod.set_global_option('DF_BASIS_SCF', 'AUG-CC-PVDZ-RI')
#    PsiMod.set_global_option('REFERENCE', 'UHF')
#    PsiMod.set_global_option('SCF_TYPE', 'DF')
#    energy('scf', **kwargs)
#
#    PsiMod.set_global_option('dft_functional', 'wB97')
#    PsiMod.set_global_option('dft_order_spherical', 25)
#    PsiMod.set_global_option('dft_num_radial', 35)
#    PsiMod.set_global_option('omega_procedure', 'ip')
#    PsiMod.set_global_option('maxiter', 50)
#    PsiMod.set_global_option('d_convergence', 5)
#    PsiMod.set_global_option('e_convergence', 7)
#    PsiMod.plugin("plugin_omega.so")
#
#    return PsiMod.get_variable('SCF TOTAL ENERGY')


# Integration with driver routines
#procedures['energy']['mp2.5'] = run_mp2_5
procedures['energy']['sherrillgroup_gold_standard'] = sherrillgroup_gold_standard
#procedures['energy']['plugin_omega'] = run_plugin_omega
