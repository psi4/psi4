#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
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
# @END LICENSE
#

from __future__ import print_function
from __future__ import absolute_import
import re
import math
from collections import defaultdict

from .exceptions import *
from . import qcformat
#import molpro_basissets
from . import options
from .pdict import PreservingDict


def harvest_output(outtext):
    """Function to separate portions of a Psi4 output file *outtext*.

    """
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

#    # Process NRE
#    mobj = re.search(r'^\s+' + r'(?:Nuclear Repulsion Energy =)' + r'\s+' + NUMBER + r'\s+hartrees\s*$',
#        outtext, re.MULTILINE)
#    if mobj:
#        print('matched nre')
#        psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)

    # Process HF  UNTESTED
    mobj = re.search(
        r'^\s+' + r'(?:Nuclear Repulsion Energy =)' + r'\s+' + NUMBER + r'\s+hartrees\s*' +
        r'(?:.*?)' +
        r'(?:Hartree-Fock SCF calculation)' +
        r'(?:.*?)' +
        r'^\s+\d+\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + 'Convergence criterion met' + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched hf')
        psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)
        psivar['HF TOTAL ENERGY'] = mobj.group(2)

    # Process DFT-D2  UNTESTED
#    mobj = re.search(
#        r'^\s+' + r'(?:Nuclear Repulsion Energy =)' + r'\s+' + NUMBER + r'\s+hartrees\s*' +
#        r'(?:.*?)' +
#        r'(?:HF-DFT SCF calculation)' +
#        r'(?:.*?)' +
#        r'^\s+' + r'(?:Empirical dispersion =)' + r'\s+' + NUMBER + r'\s+hartree\s*' +
#        r'(?:.*?)' +
#        r'^\s+\d+\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + 'Convergence criterion met' + r'\s*$',
#        outtext, re.MULTILINE | re.DOTALL)
#    if mobj:
#        print('matched dft-d2')
#        psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)
#        psivar['DISPERSION CORRECTION ENERGY'] = mobj.group(2)
#        psivar['DFT TOTAL ENERGY'] = mobj.group(3)
#        psivar['DFT FUNCTIONAL TOTAL ENERGY'] = mobj.group(3) - mboj.group(2)

    # Process DFT-D3  UNTESTED
#    mobj = re.search(
#        r'(?:grimme3)' + r'\s*' +
#        r'(?:.*?)' +
#        r'^\s+' + r'(?:Nuclear Repulsion Energy =)' + r'\s+' + NUMBER + r'\s+hartrees\s*' +
#        r'(?:.*?)' +
#        r'(?:HF-DFT SCF calculation)' +
#        r'(?:.*?)' +
#        r'^\s+\d+\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + 'Convergence criterion met' + r'\s*$',
#        outtext, re.MULTILINE | re.DOTALL)
#    if mobj:
#        print('matched dft-d3')
#        psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)
#        psivar['DISPERSION CORRECTION ENERGY'] = None
#        psivar['DFT TOTAL ENERGY'] = mobj.group(2)
#        psivar['DFT FUNCTIONAL TOTAL ENERGY'] = None

# /^((?!PART).)*$/

    # Process DFT no-D or internal-D
    mobj = re.search(
#        r'((?!grimme3).)*' + r'\s*' +  # severe negative performance impact
#        r'(?:.*?)' +
        r'^\s+' + r'(?:Nuclear Repulsion Energy =)' + r'\s+' + NUMBER + r'\s+hartrees\s*' +
        r'(?:.*?)' +
        r'(?:HF-DFT SCF calculation)' +
        r'(?:.*?)' +
        r'^\s+\d+\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + 'Convergence criterion met' + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    if mobj:
        print('matched dft')
        psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)
        #psivar['DFT TOTAL ENERGY'] = mobj.group(2)
        psivar['DFT FUNCTIONAL TOTAL ENERGY'] = mobj.group(2)
        # with negative lookahead
        #psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(2)
        #psivar['DFT TOTAL ENERGY'] = mobj.group(3)
        #psivar['DFT FUNCTIONAL TOTAL ENERGY'] = mobj.group(3)

    # Process DHDFT no-D or internal-D
    mobj = re.search(
        # negative grimme3 lookahead goes here
        #r'^\s+' + r'(?:Nuclear Repulsion Energy =)' + r'\s+' + NUMBER + r'\s+hartrees\s*' +
        #r'(?:.*?)' +
        r'(?:HF-DFT SCF calculation)' +
        r'(?:.*?)' +
        #r'^\s+\d+\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + 'Convergence criterion met' + r'\s*' +
        #r'(?:.*?)' +
        # need a not "Hartree-Fock SCF calculation" here so DFT @@@ MP2 not caught?
        r'^\s*' + r'(?:Total  (?:RI)?MP2   correlation energy =)' + r'\s+' + NUMBER + r'\s+' + r'au' + r'\s*' +
        r'^\s+' + r'(?:(?:RI)?MP2         total energy =)' + r'\s+' + NUMBER + r'\s+' + r'au' + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    if mobj:
        print('matched dhdft')
        #psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)
        #psivar['DFT TOTAL ENERGY'] = mobj.group(2)
        #psivar['DFT FUNCTIONAL TOTAL ENERGY'] = mobj.group(2)
        psivar['DOUBLE-HYBRID CORRECTION ENERGY'] = mobj.group(1)

    # Process MP2
    mobj = re.search(
        r'(?:Hartree-Fock SCF calculation)' +
        r'(?:.*?)' +
        #r'^\s+\d+\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + 'Convergence criterion met' + r'\s*' +
        #r'(?:.*?)' +
        # need a not "Hartree-Fock SCF calculation" here so DFT @@@ MP2 not caught?
        r'^\s*' + r'(?:Total  RIMP2   correlation energy =)' + r'\s+' + NUMBER + r'\s+' + r'au' + r'\s*' +
        r'^\s+' + r'(?:RIMP2         total energy =)' + r'\s+' + NUMBER + r'\s+' + r'au' + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL | re.IGNORECASE)
    if mobj:
        print('matched mp2')
        #psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)
        #psivar['DFT TOTAL ENERGY'] = mobj.group(2)
        #psivar['DFT FUNCTIONAL TOTAL ENERGY'] = mobj.group(2)
        psivar['MP2 CORRELATION ENERGY'] = mobj.group(1)
        #psivar['DOUBLE-HYBRID CORRECTION ENERGY'] = mobj.group(1)
        print(psivar)

# TODO: need to split on 'Q-Chem begins' or 'Quantum Leap' or something

#    # Process DFT no-D or internal-D WORKS BUT LOOKAHEAD VERY SLOW
#    mobj = re.search(
#        r'((?!grimme3).)*' + r'\s*' +  # severe negative performance impact
#           TODO note neg lookahead insufficient since option could be negated
#        r'(?:.*?)' +
#        r'^\s+' + r'(?:Nuclear Repulsion Energy =)' + r'\s+' + NUMBER + r'\s+hartrees\s*' +
#        r'(?:.*?)' +
#        r'(?:HF-DFT SCF calculation)' +
#        r'(?:.*?)' +
#        r'^\s+\d+\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + 'Convergence criterion met' + r'\s*$',
#        outtext, re.MULTILINE | re.DOTALL | re.IGNORECASE)
#    if mobj:
#        print('matched dft')
#        psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(2)
#        psivar['DFT TOTAL ENERGY'] = mobj.group(3)
#        psivar['DFT FUNCTIONAL TOTAL ENERGY'] = mobj.group(3)

#    # Process PsiVariables
#    mobj = re.search(r'^(?:  Variable Map:)\s*' +
#        r'^\s*(?:-+)\s*' +
#        r'^(.*?)' +
#        r'^(?:\s*?)$',
#        outtext, re.MULTILINE | re.DOTALL)
#
#    if mobj:
#        for pv in mobj.group(1).split('\n'):
#            submobj = re.search(r'^\s+' + r'"(.+?)"' + r'\s+=>\s+' + NUMBER + r'\s*$', pv)
#            if submobj:
#                psivar['%s' % (submobj.group(1))] = submobj.group(2)

    # Process Completion
    mobj = re.search(r'Thank you very much for using Q-Chem.  Have a nice day.',
        outtext, re.MULTILINE)
    if mobj:
        psivar['SUCCESS'] = True

    return psivar, psivar_coord, psivar_grad


def muster_memory(mem):
    """Transform input *mem* in MB into psi4-type options.

    """
    text = ''

    # prepare memory keywords to be set as c-side keywords
    options = defaultdict(lambda: defaultdict(dict))
    options['QCHEM']['QCHEM_MEM_TOTAL']['value'] = int(mem)
    #options['QCHEM']['QCHEM_CC_MEMORY']['value'] = int(mem)
    #options['QCHEM']['QCHEM_MEM_STATIC']['value'] = int(mem)

    for item in options['QCHEM']:
        options['QCHEM'][item]['clobber'] = True
    return text, options


def muster_basis(bas):
    """Transform input *mem* in MB into psi4-type options.

    """
    text = ''

    # prepare memory keywords to be set as c-side keywords
    options = defaultdict(lambda: defaultdict(dict))
    options['QCHEM']['QCHEM_BASIS']['value'] = bas

    for item in options['QCHEM']:
        options['QCHEM'][item]['clobber'] = True
    return text, options


class Infile(qcformat.InputFormat2):

    def __init__(self, mem, mol, mtd, der, opt):
        qcformat.InputFormat2.__init__(self, mem, mol, mtd, der, opt)

#        #print self.method, self.molecule.nactive_fragments()
#        if 'sapt' in self.method and self.molecule.nactive_fragments() != 2:
#            raise FragmentCountError("""Requested molecule has %d, not 2, fragments.""" % (self.molecule.nactive_fragments()))
#
##        # memory in MB --> MW
##        self.memory = int(math.ceil(mem / 8.0))
##        # auxiliary basis sets
##        [self.unaugbasis, self.augbasis, self.auxbasis] = self.corresponding_aux_basis()

    def format_infile_string(self):
        """

        """
        # Handle memory and comment
        cmtcmd = """$comment\n%s\n$end\n\n""" % (self.molecule.tagline)
        memcmd, memkw = muster_memory(self.memory)

        # Handle molecule and basis set
        molcmd, molkw = self.molecule.format_molecule_for_qchem(mixedbas=False)
        # TODO mixedbas=True once handling basis sets

        # not translating basis at present
        _bascmd, baskw = muster_basis(self.basis)

        # format global convergence directions
        _cdscmd, cdskw = muster_cdsgroup_options()

        # Handle calc type and quantum chemical method
        mdccmd, mdckw = procedures['energy'][self.method](self.method, self.dertype)

## make options from imdb only user options (currently non-existent). set basis and castup from here.
        # Handle driver vs input/default keyword reconciliation
        userkw = self.options  # p4util.prepare_options_for_modules()
        userkw = options.reconcile_options2(userkw, memkw)
        userkw = options.reconcile_options2(userkw, molkw)
        userkw = options.reconcile_options2(userkw, baskw)
        #userkw = qcdb.options.reconcile_options(userkw, psikw)
        userkw = options.reconcile_options2(userkw, cdskw)
        userkw = options.reconcile_options2(userkw, mdckw)

        # Handle conversion of psi4 keyword structure into cfour format
        optcmd = options.prepare_options_for_qchem(userkw)

        # Handle text to be passed untouched to psi4
        litcmd = ''

        # Assemble infile pieces
        return cmtcmd + memcmd + molcmd + optcmd + mdccmd + litcmd

#'hf'
#'df-hf'
#'b3lyp'
#'blyp'
#'bp86'
#'fno-ccsd(t)'
#'df-ccsd(t)'
#'fno-df-ccsd(t)'
#'df-b97-d'
#'df-b97-d3'
#'pbe0-2'
#'dsd-pbep86'
#'wb97x-2'
#'DLdf+d'
#'DLdf+d09'
#'df-b3lyp'
#'df-b3lyp-d'
#'df-b3lyp-d3'
#'df-wb97x-d'


def muster_cdsgroup_options():
    text = ''
    options = defaultdict(lambda: defaultdict(dict))
#    options['GLOBALS']['E_CONVERGENCE']['value'] = 8
#    options['SCF']['GUESS']['value'] = 'sad'
#    options['SCF']['MAXITER']['value'] = 200
    options['QCHEM']['QCHEM_MEM_STATIC']['value'] = 512
    options['QCHEM']['QCHEM_XC_GRID']['value'] = '000100000302'
    options['QCHEM']['QCHEM_THRESH']['value'] = 12
    options['QCHEM']['QCHEM_SCF_CONVERGENCE']['value'] = 7
    #options['QCHEM']['QCHEM_INTEGRALS_BUFFER']['value'] = 512
    options['QCHEM']['QCHEM_MAX_SCF_CYCLES']['value'] = 200

    options['QCHEM']['QCHEM_SYM_IGNORE']['value'] = True
    options['QCHEM']['QCHEM_SYMMETRY']['value'] = False
    options['QCHEM']['QCHEM_INTEGRALS_BUFFER']['value'] = 512

    return text, options


def muster_modelchem(name, dertype):
    """Transform calculation method *name* and derivative level *dertype*
    into options for cfour. While deliberately requested pieces,
    generally |cfour__cfour_deriv_level| and |cfour__cfour_calc_level|,
    are set to complain if contradicted ('clobber' set to True), other
    'recommended' settings, like |cfour__cfour_cc_program|, can be
    countermanded by keywords in input file ('clobber' set to False).
    Occasionally, want these pieces to actually overcome keywords in
    input file ('superclobber' set to True).

    """
    text = ''
    lowername = name.lower()
    options = defaultdict(lambda: defaultdict(dict))

    if dertype == 0:
        options['QCHEM']['QCHEM_JOBTYPE']['value'] = 'SP'
#        text += """energy('"""
    else:
        raise ValidationError("""Requested Psi4 dertype %d is not available.""" % (dertype))

    if lowername == 'wb97x-v':
        options['QCHEM']['QCHEM_EXCHANGE']['value'] = 'omegaB97X-V'

#        text += """mp2')\n\n"""
#
#    elif lowername == 'df-mp2':
#        options['GLOBALS']['FREEZE_CORE']['value'] = True
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        options['MP2']['MP2_TYPE']['value'] = 'df'
#        text += """mp2')\n\n"""
#
#    elif lowername == 'sapt0':
#        options['GLOBALS']['FREEZE_CORE']['value'] = True
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        text += """sapt0')\n\n"""
#
#    elif lowername == 'sapt2+':
#        options['GLOBALS']['FREEZE_CORE']['value'] = True
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        options['SAPT']['NAT_ORBS_T2']['value'] = True
#        options['SAPT']['NAT_ORBS_T3']['value'] = True
#        options['SAPT']['NAT_ORBS_V4']['value'] = True
#        options['SAPT']['OCC_TOLERANCE']['value'] = 1.0e-6
#        text += """sapt2+')\n\n"""
#
#    elif lowername == 'sapt2+(3)':
#        options['GLOBALS']['FREEZE_CORE']['value'] = True
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        options['SAPT']['NAT_ORBS_T2']['value'] = True
#        options['SAPT']['NAT_ORBS_T3']['value'] = True
#        options['SAPT']['NAT_ORBS_V4']['value'] = True
#        options['SAPT']['OCC_TOLERANCE']['value'] = 1.0e-6
#        text += """sapt2+(3)')\n\n"""
#
#    elif lowername == 'sapt2+3(ccd)':
#        options['GLOBALS']['FREEZE_CORE']['value'] = True
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        options['SAPT']['NAT_ORBS_T2']['value'] = True
#        options['SAPT']['NAT_ORBS_T3']['value'] = True
#        options['SAPT']['NAT_ORBS_V4']['value'] = True
#        options['SAPT']['OCC_TOLERANCE']['value'] = 1.0e-6
#        options['SAPT']['DO_MBPT_DISP']['value'] = True
#        text += """sapt2+3(ccd)')\n\n"""
#
#    elif lowername == 'df-b97-d3':
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        options['SCF']['DFT_SPHERICAL_POINTS']['value'] = 302
#        options['SCF']['DFT_RADIAL_POINTS']['value'] = 100
#        text += """b97-d3')\n\n"""
#
#    elif lowername == 'df-wb97x-d':
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        options['SCF']['DFT_SPHERICAL_POINTS']['value'] = 302
#        options['SCF']['DFT_RADIAL_POINTS']['value'] = 100
#        text += """wb97x-d')\n\n"""
#
#    elif lowername == 'df-b3lyp-d3':
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        options['SCF']['DFT_SPHERICAL_POINTS']['value'] = 302
#        options['SCF']['DFT_RADIAL_POINTS']['value'] = 100
#        text += """b3lyp-d3')\n\n"""
#
#    elif lowername == 'dfdf-b2plyp-d3':
#        options['GLOBALS']['FREEZE_CORE']['value'] = True
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        options['DFMP2']['MP2_TYPE']['value'] = 'df'
#        options['SCF']['DFT_SPHERICAL_POINTS']['value'] = 302
#        options['SCF']['DFT_RADIAL_POINTS']['value'] = 100
#        text += """b2plyp-d3')\n\n"""
#
#    elif lowername == 'df-wpbe':
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        options['SCF']['DFT_SPHERICAL_POINTS']['value'] = 302
#        options['SCF']['DFT_RADIAL_POINTS']['value'] = 100
#        text += """wpbe')\n\n"""
#
#    elif lowername == 'ccsd-polarizability':
#        options['GLOBALS']['FREEZE_CORE']['value'] = True
#        text = """property('ccsd', properties=['polarizability'])\n\n"""
#
#    elif lowername == 'mrccsdt(q)':
#        options['SCF']['SCF_TYPE']['value'] = 'pk'
#        options['GLOBALS']['FREEZE_CORE']['value'] = True
#        options['GLOBALS']['NAT_ORBS']['value'] = True  # needed by mrcc but not recognized by mrcc
#        options['FNOCC']['OCC_TOLERANCE']['value'] = 6
#        text += """mrccsdt(q)')\n\n"""
#
#    elif lowername == 'c4-ccsdt(q)':
#        options['CFOUR']['CFOUR_SCF_CONV']['value'] = 11
#        options['CFOUR']['CFOUR_CC_CONV']['value'] = 10
#        options['CFOUR']['CFOUR_FROZEN_CORE']['value'] = True
#        text += """c4-ccsdt(q)')\n\n"""
#
#    elif lowername == 'df-m05-2x':
#        options['SCF']['SCF_TYPE']['value'] = 'df'
#        options['SCF']['DFT_SPHERICAL_POINTS']['value'] = 302
#        options['SCF']['DFT_RADIAL_POINTS']['value'] = 100
#        text += """m05-2x')\n\n"""

    else:
        raise ValidationError("""Requested Psi4 computational methods %d is not available.""" % (lowername))

#    # Set clobbering
#    if 'CFOUR_DERIV_LEVEL' in options['CFOUR']:
#        options['CFOUR']['CFOUR_DERIV_LEVEL']['clobber'] = True
#        options['CFOUR']['CFOUR_DERIV_LEVEL']['superclobber'] = True
#    if 'CFOUR_CALC_LEVEL' in options['CFOUR']:
#        options['CFOUR']['CFOUR_CALC_LEVEL']['clobber'] = True
#        options['CFOUR']['CFOUR_CALC_LEVEL']['superclobber'] = True
#    if 'CFOUR_CC_PROGRAM' in options['CFOUR']:
#        options['CFOUR']['CFOUR_CC_PROGRAM']['clobber'] = False

    return text, options

procedures = {
    'energy': {
        'wb97x-v'       : muster_modelchem,
    }
}

qcmtdIN = procedures['energy']


def psi4_list():
    """Return an array of Psi4 methods with energies.

    """
    return sorted(procedures['energy'].keys())
