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

import re
import math
import collections

from exceptions import *
import qcformat
#import molpro_basissets
import options
from pdict import PreservingDict


def harvest_output(outtext):
    """Function to separate portions of a PSI4 output file *outtext*.

    """
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    # Process PsiVariables
    mobj = re.search(r'^(?:  Variable Map:)\s*' +
        r'^\s*(?:-+)\s*' +
        r'^(.*?)' +
        r'^(?:\s*?)$',
        outtext, re.MULTILINE | re.DOTALL)

    if mobj:
        for pv in mobj.group(1).split('\n'):
            submobj = re.search(r'^\s+' + r'"(.+?)"' + r'\s+=>\s+' + NUMBER + r'\s*$', pv)
            if submobj:
                psivar['%s' % (submobj.group(1))] = submobj.group(2)

    # Process Completion
    mobj = re.search(r'PSI4 exiting successfully. Buy a developer a beer!',
        outtext, re.MULTILINE)
    if mobj:
        psivar['SUCCESS'] = True

    return psivar, psivar_coord, psivar_grad


class Infile(qcformat.InputFormat2):

    def __init__(self, mem, mol, mtd, der, opt):
        qcformat.InputFormat2.__init__(self, mem, mol, mtd, der, opt)

        print self.method, self.molecule.nactive_fragments()
        if 'sapt' in self.method and self.molecule.nactive_fragments() != 2:
            raise FragmentCountError("""Requested molecule has %d, not 2, fragments.""" % (self.molecule.nactive_fragments()))

#        # memory in MB --> MW
#        self.memory = int(math.ceil(mem / 8.0))
#        # auxiliary basis sets
#        [self.unaugbasis, self.augbasis, self.auxbasis] = self.corresponding_aux_basis()

    def format_infile_string(self):
        """

        """
        # Handle memory and comment
        memcmd, memkw = """# %s\n\nmemory %d mb\n\n""" % (self.molecule.tagline, self.memory), {}

        # Handle molecule and basis set
        molcmd, molkw = self.molecule.format_molecule_for_psi4(), {}

        # format global convergence directions
#        text += self.format_global_parameters()
        _cdscmd, cdskw = muster_cdsgroup_options()

        # Handle calc type and quantum chemical method
        mdccmd, mdckw = procedures['energy'][self.method](self.method, self.dertype)

#        # format options
#        optcmd = qcdb.options.prepare_options_for_psi4(mdckw)

# make options from imdb only user options (currently non-existent). set basis and castup from here.
        # Handle driver vs input/default keyword reconciliation
        userkw = self.options
#        userkw = p4util.prepare_options_for_modules()
        #userkw = qcdb.options.reconcile_options(userkw, memkw)
        #userkw = qcdb.options.reconcile_options(userkw, molkw)
        #userkw = qcdb.options.reconcile_options(userkw, baskw)
        #userkw = qcdb.options.reconcile_options(userkw, psikw)
        userkw = options.reconcile_options2(userkw, cdskw)
        userkw = options.reconcile_options2(userkw, mdckw)

        # Handle conversion of psi4 keyword structure into cfour format
        optcmd = options.prepare_options_for_psi4(userkw)

        # Handle text to be passed untouched to psi4
        litcmd = """\nprint_variables()\n\n"""

        # Assemble infile pieces
        return memcmd + molcmd + optcmd + mdccmd + litcmd

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
    options = collections.defaultdict(lambda: collections.defaultdict(dict))
    options['GLOBALS']['E_CONVERGENCE']['value'] = 8
    options['SCF']['GUESS']['value'] = 'sad'
    options['SCF']['MAXITER']['value'] = 200

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
    options = collections.defaultdict(lambda: collections.defaultdict(dict))

    if dertype == 0:
        text += """energy('"""
    else:
        raise ValidationError("""Requested Psi4 dertype %d is not available.""" % (dertype))

    if lowername == 'mp2':
        options['GLOBALS']['FREEZE_CORE']['value'] = True
        options['SCF']['SCF_TYPE']['value'] = 'direct'
        options['MP2']['MP2_TYPE']['value'] = 'conv'
        text += """mp2')\n\n"""

    elif lowername == 'df-mp2':
        options['GLOBALS']['FREEZE_CORE']['value'] = True
        options['SCF']['SCF_TYPE']['value'] = 'df'
        options['MP2']['MP2_TYPE']['value'] = 'df'
        text += """mp2')\n\n"""

    elif lowername == 'sapt0':
        options['GLOBALS']['FREEZE_CORE']['value'] = True
        options['SCF']['SCF_TYPE']['value'] = 'df'
        text += """sapt0')\n\n"""

    elif lowername == 'sapt2+':
        options['GLOBALS']['FREEZE_CORE']['value'] = True
        options['SCF']['SCF_TYPE']['value'] = 'df'
        options['SAPT']['NAT_ORBS_T2']['value'] = True
        options['SAPT']['NAT_ORBS_T3']['value'] = True
        options['SAPT']['NAT_ORBS_V4']['value'] = True
        options['SAPT']['OCC_TOLERANCE']['value'] = 1.0e-6
        text += """sapt2+')\n\n"""

    elif lowername == 'sapt2+(3)':
        options['GLOBALS']['FREEZE_CORE']['value'] = True
        options['SCF']['SCF_TYPE']['value'] = 'df'
        options['SAPT']['NAT_ORBS_T2']['value'] = True
        options['SAPT']['NAT_ORBS_T3']['value'] = True
        options['SAPT']['NAT_ORBS_V4']['value'] = True
        options['SAPT']['OCC_TOLERANCE']['value'] = 1.0e-6
        text += """sapt2+(3)')\n\n"""

    elif lowername == 'sapt2+3(ccd)':
        options['GLOBALS']['FREEZE_CORE']['value'] = True
        options['SCF']['SCF_TYPE']['value'] = 'df'
        options['SAPT']['NAT_ORBS_T2']['value'] = True
        options['SAPT']['NAT_ORBS_T3']['value'] = True
        options['SAPT']['NAT_ORBS_V4']['value'] = True
        options['SAPT']['OCC_TOLERANCE']['value'] = 1.0e-6
        options['SAPT']['DO_MBPT_DISP']['value'] = True
        text += """sapt2+3(ccd)')\n\n"""

    elif lowername == 'df-b97-d3':
        options['SCF']['SCF_TYPE']['value'] = 'df'
        options['SCF']['DFT_SPHERICAL_POINTS']['value'] = 302
        options['SCF']['DFT_RADIAL_POINTS']['value'] = 100
        text += """b97-d3')\n\n"""

    elif lowername == 'df-wb97x-d':
        options['SCF']['SCF_TYPE']['value'] = 'df'
        options['SCF']['DFT_SPHERICAL_POINTS']['value'] = 302
        options['SCF']['DFT_RADIAL_POINTS']['value'] = 100
        text += """wb97x-d')\n\n"""

    elif lowername == 'df-b3lyp-d3':
        options['SCF']['SCF_TYPE']['value'] = 'df'
        options['SCF']['DFT_SPHERICAL_POINTS']['value'] = 302
        options['SCF']['DFT_RADIAL_POINTS']['value'] = 100
        text += """b3lyp-d3')\n\n"""

    elif lowername == 'ccsd-polarizability':
        options['GLOBALS']['FREEZE_CORE']['value'] = True
        text = """property('ccsd', properties=['polarizability'])\n\n"""

    else:
        raise ValidationError("""Requested Cfour computational methods %d is not available.""" % (lowername))

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
        'df-b97-d3'     : muster_modelchem,
        'df-wb97x-d'    : muster_modelchem,
        'df-b3lyp-d3'   : muster_modelchem,
        'mp2'           : muster_modelchem,
        'df-mp2'        : muster_modelchem,
        'sapt0'         : muster_modelchem,
        'sapt2+'        : muster_modelchem,
        'sapt2+(3)'     : muster_modelchem,
        'sapt2+3(ccd)'  : muster_modelchem,
        'ccsd-polarizability'  : muster_modelchem,
    }
}

qcmtdIN = procedures['energy']


def psi4_list():
    """Return an array of Psi4 methods with energies.

    """
    return sorted(procedures['energy'].keys())
