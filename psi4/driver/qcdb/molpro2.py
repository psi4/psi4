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

from __future__ import absolute_import
from __future__ import print_function
import re
import math
from decimal import Decimal
from collections import defaultdict
from .exceptions import *
from .pdict import PreservingDict
from . import qcformat
from . import molpro_basissets
from . import options


def harvest_output(outtext):
    """Function to read MRCC output file *outtext* and parse important
    quantum chemical information from it in

    """
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"


    # <<< Process NRE >>>
    mobj = re.search(
        r'^\s*' + r'(?:NUCLEAR REPULSION ENERGY)' + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        #print('matched nrc')
        psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)

    # <<< Process SCF >>>

    #mobj = re.search(
    #    r'^\s*' + r'(?:Energy of reference determinant (?:\[au\]|/au/):)' + r'\s+' + NUMBER + r'\s*$',
    #    outtext, re.MULTILINE)
    #if mobj:
    #    print('matched scf')
    #    psivar['SCF TOTAL ENERGY'] = mobj.group(1)

    # <<< Process MP2 >>>

    mobj = re.search(
        r'^\s*' + r'Reference energy[:]?\s+' + NUMBER + r'\s*' +
        r'^\s*' + r'MP2 singlet pair energy[:]?\s+' + NUMBER + r'\s*' +
        r'^\s*' + r'MP2 triplet pair energy[:]?\s+' + NUMBER + r'\s*' +
        r'^\s*' + r'MP2 correlation energy[:]?\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        #print('matched mp2')
        psivar['HF TOTAL ENERGY'] = mobj.group(1)
        psivar['MP2 CORRELATION ENERGY'] = mobj.group(4)
        psivar['MP2 TOTAL ENERGY'] = Decimal(mobj.group(1)) + Decimal(mobj.group(4))
        psivar['MP2 SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(3)) * \
            Decimal(2) / Decimal(3)
        psivar['MP2 OPPOSITE-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(4)) - \
            psivar['MP2 SAME-SPIN CORRELATION ENERGY']

    # <<< Process SAPT-like >>>

    mobj = re.search(
        #r'^\s+' + r'E1pol\s+' + NUMBER + r'\s+\(\s*' + NUMBER + r')\s+' + NUMBER + r'\s+' + NUMBER + '\s*'
        #r'^\s+' + r'E1exch\s+' + NUMBER + r'\s+\(\s*' + NUMBER + r')\s+' + NUMBER + r'\s+' + NUMBER + '\s*'
        #r'^\s+' + r'E1exch\(S2\)\s+' + NUMBER + r'\s+\(\s*' + NUMBER + r')\s+' + NUMBER + r'\s+' + NUMBER + '\s*'
        #r'^\s+' + r'E2ind\(unc\)\s+' + NUMBER + r'\s+\(\s*' + NUMBER + r')\s+' + NUMBER + r'\s+' + NUMBER + '\s*'
        #r'^\s+' + r'E2ind\s+' + NUMBER + r'\s+\(\s*' + NUMBER + r')\s+' + NUMBER + r'\s+' + NUMBER + '\s*'
        #r'^\s+' + r'E2ind-exch\s+' + NUMBER + r'\s+\(\s*' + NUMBER + r')\s+' + NUMBER + r'\s+' + NUMBER + '\s*'
        #r'^\s+' + r'E2disp\(unc\)\s+' + NUMBER + r'\s+\(\s*' + NUMBER + r')\s+' + NUMBER + r'\s+' + NUMBER + '\s*'
#        r'^\s+' + r'E2disp\s+' + NUMBER + r'\s+\(\s*' + NUMBER + r'\)\s+' + NUMBER + r'\s+' + NUMBER + '\s*',
        r'^\s+' + r'E2disp\s+' + NUMBER + r'.*$',
        #r'^\s+' + r'E2disp-exch\(unc\)\s+' + NUMBER + r'\s+\(\s*' + NUMBER + r')\s+' + NUMBER + r'\s+' + NUMBER + '\s*'
        #r'^\s+' + r'E2disp-exc\s+' + NUMBER + r'\s+\(\s*' + NUMBER + r')\s+' + NUMBER + r'\s+' + NUMBER + '\s*'
        outtext, re.MULTILINE)
    if mobj:
        #print('matched sapt-like')
        psivar['MP2C DISP20 ENERGY'] = Decimal(mobj.group(1)) / Decimal(1000)

    # <<< Process SCF-F12 >>>

    mobj = re.search(
        r'^\s+' + r'CABS-singles contribution of\s+' + NUMBER + r'\s+patched into reference energy.\s*' +
        r'^\s+' + r'New reference energy\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        #print('matched scff12')
        psivar['SCF TOTAL ENERGY'] = Decimal(mobj.group(2)) - Decimal(mobj.group(1))
        psivar['HF-CABS TOTAL ENERGY'] = mobj.group(2)

    # <<< Process MP2-F12 >>>

# DF-MP2-F12 correlation energies:
# --------------------------------
# Approx.                                    Singlet             Triplet             Ecorr            Total Energy
# DF-MP2                                -0.261035854033     -0.140514056591     -0.401549910624   -112.843952380305
# DF-MP2-F12/3*C(DX,FIX)                -0.367224875485     -0.163178266500     -0.530403141984   -112.972805611666
# DF-MP2-F12/3*C(FIX)                   -0.358294348708     -0.164988061549     -0.523282410258   -112.965684879939
# DF-MP2-F12/3C(FIX)                    -0.357375628783     -0.165176490386     -0.522552119169   -112.964954588851
#
# DF-MP2-F12 correlation energies:
# ================================
# Approx.                                    Singlet             Triplet             Ecorr            Total Energy
# DF-MP2                                -0.357960885582     -0.185676627667     -0.543637513249   -132.841755020796
# DF-MP2-F12/3*C(DX,FIX)                -0.381816069559     -0.188149510095     -0.569965579654   -132.868083087202
# DF-MP2-F12/3*C(FIX)                   -0.379285470419     -0.187468208608     -0.566753679027   -132.864871186575
# DF-MP2-F12/3C(FIX)                    -0.379246010149     -0.187531433611     -0.566777443760   -132.864894951307

    mobj = re.search(
        r'^\s*' + r'DF-MP2-F12 correlation energies:\s*' +
        r'^\s*(?:[=-]+)\s*' +
        r'^\s+' + r'Approx.\s+Singlet\s+Triplet\s+Ecorr\s+Total Energy\s*' + 
        r'^\s+' + r'DF-MP2\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' + 
        r'^\s+' + r'DF-MP2-F12/3\*C\(DX,FIX\)\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'DF-MP2-F12/3\*C\(FIX\)\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'DF-MP2-F12/3C\(FIX\)\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        #print('matched mp2f12')
        psivar['MP2 CORRELATION ENERGY'] = mobj.group(3)
        psivar['MP2 SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(2)) * \
            Decimal(2) / Decimal(3)
        psivar['MP2 OPPOSITE-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(3)) - \
            psivar['MP2 SAME-SPIN CORRELATION ENERGY']
        psivar['MP2 TOTAL ENERGY'] = mobj.group(4)
        psivar['MP2-F12 CORRELATION ENERGY'] = mobj.group(15)
        psivar['MP2-F12 SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(14)) * \
            Decimal(2) / Decimal(3)
        psivar['MP2-F12 OPPOSITE-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(15)) - \
            psivar['MP2-F12 SAME-SPIN CORRELATION ENERGY']
        psivar['MP2-F12 TOTAL ENERGY'] = mobj.group(16)

    # <<< Process CC >>>

    mobj = re.search(
        r'^\s*' + r'CCSD triplet pair energy\s+' + NUMBER + '\s*' +
        r'^\s*' + r'CCSD correlation energy\s+' + NUMBER + '\s*' +
        r'^\s*' + r'Triples \(T\) contribution\s+' + NUMBER + '\s*$',
        outtext, re.MULTILINE)
    if mobj:
        #print('matched ccsd(t)')
        psivar['CCSD CORRELATION ENERGY'] = mobj.group(2)
        psivar['CCSD SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(1)) * \
            Decimal(2) / Decimal(3)
        psivar['CCSD OPPOSITE-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(2)) - \
            psivar['CCSD SAME-SPIN CORRELATION ENERGY']
        psivar['CCSD TOTAL ENERGY'] = Decimal(mobj.group(2)) + psivar['HF TOTAL ENERGY']
        psivar['(T) CORRECTION ENERGY'] = mobj.group(3)
        psivar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(2)) + Decimal(mobj.group(3))
        psivar['CCSD(T) TOTAL ENERGY'] = psivar['CCSD(T) CORRELATION ENERGY'] + psivar['HF TOTAL ENERGY']

    # <<< Process CC-F12 >>>

    mobj = re.search(
        r'^\s*' + r'CCSD-F12a triplet pair energy\s+' + NUMBER + '\s*' +
        r'^\s*' + r'CCSD-F12a correlation energy\s+' + NUMBER + '\s*' +
        r'^\s*' + r'Triples \(T\) contribution\s+' + NUMBER + '\s*$',
        outtext, re.MULTILINE)
    if mobj:
        #print('matched ccsd(t)-f12a')
        psivar['CCSD-F12A CORRELATION ENERGY'] = mobj.group(2)
        psivar['CCSD-F12A SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(1)) * \
            Decimal(2) / Decimal(3)
        psivar['CCSD-F12A OPPOSITE-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(2)) - \
            psivar['CCSD-F12A SAME-SPIN CORRELATION ENERGY']
        psivar['CCSD-F12A TOTAL ENERGY'] = Decimal(mobj.group(2)) + psivar['HF-CABS TOTAL ENERGY']
        psivar['(T)-F12AB CORRECTION ENERGY'] = mobj.group(3)
        psivar['CCSD(T)-F12A CORRELATION ENERGY'] = Decimal(mobj.group(2)) + Decimal(mobj.group(3))
        psivar['CCSD(T)-F12A TOTAL ENERGY'] = psivar['CCSD(T)-F12A CORRELATION ENERGY'] + psivar['HF-CABS TOTAL ENERGY']
        psivar['(T*)-F12AB CORRECTION ENERGY'] = Decimal(mobj.group(3)) * \
            psivar['MP2-F12 CORRELATION ENERGY'] / psivar['MP2 CORRELATION ENERGY']
        psivar['CCSD(T*)-F12A CORRELATION ENERGY'] = Decimal(mobj.group(2)) + psivar['(T*)-F12AB CORRECTION ENERGY']
        psivar['CCSD(T*)-F12A TOTAL ENERGY'] = psivar['CCSD(T*)-F12A CORRELATION ENERGY'] + psivar['HF-CABS TOTAL ENERGY']

    mobj = re.search(
        r'^\s*' + r'CCSD-F12b triplet pair energy\s+' + NUMBER + '\s*' +
        r'^\s*' + r'CCSD-F12b correlation energy\s+' + NUMBER + '\s*' +
        r'^\s*' + r'Triples \(T\) contribution\s+' + NUMBER + '\s*$',
        outtext, re.MULTILINE)
    if mobj:
        #print('matched ccsd(t)-f12b')
        psivar['CCSD-F12B CORRELATION ENERGY'] = mobj.group(2)
        psivar['CCSD-F12B SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(1)) * \
            Decimal(2) / Decimal(3)
        psivar['CCSD-F12B OPPOSITE-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(2)) - \
            psivar['CCSD-F12B SAME-SPIN CORRELATION ENERGY']
        psivar['CCSD-F12B TOTAL ENERGY'] = Decimal(mobj.group(2)) + psivar['HF-CABS TOTAL ENERGY']
        psivar['(T)-F12AB CORRECTION ENERGY'] = mobj.group(3)
        psivar['CCSD(T)-F12B CORRELATION ENERGY'] = Decimal(mobj.group(2)) + Decimal(mobj.group(3))
        psivar['CCSD(T)-F12B TOTAL ENERGY'] = psivar['CCSD(T)-F12B CORRELATION ENERGY'] + psivar['HF-CABS TOTAL ENERGY']
        psivar['(T*)-F12AB CORRECTION ENERGY'] = Decimal(mobj.group(3)) * \
            psivar['MP2-F12 CORRELATION ENERGY'] / psivar['MP2 CORRELATION ENERGY']
        psivar['CCSD(T*)-F12B CORRELATION ENERGY'] = Decimal(mobj.group(2)) + psivar['(T*)-F12AB CORRECTION ENERGY']
        psivar['CCSD(T*)-F12B TOTAL ENERGY'] = psivar['CCSD(T*)-F12B CORRELATION ENERGY'] + psivar['HF-CABS TOTAL ENERGY']

    mobj = re.search(
        r'^\s*' + r'CCSD-F12c triplet pair energy\s+' + NUMBER + '\s*' +
        r'^\s*' + r'CCSD-F12c correlation energy\s+' + NUMBER + '\s*' +
        r'^\s*' + r'Triples \(T\) contribution\s+' + NUMBER + '\s*$',
        outtext, re.MULTILINE)
    if mobj:
        #print('matched ccsd(t)-f12c')
        psivar['CCSD-F12C CORRELATION ENERGY'] = mobj.group(2)
        psivar['CCSD-F12C SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(1)) * \
            Decimal(2) / Decimal(3)
        psivar['CCSD-F12C OPPOSITE-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(2)) - \
            psivar['CCSD-F12C SAME-SPIN CORRELATION ENERGY']
        psivar['CCSD-F12C TOTAL ENERGY'] = Decimal(mobj.group(2)) + psivar['HF-CABS TOTAL ENERGY']
        psivar['(T)-F12C CORRECTION ENERGY'] = mobj.group(3)
        psivar['CCSD(T)-F12C CORRELATION ENERGY'] = Decimal(mobj.group(2)) + Decimal(mobj.group(3))
        psivar['CCSD(T)-F12C TOTAL ENERGY'] = psivar['CCSD(T)-F12C CORRELATION ENERGY'] + psivar['HF-CABS TOTAL ENERGY']
        psivar['(T*)-F12C CORRECTION ENERGY'] = Decimal(mobj.group(3)) * \
            psivar['MP2-F12 CORRELATION ENERGY'] / psivar['MP2 CORRELATION ENERGY']
        psivar['CCSD(T*)-F12C CORRELATION ENERGY'] = Decimal(mobj.group(2)) + psivar['(T*)-F12C CORRECTION ENERGY']
        psivar['CCSD(T*)-F12C TOTAL ENERGY'] = psivar['CCSD(T*)-F12C CORRELATION ENERGY'] + psivar['HF-CABS TOTAL ENERGY']

    # Process Completion
    mobj = re.search(
        r'^\s*' + r'Variable memory released' + r'\s+$',
        outtext, re.MULTILINE)
    if mobj:
        psivar['SUCCESS'] = True

    # Process CURRENT energies (TODO: needs better way)
    if 'HF TOTAL ENERGY' in psivar:
        psivar['CURRENT REFERENCE ENERGY'] = psivar['HF TOTAL ENERGY']
        psivar['CURRENT ENERGY'] = psivar['HF TOTAL ENERGY']

    if 'HF-CABS TOTAL ENERGY' in psivar:
        psivar['CURRENT REFERENCE ENERGY'] = psivar['HF-CABS TOTAL ENERGY']
        psivar['CURRENT ENERGY'] = psivar['HF-CABS TOTAL ENERGY']

    if 'MP2 TOTAL ENERGY' in psivar and 'MP2 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP2 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP2 TOTAL ENERGY']

    if 'MP2-F12 TOTAL ENERGY' in psivar and 'MP2-F12 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP2-F12 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP2-F12 TOTAL ENERGY']

    if 'CCSD TOTAL ENERGY' in psivar and 'CCSD CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD TOTAL ENERGY']

    if 'CCSD-F12A TOTAL ENERGY' in psivar and 'CCSD-F12A CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD-F12A CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD-F12A TOTAL ENERGY']

    if 'CCSD-F12B TOTAL ENERGY' in psivar and 'CCSD-F12B CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD-F12B CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD-F12B TOTAL ENERGY']

    if 'CCSD-F12C TOTAL ENERGY' in psivar and 'CCSD-F12C CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD-F12C CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD-F12C TOTAL ENERGY']

    if 'CCSD(T) TOTAL ENERGY' in psivar and 'CCSD(T) CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD(T) CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD(T) TOTAL ENERGY']

    if 'CCSD(T)-F12A TOTAL ENERGY' in psivar and 'CCSD(T)-F12A CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD(T)-F12A CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD(T)-F12A TOTAL ENERGY']

    if 'CCSD(T)-F12B TOTAL ENERGY' in psivar and 'CCSD(T)-F12B CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD(T)-F12B CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD(T)-F12B TOTAL ENERGY']

    if 'CCSD(T)-F12C TOTAL ENERGY' in psivar and 'CCSD(T)-F12C CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD(T)-F12C CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD(T)-F12C TOTAL ENERGY']

    return psivar, psivar_coord, psivar_grad


class Infile(qcformat.InputFormat2):

    def __init__(self, mem, mol, mtd, der, opt):
        qcformat.InputFormat2.__init__(self, mem, mol, mtd, der, opt)

        #print self.method, self.molecule.nactive_fragments()
        if ('sapt' in self.method or 'mp2c' in self.method) and self.molecule.nactive_fragments() != 2:
            raise FragmentCountError("""Requested molecule has %d, not 2, fragments.""" % (self.molecule.nactive_fragments()))

#        # memory in MB --> MW
#        self.memory = int(math.ceil(mem / 8.0))
        # auxiliary basis sets
        [self.unaugbasis, self.augbasis, self.auxbasis] = self.corresponding_aux_basis()

    def muster_basis_options(self):
        text = ''
        lowername = self.method.lower()
        options = defaultdict(lambda: defaultdict(dict))
    
        options['BASIS']['ORBITAL']['value'] = self.basis
    
        # this f12 basis setting may be totally messed up
        if self.method in ['ccsd(t)-f12-optri']:
            if self.basis == 'cc-pvdz-f12':
                options['BASIS']['JKFIT']['value'] = 'aug-cc-pvtz/jkfit'
                options['BASIS']['JKFITC']['value'] = self.basis + '/optri'
                options['BASIS']['MP2FIT']['value'] = 'aug-cc-pvtz/mp2fit'
        elif self.method in ['ccsd(t)-f12-cabsfit']:
            if self.unaugbasis and self.auxbasis:
                #options['BASIS']['JKFIT']['value'] = self.auxbasis + '/jkfit'
                #options['BASIS']['JKFITB']['value'] = self.unaugbasis + '/jkfit'
                #options['BASIS']['MP2FIT']['value'] = self.auxbasis + '/mp2fit'
                #options['BASIS']['DFLHF']['value'] = self.auxbasis + '/jkfit'
                options['BASIS']['JKFITC']['value'] = 'aug-cc-pv5z/mp2fit'
            else:
                raise ValidationError("""Auxiliary basis not predictable from orbital basis '%s'""" % (self.basis))
        elif ('df-' in self.method) or ('f12' in self.method) or (self.method in ['mp2c', 'dft-sapt', 'dft-sapt-pbe0acalda']):
            if self.unaugbasis and self.auxbasis:
                options['BASIS']['JKFIT']['value'] = self.auxbasis + '/jkfit'
                options['BASIS']['JKFITB']['value'] = self.unaugbasis + '/jkfit'
                options['BASIS']['MP2FIT']['value'] = self.auxbasis + '/mp2fit'
                options['BASIS']['DFLHF']['value'] = self.auxbasis + '/jkfit'
            else:
                raise ValidationError("""Auxiliary basis not predictable from orbital basis '%s'""" % (self.basis))
        return text, options

    def prepare_basis_for_molpro(self):
        text = ''
    
        for opt, val in self.options['BASIS'].items():
                #print opt, val['value']
                #print molpro_basissets.altbasis.keys()
                if not text:
                    text += """basis={\n"""
                try:
                    # jaxz, maxz, etc.
                    for line in molpro_basissets.altbasis[val['value']]:
                        text += """%s\n""" % (line)
                    text += '\n'
                except KeyError:
                    # haxz
                    if val['value'].startswith('heavy-aug-'):
                        text += """set,%s; default,%s,H=%s\n""" % (opt.lower(), self.augbasis, self.unaugbasis)
                    # xz, axz, 6-31g*
                    else:
                        text += """set,%s; default,%s\n""" % (opt.lower(), val['value'])
    
        if text:
            text += """}\n\n"""

        return text

    def format_infile_string(self):
        """

        """
        # Handle memory and comment
        memcmd, _memkw = """***, %s\nmemory,%d,m\n""" % (self.molecule.tagline, int(math.ceil(self.memory / 8.0))), {}

        # Handle molecule and basis set
        molcmd, _molkw = self.molecule.format_molecule_for_molpro(), {}



        # format global convergence directions
#        text += self.format_global_parameters()
        _cdscmd, cdskw = muster_cdsgroup_options(self.method)

        # Handle calc type and quantum chemical method
        mdccmd, mdckw, mdcls = procedures['energy'][self.method](self.method, self.dertype, self.molecule)
        _bascmd, baskw = self.muster_basis_options()


#        # format options
#        optcmd = qcdb.options.prepare_options_for_psi4(mdckw)

# make options from imdb only user options (currently non-existent). set basis and castup from here.
        # Handle driver vs input/default keyword reconciliation
        userkw = self.options
#        userkw = p4util.prepare_options_for_modules()
        #userkw = qcdb.options.reconcile_options(userkw, memkw)
        #userkw = qcdb.options.reconcile_options(userkw, molkw)
        userkw = options.reconcile_options2(userkw, cdskw)
        userkw = options.reconcile_options2(userkw, baskw)
        #userkw = qcdb.options.reconcile_options(userkw, psikw)
        userkw = options.reconcile_options2(userkw, mdckw)

        # Handle conversion of psi4 keyword structure into cfour format
        #optcmdB = options.prepare_options_for_psi4(userkw)
        optcmd = prepare_options_for_molpro(userkw, mdcls)
        bascmd, _baskw = self.prepare_basis_for_molpro(), {} #self.options['BASIS']), {}

        # Handle text to be passed untouched
        litcmd = """\nshow[1,20f20.12],ee*,ce*,te*\nshow[1,60f20.12],_E*\n\n"""


        # Assemble infile pieces
        return memcmd + molcmd + bascmd + optcmd + mdccmd + litcmd


def muster_cdsgroup_options(name):
    text = ''
    lowername = name.lower()
    options = defaultdict(lambda: defaultdict(dict))

    options['GTHRESH']['ZERO']['value'] = 1.0e-14
    options['GTHRESH']['ONEINT']['value'] = 1.0e-14
    options['GTHRESH']['TWOINT']['value'] = 1.0e-14
    options['GTHRESH']['ENERGY']['value'] = 1.0e-9

    if name in ['mp2c', 'dft-sapt-shift', 'dft-sapt', 'dft-sapt-pbe0ac', 'dft-sapt-pbe0acalda']:
        options['GTHRESH']['ENERGY']['value'] = 1.0e-8
        options['GTHRESH']['ORBITAL']['value'] = 1.0e-8
        options['GTHRESH']['GRID']['value'] = 1.0e-8
    elif name in ['b3lyp', 'b3lyp-d', 'df-b3lyp', 'df-b3lyp-d']:
        options['GTHRESH']['ENERGY']['value'] = 1.0e-8
        options['GTHRESH']['ORBITAL']['value'] = 1.0e-7
        options['GTHRESH']['GRID']['value'] = 1.0e-8
    else:
        pass

    return text, options





def prepare_options_for_molpro(options, proc):
    """Function to take the full snapshot of the liboptions object
    encoded in dictionary *options*, find the options directable toward
    Cfour (options['CFOUR']['CFOUR_**']) that aren't default, then write
    a CFOUR deck with those options.
    Note that unlike the cfour version, this uses complete options deck.

    """
    text = ''

    if len(options['GTHRESH']) > 0:
        text += 'gthresh'
        for opt, val in options['GTHRESH'].items():
            text += """,%s=%s""" % (opt, val['value'])
        text += '\n\n'

    for item in proc:
        if len(options[item.upper()]) > 0:
            text += """{%s%s}\n""" % (item, options[item.upper()]['OPTIONS']['value'])
        else:
            text += """%s\n""" % (item)

    if text:
        text += '\n'

    return text


def muster_modelchem(name, dertype, mol):
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
    proc = []

    if dertype == 0:
        pass
    else:
        raise ValidationError("""Requested Psi4 dertype %d is not available.""" % (dertype))

    if lowername == 'mp2':
        pass
        options['GLOBALS']['FREEZE_CORE']['value'] = True
        options['SCF']['SCF_TYPE']['value'] = 'direct'
        options['MP2']['MP2_TYPE']['value'] = 'conv'
        text += """mp2')\n\n"""

    elif lowername == 'ccsd(t)-f12':
        proc.append('rhf')
        proc.append('ccsd(t)-f12')
        options['CCSD(T)-F12']['OPTIONS']['value'] = ',df_basis=mp2fit,df_basis_exch=jkfitb,ri_basis=jkfitb'

    elif lowername == 'ccsd(t)-f12c':
        proc.append('rhf')
        proc.append('ccsd(t)-f12c')
        options['CCSD(T)-F12C']['OPTIONS']['value'] = ',df_basis=mp2fit,df_basis_exch=jkfitb,ri_basis=jkfitb'

    elif lowername == 'ccsd(t)-f12-optri':
        proc.append('rhf')
        proc.append('ccsd(t)-f12')
        options['CCSD(T)-F12']['OPTIONS']['value'] = ',df_basis=mp2fit,df_basis_exch=jkfit,ri_basis=jkfitc'

    elif lowername == 'ccsd(t)-f12-cabsfit':
        proc.append('rhf')
        proc.append('ccsd(t)-f12')
        options['CCSD(T)-F12']['OPTIONS']['value'] = ',df_basis=jkfitc,df_basis_exch=jkfitc,ri_basis=jkfitc'

    elif lowername == 'mp2c':
        proc.append('gdirect')
        proc.append(mol.extract_fragments(1, 2).format_molecule_for_molpro())
        proc.append('df-hf,')
        proc.append('df-ks,')
        proc.append('sapt; monomerA')
        options['DF-HF,']['OPTIONS']['value'] = """basis=jkfit,locorb=0; start,atdens; save,1101.2"""
        options['DF-KS,']['OPTIONS']['value'] = """lhf,df_basis=dflhf,basis_coul=jkfitb,basis_exch=jkfitb; dftfac,1.0; start,1101.2; save,2101.2"""

        proc.append(mol.extract_fragments(2, 1).format_molecule_for_molpro())
        proc.append('df-hf')
        proc.append('df-ks')
        proc.append('sapt; monomerB')
        options['DF-HF']['OPTIONS']['value'] = """,basis=jkfit,locorb=0; start,atdens; save,1102.2"""
        options['DF-KS']['OPTIONS']['value'] = """,lhf,df_basis=dflhf,basis_coul=jkfitb,basis_exch=jkfitb; dftfac,1.0; start,1102.2; save,2102.2"""

        proc.append(mol.format_molecule_for_molpro())
        proc.append('sapt; intermol')
        options['SAPT; INTERMOL']['OPTIONS']['value'] = """,saptlevel=3,ca=2101.2,cb=2102.2,icpks=0,fitlevel=3,nlexfac=0.0,cfac=0.0; dfit,basis_coul=jkfit,basis_exch=jkfit,cfit_scf=3"""

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

    return text, options, proc

procedures = {
    'energy': {
        'mp2c'           : muster_modelchem,
        'ccsd(t)-f12'    : muster_modelchem,
        'ccsd(t)-f12c'   : muster_modelchem,
        'ccsd(t)-f12-optri' : muster_modelchem,
        'ccsd(t)-f12-cabsfit' : muster_modelchem,
        #'sapt0'         : muster_modelchem,
        #'sapt2+'        : muster_modelchem,
        #'sapt2+(3)'     : muster_modelchem,
        #'sapt2+3(ccd)'  : muster_modelchem,
    }
}

qcmtdIN = procedures['energy']


def psi4_list():
    """Return an array of Psi4 methods with energies.

    """
    return procedures['energy'].keys()
