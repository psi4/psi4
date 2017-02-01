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
import struct
from collections import defaultdict
from decimal import Decimal
from .pdict import PreservingDict
from .periodictable import *
from .physconst import *
from .exceptions import *
from .molecule import Molecule
from .orient import OrientMols
from .options import conv_float2negexp


def harvest_output(outtext):
    """Function to separate portions of a CFOUR output file *outtest*,
    divided by xjoda.

    """
    pass_psivar = []
    pass_coord = []
    pass_grad = []

    for outpass in re.split(r'--invoking executable xjoda', outtext, re.MULTILINE):
        psivar, c4coord, c4grad = harvest_outfile_pass(outpass)
        pass_psivar.append(psivar)
        pass_coord.append(c4coord)
        pass_grad.append(c4grad)

        #print '\n\nXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n'
        #print outpass
        #print psivar, c4coord, c4grad
        #print psivar, c4grad
        #print '\n\nxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n'

    retindx = -1 if pass_coord[-1] else -2

#    print '    <<<  C4 PSIVAR  >>>'
#    for item in pass_psivar[retindx]:
#        print('       %30s %16.8f' % (item, pass_psivar[retindx][item]))
#    print '    <<<  C4 COORD   >>>'
#    for item in pass_coord[retindx]:
#        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    print '    <<<   C4 GRAD   >>>'
#    for item in pass_grad[retindx]:
#        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

    return pass_psivar[retindx], pass_coord[retindx], pass_grad[retindx]


def harvest_outfile_pass(outtext):
    """Function to read CFOUR output file *outtext* and parse important
    quantum chemical information from it in

    """
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None

#    TODO: BCC
#          CI
#          QCISD(T)
#          other ROHF tests
#          vcc/ecc

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    # Process NRE
    mobj = re.search(r'^\s+' + r'(?:Nuclear repulsion energy :)' + r'\s+' + NUMBER + r'\s+a\.u\.\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched nre')
        psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)

    # Process SCF
    mobj = re.search(
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched scf1')
        psivar['SCF TOTAL ENERGY'] = mobj.group(1)

    mobj = re.search(
        r'^\s+' + r'(?:E\(SCF\)=)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched scf2')
        psivar['SCF TOTAL ENERGY'] = mobj.group(1)

    mobj = re.search(
        r'^\s+' + r'(?:SCF has converged.)' + r'\s*$' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched scf3')
        psivar['SCF TOTAL ENERGY'] = mobj.group(1)

    # Process MP2
    mobj = re.search(
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched mp2r')
        psivar['MP2 SAME-SPIN CORRELATION ENERGY'] = 2 * Decimal(mobj.group(1))
        psivar['MP2 OPPOSITE-SPIN CORRELATION ENERGY'] = mobj.group(2)
        psivar['MP2 CORRELATION ENERGY'] = 2 * Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        psivar['MP2 TOTAL ENERGY'] = mobj.group(4)

    mobj = re.search(
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(BB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched mp2u')
        psivar['MP2 SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        psivar['MP2 OPPOSITE-SPIN CORRELATION ENERGY'] = mobj.group(3)
        psivar['MP2 CORRELATION ENERGY'] = Decimal(mobj.group(1)) + \
            Decimal(mobj.group(2)) + Decimal(mobj.group(3))
        psivar['MP2 TOTAL ENERGY'] = mobj.group(5)

    mobj = re.search(
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(BB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(SINGLE\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched mp2ro')
        psivar['MP2 SAME-SPIN CORRELATION ENERGY'] = Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        psivar['MP2 OPPOSITE-SPIN CORRELATION ENERGY'] = mobj.group(3)
        psivar['MP2 SINGLES ENERGY'] = mobj.group(4)
        psivar['MP2 CORRELATION ENERGY'] = Decimal(mobj.group(1)) + \
            Decimal(mobj.group(2)) + Decimal(mobj.group(3)) + Decimal(mobj.group(4))
        psivar['MP2 TOTAL ENERGY'] = mobj.group(6)

    # Process MP3
    mobj = re.search(
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp3r')
        dmp2 = Decimal(mobj.group(1))
        dmp3 = Decimal(mobj.group(3))
        psivar['MP2 CORRELATION ENERGY'] = dmp2
        psivar['MP2 TOTAL ENERGY'] = mobj.group(2)
        psivar['MP3 CORRELATION ENERGY'] = dmp2 + dmp3
        psivar['MP3 TOTAL ENERGY'] = mobj.group(4)
        psivar['MP2.5 CORRELATION ENERGY'] = dmp2 + Decimal('0.500000000000') * dmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']

    mobj = re.search(
        r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp3ro')
        dmp2 = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp3 = Decimal(mobj.group(5)) + Decimal(mobj.group(7))
        psivar['MP3 CORRELATION ENERGY'] = dmp2 + dmp3
        psivar['MP3 TOTAL ENERGY'] = mobj.group(8)
        psivar['MP2.5 CORRELATION ENERGY'] = dmp2 + Decimal('0.500000000000') * dmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']

    # Process MP4
    mobj = re.search(
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:Q-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp4r')
        dmp2 = Decimal(mobj.group(1))
        dmp3 = Decimal(mobj.group(3))
        dmp4sdq = Decimal(mobj.group(5)) + Decimal(mobj.group(7)) + Decimal(mobj.group(9))
        psivar['MP2 CORRELATION ENERGY'] = dmp2
        psivar['MP2 TOTAL ENERGY'] = mobj.group(2)
        psivar['MP3 CORRELATION ENERGY'] = dmp2 + dmp3
        psivar['MP3 TOTAL ENERGY'] = mobj.group(4)
        psivar['MP2.5 CORRELATION ENERGY'] = dmp2 + Decimal('0.500000000000') * dmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']
        psivar['MP4(SDQ) CORRELATION ENERGY'] = dmp2 + dmp3 + dmp4sdq
        psivar['MP4(SDQ) TOTAL ENERGY'] = mobj.group(10)

    mobj = re.search(
        r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:L-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:NL-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp4ro')
        dmp2 = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp3 = Decimal(mobj.group(5)) + Decimal(mobj.group(7))
        dmp4sdq = Decimal(mobj.group(9)) + Decimal(mobj.group(11))
        psivar['MP2 CORRELATION ENERGY'] = dmp2
        psivar['MP2 TOTAL ENERGY'] = mobj.group(4)
        psivar['MP3 CORRELATION ENERGY'] = dmp2 + dmp3
        psivar['MP3 TOTAL ENERGY'] = mobj.group(8)
        psivar['MP2.5 CORRELATION ENERGY'] = dmp2 + Decimal('0.500000000000') * dmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']
        psivar['MP4(SDQ) CORRELATION ENERGY'] = dmp2 + dmp3 + dmp4sdq
        psivar['MP4(SDQ) TOTAL ENERGY'] = mobj.group(12)

    mobj = re.search(
        r'^\s+' + r'(?:D-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:Q-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:T-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp4tr')
        dmp4sdq = Decimal(mobj.group(1)) + Decimal(mobj.group(3)) + Decimal(mobj.group(5))
        dmp4t = Decimal(mobj.group(7))
        psivar['MP4(SDQ) CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY'] + dmp4sdq
        psivar['MP4(SDQ) TOTAL ENERGY'] = mobj.group(6)
        psivar['MP4(T) CORRECTION ENERGY'] = dmp4t
        psivar['MP4(SDTQ) CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY'] + dmp4sdq + dmp4t
        psivar['MP4(SDTQ) TOTAL ENERGY'] = mobj.group(8)
        psivar['MP4 CORRELATION ENERGY'] = psivar['MP4(SDTQ) CORRELATION ENERGY']
        psivar['MP4 TOTAL ENERGY'] = psivar['MP4(SDTQ) TOTAL ENERGY']

    mobj = re.search(
        r'^\s+' + r'(?:L-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:NL-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:WT12-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:T-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp4tro')
        dmp4sdq = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp4t = Decimal(mobj.group(5)) + Decimal(mobj.group(7))  # TODO: WT12 with T, not SDQ?
        psivar['MP4(SDQ) CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY'] + dmp4sdq
        psivar['MP4(SDQ) TOTAL ENERGY'] = mobj.group(4)
        psivar['MP4(T) CORRECTION ENERGY'] = dmp4t
        psivar['MP4(SDTQ) CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY'] + dmp4sdq + dmp4t
        psivar['MP4(SDTQ) TOTAL ENERGY'] = mobj.group(8)
        psivar['MP4 CORRELATION ENERGY'] = psivar['MP4(SDTQ) CORRELATION ENERGY']
        psivar['MP4 TOTAL ENERGY'] = psivar['MP4(SDTQ) TOTAL ENERGY']

    # Process CC Iterations
    mobj = re.search(
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+DIIS\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s*(?:A miracle (?:has come|come) to pass. The CC iterations have converged.)\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched cc with full %s iterating %s' % (mobj.group('fullCC'), mobj.group('iterCC')))
        psivar['%s CORRELATION ENERGY' % (mobj.group('iterCC'))] = mobj.group(3)
        psivar['%s TOTAL ENERGY' % (mobj.group('iterCC'))] = mobj.group(4)

    # Process CC(T)
    mobj = re.search(
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E\(CCSD\))' + r'\s+=\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E\(CCSD\(T\)\))' + r'\s+=\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched ccsd(t) vcc')
        psivar['SCF TOTAL ENERGY'] = mobj.group(1)
        psivar['CCSD TOTAL ENERGY'] = mobj.group(2)
        psivar['(T) CORRECTION ENERGY'] = Decimal(mobj.group(3)) - Decimal(mobj.group(2))
        psivar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(3)) - Decimal(mobj.group(1))
        psivar['CCSD(T) TOTAL ENERGY'] = mobj.group(3)

    mobj = re.search(
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s*' + NUMBER + r'\s+a\.u\.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:CCSD energy)' + r'\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:Total perturbative triples energy:)' + r'\s+' + NUMBER + r'\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s+' + r'(?:CCSD\(T\) energy)' + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched ccsd(t) ecc')
        psivar['SCF TOTAL ENERGY'] = mobj.group(1)
        psivar['CCSD TOTAL ENERGY'] = mobj.group(2)
        psivar['(T) CORRECTION ENERGY'] = mobj.group(3)
        psivar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(4)) - Decimal(mobj.group(1))
        psivar['CCSD(T) TOTAL ENERGY'] = mobj.group(4)

    mobj = re.search(
        r'^\s+' + r'(?:CCSD energy)' + r'\s+' + NUMBER + r'\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s+' + r'(?:CCSD\(T\) energy)' + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched ccsd(t) lamb')
        psivar['CCSD TOTAL ENERGY'] = mobj.group(1)
        psivar['(T) CORRECTION ENERGY'] = Decimal(mobj.group(2)) - Decimal(mobj.group(1))
        psivar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(2)) - psivar['SCF TOTAL ENERGY']
        psivar['CCSD(T) TOTAL ENERGY'] = mobj.group(2)

    # Process SCS-CC
    mobj = re.search(
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s*' + r'(?:@CCENRG-I, Correlation energies.)' + r'\s+(?:ECCAA)\s+' + NUMBER + r'\s*' +
        r'^\s+(?:ECCBB)\s+' + NUMBER + '\s*' +
        r'^\s+(?:ECCAB)\s+' + NUMBER + '\s*' +
        r'^\s+(?:Total)\s+' + NUMBER + '\s*',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:  # PRINT=2 to get SCS-CC components
        print('matched scscc')
        psivar['%s SAME-SPIN CORRELATION ENERGY' % (mobj.group('iterCC'))] = Decimal(mobj.group(3)) + Decimal(mobj.group(4))
        psivar['%s OPPOSITE-SPIN CORRELATION ENERGY' % (mobj.group('iterCC'))] = mobj.group(5)
        psivar['%s CORRELATION ENERGY' % (mobj.group('iterCC'))] = mobj.group(6)

    mobj = re.search(
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'Amplitude equations converged in' + r'\s*\d+\s*' + r'iterations.\s*' +
        r'^\s+' + r'The AA contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The BB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The AB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The total correlation energy is\s+' + NUMBER + r'\s+a.u.\s*' +
        r'(?:.*?)' +
        #r'^\s+' + r'The CC iterations have converged.' + r'\s*$',
        r'^\s+' + r'(?:A miracle come to pass. )?' + r'The CC iterations have converged.' + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:  # PRINT=2 to get SCS components
        print('matched scscc2')
        psivar['%s SAME-SPIN CORRELATION ENERGY' % (mobj.group('iterCC'))] = Decimal(mobj.group(3)) + Decimal(mobj.group(4))
        psivar['%s OPPOSITE-SPIN CORRELATION ENERGY' % (mobj.group('iterCC'))] = mobj.group(5)
        psivar['%s CORRELATION ENERGY' % (mobj.group('iterCC'))] = mobj.group(6)

    # Process gradient
    mobj = re.search(
        r'\s+' + r'Molecular gradient' + r'\s*' +
        r'\s+' + r'------------------' + r'\s*' +
        r'\s+' + r'\n' +
        r'(?:(?:\s+[A-Z]+\s*#\d+\s+[xyz]\s+[-+]?\d+\.\d+\s*\n)+)' +  # optional, it seems
        r'\n\n' +  # optional, it seems
        r'((?:\s+[A-Z]+\s*#\d+\s+\d?\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)' +
        r'\n\n' +
        r'\s+' + 'Molecular gradient norm',
        outtext, re.MULTILINE)
    if mobj:
        print('matched molgrad')
        atoms = []
        psivar_grad = []
        for line in mobj.group(1).splitlines():
            lline = line.split()
            atoms.append(lline[0])
            #psivar_gradient.append([Decimal(lline[-3]), Decimal(lline[-2]), Decimal(lline[-1])])
            psivar_grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])

    # Process geometry
    mobj = re.search(
#        r'\s+(?:-+)\s*' +
#        r'^\s+' + r'Z-matrix   Atomic            Coordinates (in bohr)' + r'\s*' +
        r'^\s+' + r'Symbol    Number           X              Y              Z' + r'\s*' +
        r'^\s+(?:-+)\s*' +
        r'((?:\s+[A-Z]+\s+[0-9]+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)' +
        r'^\s+(?:-+)\s*',
        outtext, re.MULTILINE)
    if mobj:
        print('matched geom')
        molxyz = '%d bohr\n\n' % len(mobj.group(1).splitlines())
        for line in mobj.group(1).splitlines():
            lline = line.split()
            molxyz += '%s %16s %16s %16s\n' % (lline[0], lline[-3], lline[-2], lline[-1])
        # Rather a dinky Molecule as no ghost, charge, or multiplicity
        psivar_coord = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)

    # Process atom geometry
    mobj = re.search(
        r'^\s+' + r'@GETXYZ-I,     1 atoms read from ZMAT.' + r'\s*' +
        r'^\s+' + r'[0-9]+\s+([A-Z]+)\s+[0-9]+\s+' + NUMBER + r'\s*',
        outtext, re.MULTILINE)
    if mobj:
        print('matched atom')
        # Dinky Molecule
        molxyz = '1 bohr\n\n%s 0.0 0.0 0.0\n' % (mobj.group(1))
        psivar_coord = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)

    # Process error codes
    mobj = re.search(
        r'^\s*' + r'--executable ' + r'(\w+)' + r' finished with status' + r'\s+' + r'([1-9][0-9]*)',
        outtext, re.MULTILINE)
    if mobj:
        print('matched error')
        psivar['CFOUR ERROR CODE'] = mobj.group(2)

    # Process CURRENT energies (TODO: needs better way)
    if 'SCF TOTAL ENERGY' in psivar:
        psivar['CURRENT REFERENCE ENERGY'] = psivar['SCF TOTAL ENERGY']
        psivar['CURRENT ENERGY'] = psivar['SCF TOTAL ENERGY']

    if 'MP2 TOTAL ENERGY' in psivar and 'MP2 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP2 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP2 TOTAL ENERGY']

    if 'MP3 TOTAL ENERGY' in psivar and 'MP3 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP3 TOTAL ENERGY']

    if 'MP4 TOTAL ENERGY' in psivar and 'MP4 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP4 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP4 TOTAL ENERGY']

#    if ('%s TOTAL ENERGY' % (mobj.group('fullCC')) in psivar) and \
#       ('%s CORRELATION ENERGY' % (mobj.group('fullCC')) in psivar):
#        psivar['CURRENT CORRELATION ENERGY'] = psivar['%s CORRELATION ENERGY' % (mobj.group('fullCC')]
#        psivar['CURRENT ENERGY'] = psivar['%s TOTAL ENERGY' % (mobj.group('fullCC')]

    if 'CC2 TOTAL ENERGY' in psivar and 'CC2 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CC2 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CC2 TOTAL ENERGY']

    if 'CCSD TOTAL ENERGY' in psivar and 'CCSD CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD TOTAL ENERGY']

    if 'CCSD(T) TOTAL ENERGY' in psivar and 'CCSD(T) CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD(T) CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD(T) TOTAL ENERGY']

    if 'CC3 TOTAL ENERGY' in psivar and 'CC3 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CC3 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CC3 TOTAL ENERGY']

    if 'CCSDT TOTAL ENERGY' in psivar and 'CCSDT CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSDT CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSDT TOTAL ENERGY']

    return psivar, psivar_coord, psivar_grad


def harvest(p4Mol, c4out, **largs):
    """Parses all the pieces of output from Cfour: the stdout in
    *c4out* and the contents of various scratch files like GRD stored
    in their namesake keys in *largs*. Since all Cfour output uses
    its own orientation and atom ordering for the given molecule,
    a qcdb.Molecule *p4Mol*, if supplied, is used to transform the
    Cfour output back into consistency with *p4Mol*.

    """
    # Collect results from output file and subsidiary files
    outPsivar, outMol, outGrad = harvest_output(c4out)

    if 'GRD' in largs:
        grdMol, grdGrad = harvest_GRD(largs['GRD'])
    else:
        grdMol, grdGrad = None, None

    if 'FCMFINAL' in largs:
        fcmHess = harvest_FCM(largs['FCMFINAL'])
    else:
        fcmHess = None

    if 'DIPOL' in largs:
        dipolDip = harvest_DIPOL(largs['DIPOL'])
    else:
        dipolDip = None

    # Reconcile the coordinate information: several cases
    #   Case                            p4Mol   GRD      Check consistency           Apply orientation?     ReturnMol (1-19-2014)
    #   sp with mol thru cfour {}       None    None              outMol             N.C.                   outMol
    #   opt with mol thru cfour {}      None    grdMol            outMol && grdMol   N.C.                   grdMol
    #   sp with mol thru molecule {}    p4Mol   None     p4Mol && outMol             p4Mol <-- outMol       p4Mol (same as input arg)
    #   opt with mol thru molecule {}   p4Mol   grdMol   p4Mol && outMol && grdMol   p4Mol <-- grdMol       p4Mol (same as input arg)

    if outMol:
        if grdMol:
            if abs(outMol.nuclear_repulsion_energy() - grdMol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValidationError("""Cfour outfile (NRE: %f) inconsistent with Cfour GRD (NRE: %f).""" % \
                        (outMol.nuclear_repulsion_energy(), grdMol.nuclear_repulsion_energy()))
        if p4Mol:
            if abs(outMol.nuclear_repulsion_energy() - p4Mol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValidationError("""Cfour outfile (NRE: %f) inconsistent with Psi4 input (NRE: %f).""" % \
                    (outMol.nuclear_repulsion_energy(), p4Mol.nuclear_repulsion_energy()))
    else:
        raise ValidationError("""No coordinate information extracted from Cfour output.""")

#    print '    <<<   [1] P4-MOL   >>>'
#    if p4Mol:
#        p4Mol.print_out_in_bohr()
#    print '    <<<   [2] C4-OUT-MOL   >>>'
#    if outMol:
#        outMol.print_out_in_bohr()
#    print '    <<<   [3] C4-GRD-MOL   >>>'
#    if grdMol:
#        grdMol.print_out_in_bohr()

    # Set up array reorientation object
    if p4Mol and grdMol:
        p4c4 = OrientMols(p4Mol, grdMol)
        oriCoord = p4c4.transform_coordinates2(grdMol)
        oriGrad = p4c4.transform_gradient(grdGrad)
        oriDip = None if dipolDip is None else p4c4.transform_vector(dipolDip)
    elif p4Mol and outMol:
        p4c4 = OrientMols(p4Mol, outMol)
        oriCoord = p4c4.transform_coordinates2(outMol)
        oriGrad = None
        oriDip = None if dipolDip is None else p4c4.transform_vector(dipolDip)
    elif outMol:
        oriCoord = None
        oriGrad = None
        oriDip = None if dipolDip is None else dipolDip

#    print p4c4
#    print '    <<<   [4] C4-ORI-MOL   >>>'
#    if oriCoord is not None:
#        for item in oriCoord:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#
#    print '    <<<   [1] C4-GRD-GRAD   >>>'
#    if grdGrad is not None:
#        for item in grdGrad:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    print '    <<<   [2] C4-ORI-GRAD   >>>'
#    if oriGrad is not None:
#        for item in oriGrad:
#            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

    retMol = None if p4Mol else grdMol

    if oriDip:
        outPsivar['CURRENT DIPOLE X'] = str(oriDip[0] * psi_dipmom_au2debye)
        outPsivar['CURRENT DIPOLE Y'] = str(oriDip[1] * psi_dipmom_au2debye)
        outPsivar['CURRENT DIPOLE Z'] = str(oriDip[2] * psi_dipmom_au2debye)

    if oriGrad:
        retGrad = oriGrad
    elif grdGrad:
        retGrad = grdGrad
    else:
        retGrad = None

    return outPsivar, retGrad, retMol


def harvest_GRD(grd):
    """Parses the contents *grd* of the Cfour GRD file into the gradient
    array and coordinate information. The coordinate info is converted
    into a rather dinky Molecule (no charge, multiplicity, or fragment),
    but this is these coordinates that govern the reading of molecule
    orientation by Cfour. Return qcdb.Molecule and gradient array.

    """
    grd = grd.splitlines()
    Nat = int(grd[0].split()[0])
    molxyz = '%d bohr\n\n' % (Nat)

    grad = []
    for at in range(Nat):
        mline = grd[at + 1].split()
        el = 'GH' if int(float(mline[0])) == 0 else z2el[int(float(mline[0]))]
        molxyz += '%s %16s %16s %16s\n' % (el, mline[-3], mline[-2], mline[-1])
        lline = grd[at + 1 + Nat].split()
        grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])
    mol = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)

    return mol, grad


def harvest_zmat(zmat):
    """Parses the contents of the Cfour ZMAT file into array and
    coordinate information. The coordinate info is converted into a
    rather dinky Molecule (no fragment, but does read charge, mult,
    unit). Return qcdb.Molecule. Written for findif zmat* where
    geometry always Cartesian and Bohr.

    """
    zmat = zmat.splitlines()[1:]  # skip comment line
    Nat = 0
    readCoord = True
    isBohr = ''
    charge = 0
    mult = 1
    molxyz = ''
    cgeom = []
    for line in zmat:
        if line.strip() == '':
            readCoord = False
        elif readCoord:
            lline = line.split()
            molxyz += line + '\n'
            Nat += 1
        else:
            if line.find('CHARGE') > -1:
                idx = line.find('CHARGE')
                charge = line[idx + 7:]
                idxc = charge.find(',')
                if idxc > -1:
                    charge = charge[:idxc]
                charge = int(charge)
            if line.find('MULTIPLICITY') > -1:
                idx = line.find('MULTIPLICITY')
                mult = line[idx + 13:]
                idxc = mult.find(',')
                if idxc > -1:
                    mult = mult[:idxc]
                mult = int(mult)
            if line.find('UNITS=BOHR') > -1:
                isBohr = ' bohr'

    molxyz = '%d%s\n%d %d\n' % (Nat, isBohr, charge, mult) + molxyz
    mol = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)

    return mol


def harvest_FCM(fcm):
    """Parses the contents *fcm* of the Cfour FCMFINAL file into a hessian array.

    """
    fcm = fcm.splitlines()
    Nat = int(fcm[0].split()[0])
    Ndof = int(fcm[0].split()[1])

    empty = True
    hess = []
    for df in range(Ndof):
        for at in range(Nat):
            lline = fcm[Ndof * at + at + 1].split()
            if empty:
                if (abs(float(lline[0])) > 1.0e-8) or \
                   (abs(float(lline[1])) > 1.0e-8) or \
                   (abs(float(lline[2])) > 1.0e-8):
                    empty = False
            fcm.append([float(lline[0]), float(lline[1]), float(lline[2])])

    return None if empty else hess


def harvest_DIPOL(dipol):
    """Parses the contents *dipol* of the Cfour DIPOL file into a dipol vector.

    """
    dipol = dipol.splitlines()
    lline = dipol[0].split()
    dip = [float(lline[0]), float(lline[1]), float(lline[2])]

    #return None if empty else dip
    return dip


def muster_memory(mem):
    """Transform input *mem* in MB into psi4-type options for cfour.

    """
    text = ''

    # prepare memory keywords to be set as c-side keywords
    options = defaultdict(lambda: defaultdict(dict))
    options['CFOUR']['CFOUR_MEMORY_SIZE']['value'] = int(mem)
    options['CFOUR']['CFOUR_MEM_UNIT']['value'] = 'MB'

    for item in options['CFOUR']:
        options['CFOUR'][item]['clobber'] = True
    return text, options

#   Ways of modifying a computation
#   global:     set global c-side option
#   local:      set local c-side option
#   kwarg:      set kwarg
#   i-local:    set global=local c-side option to an interface module
#   ro-def:     code uses default entirely specified by read_options
#   module-def: code uses default that is complex mixture of read_options settings
#   i-def:      interfaced code uses defaults not entirely expressed in read_options
#   driver-def: driver code sets complex defaults
#
#   Pure psi4 operation
#   kwarg ~= local > global > driver-def > module-def > ro-def
#
#   Interfaced psi4 operation
#   kwarg ~= i-local > local > global > driver-def > i-def

#   P4 infrastructure replacing interfaced infrastructure (mol, basis, mem) where unavoidable overlap in how things are specified (mult in mol{} vs keyword) is treated as a clobber & complain if conflict VS P4 infrastructure as an aliased/convenient leak into interfaced infrastructure (psi) and is strictly no clobber or complain.


def muster_psi4options(opt):
    """Translate psi4 keywords *opt* that have been explicitly set into
    their Cfour counterparts. Since explicitly set Cfour module keyword
    values will always be used preferentially to these inferred from
    psi4, the 'clobber' property is set to False.

    """
    text = ''
    options = defaultdict(lambda: defaultdict(dict))

    if 'GLOBALS' in opt:
        if 'PUREAM' in opt['GLOBALS']:
            options['CFOUR']['CFOUR_SPHERICAL']['value'] = \
                opt['MINTS']['PUREAM']['value']

    if 'SCF' in opt:
        if 'REFERENCE' in opt['SCF']:
            options['CFOUR']['CFOUR_REFERENCE']['value'] = \
                {'RHF': 'RHF',
                 'UHF': 'UHF',
                 'ROHF': 'ROHF'}[opt['SCF']['REFERENCE']['value']]

        if 'D_CONVERGENCE' in opt['SCF']:
            options['CFOUR']['CFOUR_SCF_CONV']['value'] = \
                conv_float2negexp(opt['SCF']['D_CONVERGENCE']['value'])

        if 'MAXITER' in opt['SCF']:
            options['CFOUR']['CFOUR_SCF_MAXCYC']['value'] = \
                opt['SCF']['MAXITER']['value']

        if 'DAMPING_PERCENTAGE' in opt['SCF']:
            options['CFOUR']['CFOUR_SCF_DAMPING']['value'] = \
                int(10 * opt['SCF']['DAMPING_PERCENTAGE']['value'])

    for item in options['CFOUR']:
        options['CFOUR'][item]['clobber'] = False
    return text, options

# Philosophy break:
#   Specification options
#   Massaging options

#   * No program's defaults should be tampered with w/o provokation

#   want all defaults applied to all programs, so p4 scf_conv is 5 and c4 scf_conv is 5
#   want separate regimes, so conv 6 covers all the p4 parts and cfour_conv = 8 covers the c4 parts
#   want mixture, so basis gets applied to c4 but others don't
#   first case, when options specified explicitly

#   [scf][d_convergence]    [cfour][cfour_scf_conv]     what happens?
#   8 from opt()            7 by default
#   6 from set {...}        7 by default                6 (guideline that psi4 format converts when clear)
#   8 from opt()            5 from set {...}            5 (local trumps)
#   6 from set {...}        5 from set {...}            5 (local trumps)
#
#   energy(name)            [cfour][cfour_calc_level]
#   c4-scf                  SCF by default
#   c4-scf                  CCSD from set {...}


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
        if lowername == 'cfour':
            pass  # permit clean operation of sandwich mode
        else:
            options['CFOUR']['CFOUR_DERIV_LEVEL']['value'] = 'ZERO'
    elif dertype == 1:
        options['CFOUR']['CFOUR_DERIV_LEVEL']['value'] = 'FIRST'
    elif dertype == 2:
        options['CFOUR']['CFOUR_DERIV_LEVEL']['value'] = 'SECOND'
    else:
        raise ValidationError("""Requested Cfour dertype %d is not available.""" % (dertype))

    if lowername == 'cfour':
        pass
    elif lowername == 'c4-scf':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'SCF'

    elif lowername == 'c4-mp2':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'MP2'

    elif lowername == 'c4-mp3':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'MP3'

    elif lowername == 'c4-mp4(sdq)':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'SDQ-MP4'

    elif lowername == 'c4-mp4':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'MP4'

    elif lowername == 'c4-cc2':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CC2'

    elif lowername == 'c4-ccsd':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSD'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'ECC'

    elif lowername == 'c4-cc3':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CC3'

    elif lowername == 'c4-ccsd(t)':
        # Can't use (T) b/c bug in xsymcor lops it off
        #options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSD(T)'
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSD[T]'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'ECC'

    elif lowername == 'c4-ccsdt':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSDT'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'ECC'

    else:
        raise ValidationError("""Requested Cfour computational methods %d is not available.""" % (lowername))

    # Set clobbering
    if 'CFOUR_DERIV_LEVEL' in options['CFOUR']:
        options['CFOUR']['CFOUR_DERIV_LEVEL']['clobber'] = True
        options['CFOUR']['CFOUR_DERIV_LEVEL']['superclobber'] = True
    if 'CFOUR_CALC_LEVEL' in options['CFOUR']:
        options['CFOUR']['CFOUR_CALC_LEVEL']['clobber'] = True
        options['CFOUR']['CFOUR_CALC_LEVEL']['superclobber'] = True
    if 'CFOUR_CC_PROGRAM' in options['CFOUR']:
        options['CFOUR']['CFOUR_CC_PROGRAM']['clobber'] = False

    return text, options


def cfour_list():
    """Return an array of Cfour methods with energies. Appended
    to procedures['energy'].

    """
    val = []
    val.append('cfour')
    val.append('c4-scf')
    val.append('c4-mp2')
    val.append('c4-mp3')
    val.append('c4-mp4(sdq)')
    val.append('c4-mp4')
    val.append('c4-cc2')
    val.append('c4-ccsd')
    val.append('c4-cc3')
    val.append('c4-ccsd(t)')
    val.append('c4-ccsdt')
    return val


def cfour_gradient_list():
    """Return an array of Cfour methods with analytical gradients.
    Appended to procedures['gradient'].

    """
    val = []
    val.append('cfour')
    val.append('c4-scf')
    val.append('c4-mp2')
    val.append('c4-mp3')
    val.append('c4-mp4(sdq)')
    val.append('c4-mp4')
    val.append('c4-cc2')
    val.append('c4-ccsd')
    val.append('c4-cc3')
    val.append('c4-ccsd(t)')
    val.append('c4-ccsdt')
    return val


def cfour_psivar_list():
    """Return a dict with keys of most Cfour methods and values of dicts
    with the PSI Variables returned by those methods. Used by cbs()
    wrapper to avoid unnecessary computations in compound methods.
    Result is appended to ``VARH``.

    """
    VARH = {}
    VARH['c4-scf'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY'}
    VARH['c4-mp2'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY'}
    VARH['c4-mp3'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                      'c4-mp2.5corl': 'MP2.5 CORRELATION ENERGY',
                        'c4-mp3corl': 'MP3 CORRELATION ENERGY'}
    VARH['c4-mp4(sdq)'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                      'c4-mp2.5corl': 'MP2.5 CORRELATION ENERGY',
                        'c4-mp3corl': 'MP3 CORRELATION ENERGY',
                   'c4-mp4(sdq)corl': 'MP4(SDQ) CORRELATION ENERGY'}
    VARH['c4-mp4'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                      'c4-mp2.5corl': 'MP2.5 CORRELATION ENERGY',
                        'c4-mp3corl': 'MP3 CORRELATION ENERGY',
                   'c4-mp4(sdq)corl': 'MP4(SDQ) CORRELATION ENERGY',
                        'c4-mp4corl': 'MP4(SDTQ) CORRELATION ENERGY'}
    VARH['c4-cc2'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                        'c4-cc2corl': 'CC2 CORRELATION ENERGY'}
    VARH['c4-ccsd'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                       'c4-ccsdcorl': 'CCSD CORRELATION ENERGY'}
    VARH['c4-cc3'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                        'c4-cc3corl': 'CC3 CORRELATION ENERGY'}
    VARH['c4-ccsd(t)'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                       'c4-ccsdcorl': 'CCSD CORRELATION ENERGY',
                    'c4-ccsd(t)corl': 'CCSD(T) CORRELATION ENERGY'}
    VARH['c4-ccsdt'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                       'c4-ccsdcorl': 'CCSD CORRELATION ENERGY',
                      'c4-ccsdtcorl': 'CCSDT CORRELATION ENERGY'}

    return VARH


#def backtransform(chgeMol, permMol, chgeGrad=None, chgeDip=None):
#def format_fjobarc(fje, fjelem, fjcoord, fjgrd, map, fjdip):
def format_fjobarc(energy, map, elem, coordinates, gradient, dipole):
    """Takes the key results from a gradient computation (*energy*,
    element Z list *elem*, *coordinates*, *gradient*,
    *dipole*, and atom ordering *map*) and writes a string *fja*
    that exactly mimics the contents of a Cfour FJOBARC file.

    """
    fja = 'TOTENERG\n'
    fja += '%15d%15d\n' % (struct.unpack("ii", struct.pack("d", energy)))
    fja += 'COORD\n'
    Nat = len(coordinates)
    flatcoord = []
    for at in range(Nat):
        for xyz in range(3):
            flatcoord.append(coordinates[map[at]][xyz])
    for idx in range(len(flatcoord)):
        if abs(flatcoord[idx]) < 1.0E-14:  # TODO
            flatcoord[idx] = 0.0
        fja += '%15d%15d' % (struct.unpack("ii", struct.pack("d", flatcoord[idx])))
        if idx % 2 == 1:
            fja += '\n'
    if len(flatcoord) % 2 == 1:
        fja += '\n'
    fja += 'MAP2ZMAT\n'
    for idx in range(Nat):
        fja += '%15d%15d' % (struct.unpack("ii", struct.pack("l", map[idx] + 1)))
        if idx % 2 == 1:
            fja += '\n'
    if Nat % 2 == 1:
        fja += '\n'
    fja += 'GRD FILE\n'
    fja += '%5d%20.10f\n' % (Nat, 0.0)
    for at in range(Nat):
        fja += '%20.10f%20.10f%20.10f%20.10f\n' % (elem[at], coordinates[at][0], coordinates[at][1], coordinates[at][2])
    for at in range(Nat):
        fja += '%20.10f%20.10f%20.10f%20.10f\n' % (elem[at], gradient[at][0], gradient[at][1], gradient[at][2])
    fja += 'DIPOL FILE\n'
    fja += '%20.10f%20.10f%20.10f\n' % (dipole[0], dipole[1], dipole[2])

    return fja


def backtransform(chgeMol, permMol, chgeGrad=None, chgeDip=None):
    """Here, *chgeMol* and *chgeGrd* need to be turned into the native Cfour
    orientation embodied by *permMol*. Currently for vpt2.

    """
    # Set up array reorientation object
    p4c4 = OrientMols(permMol, chgeMol)  # opposite than usual
    oriCoord = p4c4.transform_coordinates2(chgeMol)
    p4Elem = []
    for at in range(chgeMol.natom()):
        p4Elem.append(chgeMol.Z(at))
    oriElem = p4c4.transform_elementlist(p4Elem)
    oriElemMap = p4c4.Catommap

    oriGrad = None if chgeGrad is None else p4c4.transform_gradient(chgeGrad)
    oriDip = None if chgeDip is None else p4c4.transform_vector(chgeDip)

    if chgeGrad and chgeDip:
        return oriElemMap, oriElem, oriCoord, oriGrad, oriDip
    else:
        return oriElemMap, oriElem, oriCoord


#def backtransform_grad(p4Mol, c4Mol, p4Grd, p4Dip):
#    """Here, p4Mol and p4Grd need to be turned into the native Cfour
#    orientation embodied by c4Mol. Currently for vpt2.
#
#    """
#    # Set up array reorientation object
#    p4c4 = OrientMols(c4Mol, p4Mol)  # opposite than usual
#    oriCoord = p4c4.transform_coordinates2(p4Mol)
#    oriGrad = p4c4.transform_gradient(p4Grd)
#    p4Elem = []
#    for at in range(p4Mol.natom()):
#        p4Elem.append(p4Mol.Z(at))
#    oriElem = p4c4.transform_elementlist(p4Elem)
#    oriElemMap = p4c4.Catommap
#    oriDip = p4c4.transform_vector(p4Dip)
#
#    #print p4c4
#    #print '    <<<   Input C4 Mol   >>>'
#    #c4Mol.print_out()
#    #print '    <<<   Input P4 Mol   >>>'
#    #p4Mol.print_out()
#    #print '    <<<   Input P4 Grad   >>>'
#    #if p4Grd is not None:
#    #    for item in p4Grd:
#    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    #print '    <<<   Rotated P4 Coord   >>>'
#    #if oriCoord is not None:
#    #    for item in oriCoord:
#    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    #print '    <<<   Rotated P4 Elem   >>>'
#    #if oriElem is not None:
#    #    for item in oriElem :
#    #        print('       %16.8f' % (item))
#    #print '    <<<   Rotated P4 Dip  >>>'
#    #if oriDip is not None:
#    #    print('       %16.8f %16.8f %16.8f' % (oriDip[0], oriDip[1], oriDip[2]))
#    #print '    <<<   Rotated P4 Grad   >>>'
#    #if oriGrad is not None:
#    #    for item in oriGrad:
#    #        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#
#    return oriElemMap, oriElem, oriCoord, oriGrad, oriDip
#    #return oriElem, oriCoord, oriGrad, oriElemMap, oriDip


def jajo2mol(jajodic):
    """Returns a Molecule from entries in dictionary *jajodic* extracted
    from JAINDX and JOBARC.

    """
    map = jajodic['MAP2ZMAT']
    elem = jajodic['ATOMCHRG']
    coord = jajodic['COORD   ']
    Nat = len(elem)

    molxyz = '%d bohr\n\n' % (Nat)
    # TODO chgmult, though not really necessary for reorientation
    for at in range(Nat):
        posn = map[at] - 1
        el = 'GH' if elem[posn] == 0 else z2el[elem[posn]]
        posn *= 3
        molxyz += '%s %21.15f %21.15f %21.15f\n' % (el, coord[posn], coord[posn + 1], coord[posn + 2])
    mol = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)

    return mol
