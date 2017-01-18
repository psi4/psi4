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

"""Parent classes for quantum chemistry program input and output file
formats.
"""
from __future__ import absolute_import
from __future__ import print_function
import re


class InputFormat(object):

    def __init__(self, mem, mtd, bas, mol, sys, cast):

        # total job memory in MB
        self.memory = mem
        # computational method
        self.method = mtd.lower()
        # qcdb.Molecule object
        self.molecule = mol
        # database member index
        self.index = sys
        # orbital basis set
        self.basis = bas.lower()
        # do cast up from sto-3g basis?
        self.castup = cast

    def corresponding_aux_basis(self):
        """For Dunning basis sets, returns strings from which auxiliary
        basis sets and heavy-aug can be constructed. Note that
        valence/core-valence/etc. is conserved and X-zeta/(X+d)zeta is
        not, since this is the usual aux basis pattern.
        *augbasis* is round up to the nearest aug-cc-pVXZ
        *rootbasis* is round down to the nearest cc-pVXZ
        *auxbasis* is round up to the nearest cc-pVXZ or aug-cc-pVXZ
        """
        Dunmatch = re.compile(r'^(.*cc-)(pv|pcv|pwcv).*?([dtq56]).*z$').match(self.basis)

        if Dunmatch:
            rootbas = 'cc-' + Dunmatch.group(2) + Dunmatch.group(3) + 'z'
            augbas = 'aug-cc-' + Dunmatch.group(2) + Dunmatch.group(3) + 'z'
            if Dunmatch.group(1) == 'cc-':
                auxbas = rootbas
            else:
                auxbas = augbas
        else:
            rootbas = None
            augbas = None
            auxbas = None

        return [rootbas, augbas, auxbas]


class InputFormat2(object):

    def __init__(self, mem, mol, mtd, der, opt):

        # total job memory in MB
        self.memory = mem
        # qcdb.Molecule object
        self.molecule = mol
        # computational method
        self.method = mtd.lower()
        # computational derivative level
        self.dertype = der
        # options dictionary
        self.options = opt
        # orbital basis set
        self.basis = opt['GLOBALS']['BASIS']['value'].lower()
        # do cast up from sto-3g basis?
        self.castup = opt['SCF']['BASIS_GUESS']['value']

    def corresponding_aux_basis(self):
        """For Dunning basis sets, returns strings from which auxiliary
        basis sets and heavy-aug can be constructed. Note that
        valence/core-valence/etc. is conserved and X-zeta/(X+d)zeta is
        not, since this is the usual aux basis pattern.
        *augbasis* is round up to the nearest aug-cc-pVXZ
        *rootbasis* is round down to the nearest cc-pVXZ
        *auxbasis* is round up to the nearest cc-pVXZ or aug-cc-pVXZ
        """
        Dunmatch = re.compile(r'^(.*cc-)(pv|pcv|pwcv).*?([dtq56]).*z$').match(self.basis)

        if Dunmatch:
            rootbas = 'cc-' + Dunmatch.group(2) + Dunmatch.group(3) + 'z'
            augbas = 'aug-cc-' + Dunmatch.group(2) + Dunmatch.group(3) + 'z'
            if Dunmatch.group(1) == 'cc-':
                auxbas = rootbas
            else:
                auxbas = augbas
        else:
            rootbas = None
            augbas = None
            auxbas = None

        return [rootbas, augbas, auxbas]
