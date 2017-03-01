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

"""Module with classes to integrate MM charges into
a QM calculation.

"""
from __future__ import absolute_import
import re
import os
import math
from psi4.driver import *


class Diffuse(object):

    def __init__(self, molecule, basisname, ribasisname):

        self.molecule = molecule
        self.basisname = basisname
        self.ribasisname = ribasisname
        self.basis = None
        self.ribasis = None
        self.da = None
        self.Da = None
        self.wfn = None

    def __str__(self):

        s = '    => Diffuse <=\n\n'
        s = s + '    ' + str(self.molecule) + '\n'
        s = s + '    ' + self.basisname + '\n'
        s = s + '    ' + self.ribasisname + '\n'
        s = s + '\n'

        return s

    def fitScf(self):
        """Function to run scf and fit a system of diffuse charges to
        resulting density.

        """
        basisChanged = core.has_option_changed("BASIS")
        ribasisChanged = core.has_option_changed("DF_BASIS_SCF")
        scftypeChanged = core.has_option_changed("SCF_TYPE")

        basis = core.get_option("BASIS")
        ribasis = core.get_option("DF_BASIS_SCF")
        scftype = core.get_option("SCF_TYPE")

        core.print_out("    => Diffuse SCF (Determines Da) <=\n\n")

        core.set_global_option("BASIS", self.basisname)
        core.set_global_option("DF_BASIS_SCF", self.ribasisname)
        core.set_global_option("SCF_TYPE", "DF")
        E, ref = energy('scf', return_wfn=True, molecule=self.molecule)
        self.wfn = ref
        core.print_out("\n")

        self.fitGeneral()

        core.clean()

        core.set_global_option("BASIS", basis)
        core.set_global_option("DF_BASIS_SCF", ribasis)
        core.set_global_option("SCF_TYPE", scftype)

        if not basisChanged:
            core.revoke_option_changed("BASIS")
        if not ribasisChanged:
            core.revoke_option_changed("DF_BASIS_SCF")
        if not scftypeChanged:
            core.revoke_option_changed("SCF_TYPE")

    def fitGeneral(self):
        """Function to perform a general fit of diffuse charges
        to wavefunction density.

        """
        core.print_out("    => Diffuse Charge Fitting (Determines da) <=\n\n")
        self.Da = self.wfn.Da()
        self.basis = self.wfn.basisset()
        parser = core.Gaussian94BasisSetParser()
        self.ribasis = core.BasisSet.construct(parser, self.molecule, "DF_BASIS_SCF")

        fitter = core.DFChargeFitter()
        fitter.setPrimary(self.basis)
        fitter.setAuxiliary(self.ribasis)
        fitter.setD(self.Da)
        self.da = fitter.fit()
        self.da.scale(2.0)

    def populateExtern(self, extern):
        # Electronic Part
        extern.addBasis(self.ribasis, self.da)
        # Nuclear Part
        for A in range(0, self.molecule.natom()):
            extern.addCharge(self.molecule.Z(A), self.molecule.x(A), self.molecule.y(A), self.molecule.z(A))


class QMMM(object):

    def __init__(self):
        self.charges = []
        self.diffuses = []
        self.extern = core.ExternalPotential()

    def addDiffuse(self, diffuse):
        """Function to add a diffuse charge field *diffuse*."""
        self.diffuses.append(diffuse)

    def addChargeBohr(self, Q, x, y, z):
        """Function to add a point charge of magnitude *Q* at
        position (*x*, *y*, *z*) Bohr.

        """
        self.charges.append([Q, x, y, z])

    def addChargeAngstrom(self, Q, x, y, z):
        """Function to add a point charge of magnitude *Q* at
        position (*x*, *y*, *z*) Angstroms.

        """
        self.charges.append([Q, x / constants.bohr2angstroms, y / constants.bohr2angstroms, z / constants.bohr2angstroms])

    def __str__(self):

        s = '   ==> QMMM <==\n\n'

        s = s + '   => Charges (a.u.) <=\n\n'
        s = s + '    %11s %11s %11s %11s\n' % ('Z', 'x', 'y', 'z')
        for k in range(0, len(self.charges)):
            s = s + '    %11.7f %11.3E %11.3E %11.3E\n' % (self.charges[k][0], self.charges[k][1], self.charges[k][2], self.charges[k][3])
        s = s + '\n'

        s = s + '    => Diffuses <=\n\n'

        for k in range(0, len(self.diffuses)):
            s = s + str(self.diffuses[k])

        return s

    def populateExtern(self):
        """Function to define a charge field external to the
        molecule through point and diffuse charges.

        """
        # Charges
        for charge in self.charges:
            self.extern.addCharge(charge[0], charge[1], charge[2], charge[3])
        # Diffuses
        for diffuse in self.diffuses:
            diffuse.populateExtern(self.extern)
