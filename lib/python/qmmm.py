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

"""Module with classes to integrate MM charges into
a QM calculation.

"""
import psi4
import re
import os
import math
import p4const
from molutil import *
from driver import *


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
        basisChanged = psi4.has_option_changed("BASIS")
        ribasisChanged = psi4.has_option_changed("DF_BASIS_SCF")
        scftypeChanged = psi4.has_option_changed("SCF_TYPE")

        basis = psi4.get_option("BASIS")
        ribasis = psi4.get_option("DF_BASIS_SCF")
        scftype = psi4.get_option("SCF_TYPE")

        psi4.print_out("    => Diffuse SCF (Determines Da) <=\n\n")
        activate(self.molecule)

        psi4.set_global_option("BASIS", self.basisname)
        psi4.set_global_option("DF_BASIS_SCF", self.ribasisname)
        psi4.set_global_option("SCF_TYPE", "DF")
        energy('scf')
        psi4.print_out("\n")

        self.fitGeneral()

        psi4.clean()

        psi4.set_global_option("BASIS", basis)
        psi4.set_global_option("DF_BASIS_SCF", ribasis)
        psi4.set_global_option("SCF_TYPE", scftype)

        if not basisChanged:
            psi4.revoke_option_changed("BASIS")
        if not ribasisChanged:
            psi4.revoke_option_changed("DF_BASIS_SCF")
        if not scftypeChanged:
            psi4.revoke_option_changed("SCF_TYPE")

    def fitGeneral(self):
        """Function to perform a general fit of diffuse charges
        to wavefunction density.

        """
        psi4.print_out("    => Diffuse Charge Fitting (Determines da) <=\n\n")
        self.wfn = psi4.wavefunction()
        self.Da = self.wfn.Da()
        self.basis = self.wfn.basisset()
        parser = psi4.Gaussian94BasisSetParser()
        self.ribasis = psi4.BasisSet.construct(parser, self.molecule, "DF_BASIS_SCF")

        fitter = psi4.DFChargeFitter()
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
        self.extern = psi4.ExternalPotential()

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
        self.charges.append([Q, x / p4const.psi_bohr2angstroms, y / p4const.psi_bohr2angstroms, z / p4const.psi_bohr2angstroms])

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
