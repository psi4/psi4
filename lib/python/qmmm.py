"""Module with classes to integrate MM charges into
a QM calculation.

"""
import PsiMod
import re
import os
import input
import math
import physconst
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
        basisChanged = PsiMod.has_option_changed("BASIS")
        ribasisChanged = PsiMod.has_option_changed("DF_BASIS_SCF")
        scftypeChanged = PsiMod.has_option_changed("SCF_TYPE")

        basis = PsiMod.get_option("BASIS")
        ribasis = PsiMod.get_option("DF_BASIS_SCF")
        scftype = PsiMod.get_option("SCF_TYPE")

        PsiMod.print_out("    => Diffuse SCF (Determines Da) <=\n\n")
        activate(self.molecule)

        PsiMod.set_global_option("BASIS", self.basisname)
        PsiMod.set_global_option("DF_BASIS_SCF", self.ribasisname)
        PsiMod.set_global_option("SCF_TYPE", "DF")
        energy('scf')
        PsiMod.print_out("\n")

        self.fitGeneral()

        PsiMod.clean()

        PsiMod.set_global_option("BASIS", basis)
        PsiMod.set_global_option("DF_BASIS_SCF", ribasis)
        PsiMod.set_global_option("SCF_TYPE", scftype)

        if not basisChanged:
            PsiMod.revoke_option_changed("BASIS")
        if not ribasisChanged:
            PsiMod.revoke_option_changed("DF_BASIS_SCF")
        if not scftypeChanged:
            PsiMod.revoke_option_changed("SCF_TYPE")

    def fitGeneral(self):
        """Function to perform a general fit of diffuse charges
        to wavefunction density.

        """
        PsiMod.print_out("    => Diffuse Charge Fitting (Determines da) <=\n\n")
        self.wfn = PsiMod.reference_wavefunction()
        self.Da = self.wfn.Da()
        self.basis = self.wfn.basisset()
        parser = PsiMod.Gaussian94BasisSetParser()
        self.ribasis = PsiMod.BasisSet.construct(parser, self.molecule, "DF_BASIS_SCF")

        fitter = PsiMod.DFChargeFitter()
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
        self.extern = PsiMod.ExternalPotential()

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
        self.charges.append([Q, x / physconst.psi_bohr2angstroms, y / physconst.psi_bohr2angstroms, z / physconst.psi_bohr2angstroms])

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
