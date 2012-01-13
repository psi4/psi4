import PsiMod
import re
import os
import input
import math
from molecule import * 
from driver import * 

class Diffuse:

    def __init__(self,molecule,basisname,ribasisname):

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

        PsiMod.print_out("    => Diffuse SCF (Determines Da) <=\n\n")
        activate(self.molecule)
        PsiMod.set_global_option("BASIS", self.basisname) 
        PsiMod.set_global_option("RI_BASIS_SCF", self.ribasisname) 
        PsiMod.set_global_option("SCF_TYPE", "DF") 
        energy('scf') 
        PsiMod.print_out("\n")

        self.fitGeneral()

    def fitGeneral(self):

        PsiMod.print_out("    => Diffuse Charge Fitting (Determines da) <=\n\n");        
        self.wfn = PsiMod.reference_wavefunction() 
        self.Da  = self.wfn.Da() 
        self.basis = self.wfn.basisset()
        parser = PsiMod.Gaussian94BasisSetParser()
        self.ribasis = PsiMod.BasisSet.construct(parser,self.molecule,"RI_BASIS_SCF")

        fitter = PsiMod.DFChargeFitter()
        fitter.setPrimary(self.basis)
        fitter.setAuxiliary(self.ribasis)
        fitter.setD(self.Da)
        self.da = fitter.fit()

    def populateExtern(self, extern):
        # Electronic Part
        extern.addBasis(self.ribasis,self.da)
        # Nuclear Part
        for A in range(0,self.molecule.natom()):
            extern.addCharge(self.molecule.Z(A),self.molecule.x(A),self.molecule.y(A),self.molecule.z(A))
        
class QMMM:

    def __init__(self):
        self.charges = []
        self.diffuses = []
     
    def addDiffuse(self, diffuse):
        self.diffuses.append(diffuse)

    def addChargeBohr(self, Q, x, y, z):
        self.charges.append([Q,x,y,z])

    def addChargeAngstrom(self, Q, x, y, z):
        self.charges.append([Q,x/0.529177249,y/0.529177249,z/0.529177249])

    def __str__(self):
        
        s = '   ==> QMMM <==\n\n';

        s = s + '   => Charges (a.u.) <=\n\n';
        s = s + '    %11s %11s %11s %11s\n' % ('Z', 'x', 'y', 'z')
        for k in range(0,len(self.charges)):
            s = s + '    %11.7f %11.3E %11.3E %11.3E\n' % (self.charges[k][0], self.charges[k][1], self.charges[k][2], self.charges[k][3])
        s = s + '\n'

        s = s + '    => Diffuses <=\n\n';

        for k in range(0,len(self.diffuses)):
            s = s + str(self.diffuses[k])

        return s 
    
    def populateExtern(self, extern=None):

        if (extern == None):
            extern = PsiMod.ExternalPotential()
            extern.setName("QM/MM Field")
            activate_potential(extern)
    
        # Charges
        for charge in self.charges:
            extern.addCharge(charge[0],charge[1],charge[2],charge[3])
        # Diffuses
        for diffuse in self.diffuses:
            diffuse.populateExtern(extern)

        extern.print_out()
