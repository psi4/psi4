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
import math

from .exceptions import *
from . import qcformat
from . import molpro_basissets


class MolproIn(qcformat.InputFormat):

    def __init__(self, mem, mtd, bas, mol, sys, cast):
        qcformat.InputFormat.__init__(self, mem, mtd, bas, mol, sys, cast)

        # memory in MB --> MW
        self.memory = int(math.ceil(mem / 8.0))
        # auxiliary basis sets
        [self.unaugbasis, self.augbasis, self.auxbasis] = self.corresponding_aux_basis()

    def format_global_parameters(self):
        text = ''

        if self.method in ['mp2c', 'dft-sapt-shift', 'dft-sapt', 'dft-sapt-pbe0ac', 'dft-sapt-pbe0acalda']:
            text += """GTHRESH,ZERO=1.e-14,ONEINT=1.e-14,TWOINT=1.e-14,ENERGY=1.e-8,ORBITAL=1.e-8,GRID=1.e-8\n\n"""
        elif self.method in ['b3lyp', 'b3lyp-d', 'df-b3lyp', 'df-b3lyp-d']:
            text += """GTHRESH,ZERO=1.e-14,ONEINT=1.e-14,TWOINT=1.e-14,ENERGY=1.e-8,ORBITAL=1.e-7,GRID=1.e-8\n\n"""
        else:
            text += """GTHRESH,ZERO=1.e-14,ONEINT=1.e-14,TWOINT=1.e-14,ENERGY=1.e-9\n\n"""

        return text

    def format_basis(self):
        text = ''
        text += """basis={\n"""

        try:
            # jaxz, maxz, etc.
            for line in molpro_basissets.altbasis[self.basis]:
                text += """%s\n""" % (line)
            text += '\n'
        except KeyError:
            # haxz
            if self.basis.startswith('heavy-aug-'):
                text += """set,orbital; default,%s,H=%s\n""" % (self.basis[6:], self.unaugbasis)
            # xz, axz, 6-31g*
            else:
                text += """set,orbital; default,%s\n""" % (self.basis)

        if ('df-' in self.method) or ('f12' in self.method) or (self.method in ['mp2c', 'dft-sapt', 'dft-sapt-pbe0acalda']):
            if self.unaugbasis and self.auxbasis:

                text += """set,jkfit;   default,%s/jkfit\n""" % (self.auxbasis)
                text += """set,jkfitb;  default,%s/jkfit\n""" % (self.unaugbasis)
                text += """set,mp2fit;  default,%s/mp2fit\n""" % (self.auxbasis)
                text += """set,dflhf;   default,%s/jkfit\n""" % (self.auxbasis)
            else:
                raise ValidationError("""Auxiliary basis not predictable from orbital basis '%s'""" % (self.basis))

        text += """}\n\n"""
        return text

    def format_infile_string(self):
        text = ''

        # format comment and memory
        text += """***, %s %s\n""" % (self.index, self.molecule.tagline)
        text += """memory,%d,m\n""" % (self.memory)

        # format molecule, incl. charges and dummy atoms
        text += self.molecule.format_molecule_for_molpro()

        # format global convergence directions
        text += self.format_global_parameters()

        # format castup directions
        if self.castup is True:
            text += """basis=sto-3g\n"""
            text += """rhf\n"""
            text += '\n'

        # format basis set
        text += self.format_basis()

        # format method
        for line in qcmtdIN[self.method]:
            text += """%s\n""" % (line)
        text += """show[1,20f20.12],ee*,ce*,te*\n"""
        text += """show[1,60f20.12],_E*\n"""
        text += '\n'

        return text


qcmtdIN = {
'ccsd(t)-f12': [
    'rhf',
    'eehf=energy',
    'ccsd(t)-f12,df_basis=mp2fit,df_basis_exch=jkfitb,ri_basis=jkfitb',
    'eemp2=emp2',
    'cemp2=eemp2-eehf',
    'eemp3=emp3',
    'cemp3=eemp3-eehf',
    'eeccsd=energc',
    'ceccsd=eeccsd-eehf',
    'eeccsdt=energy',
    'ceccsdt=eeccsdt-eehf',
    'temp2=emp2_trip',
    'teccsd=ectrip'],

'ccsd(t)': [
    'rhf',
    'eehf=energy',
    'ccsd(t)',
    'eemp2=emp2',
    'cemp2=eemp2-eehf',
    'eemp3=emp3',
    'cemp3=eemp3-eehf',
    'eeccsd=energc',
    'ceccsd=eeccsd-eehf',
    'eeccsdt=energy',
    'ceccsdt=eeccsdt-eehf',
    'temp2=emp2_trip',
    'teccsd=ectrip'],

'mp3': [
    'gdirect',
    'rhf',
    'eehf=energy',
    'mp3',
    'eemp2=emp2',
    'eemp3=emp3',
    'eemp25=0.5*(eemp2+eemp3)',
    'cemp2=eemp2-eehf',
    'cemp3=eemp3-eehf',
    'cemp25=eemp25-eehf',
    'temp2=emp2_trip',
    'temp3=ectrip'],

'mp2': [
    'gdirect',
    'rhf',
    'eehf=energy',
    'mp2',
    'eemp2=emp2',
    'cemp2=eemp2-eehf',
    'temp2=emp2_trip'],

'df-hf-mp2': [
    'gdirect',
    '{df-hf,basis=jkfit}',
    'eehf=energy',
    'mp2',
    'eemp2=emp2',
    'cemp2=eemp2-eehf',
    'temp2=emp2_trip'],

'hf-df-mp2': [
    'gdirect',
    'rhf',
    'eehf=energy',
    '{df-mp2,basis_mp2=mp2fit}',
    'eemp2=emp2',
    'cemp2=eemp2-eehf',
    'temp2=emp2_trip'],

'hf': [
    'rhf',
    'eehf=energy'],

'mp2-f12': [
    'gdirect',
    'rhf',
    'eehf=energy',
    'mp2-f12',
    'eemp2=emp2',
    'cemp2=eemp2-eehf',
    'temp2=emp2_trip'],

'df-mp2-f12': [
    'gdirect',
    #'rhf',
    '{df-hf,basis=jkfit}',
    'eehf=energy',
    #'{df-mp2-f12,df_basis=mp2fit,df_basis_exch=jkfit,ri_basis=optrib}',
    '{df-mp2-f12,df_basis=mp2fit,df_basis_exch=jkfitb,ri_basis=jkfitb}',
    'eemp2=emp2',
    'cemp2=eemp2-eehf',
    'temp2=emp2_trip'],

'df-mp2': [
    'gdirect',
    '{df-hf,basis=jkfit}',
    'eehf=energy',
    '{df-mp2,basis_mp2=mp2fit}',
    'eemp2=emp2',
    'cemp2=eemp2-eehf',
    'temp2=emp2_trip'],

'df-hf': [
    'gdirect',
    '{df-hf,basis=jkfit}',
    'eehf=energy'],

'b3lyp-d': [
    'gdirect',
    'rks,b3lyp3',
    'eehf=energy',
    'dispcorr',
    'eehfd=eehf+edisp'],

'df-b3lyp-d': [
    'gdirect',
    '{df-rks,b3lyp3,basis=jkfit}',
    'eehf=energy',
    'dispcorr',
    'eehfd=eehf+edisp'],

'b3lyp': [
    'gdirect',
    'rks,b3lyp3',
    'eehf=energy'],

'df-b3lyp': [
    'gdirect',
    '{df-rks,b3lyp3,basis=jkfit}',
    'eehf=energy'],

#'mp2c': [ # this job computes one part [E_disp(TDDFT)] of the three parts of a MP2C calculation
#        # check that nfrag = 2
#         'gdirect',
#         'ga=1101.2; gb=1102.2',
#         'ca=2101.2; cb=2102.2\n',
#
#         $spin = $cgmp{MLPmol1} - 1;
#         'SET,CHARGE=$cgmp{CHGmol1}',
#         'SET,SPIN=$spin',
#         'dummy',
#         foreach $at (@monoBreal) { print $handle ",$at"; }
#         ''
#         '{df-hf,basis=jkfit,locorb=0; start,atdens; save,$ga}',
#         '{df-ks,lhf,df_basis=dflhf,basis_coul=jkfitb,basis_exch=jkfitb; dftfac,1.0; start,$ga; save,$ca}',
#         'eehfa=energy; sapt; monomerA',
#         '',
#
#         $spin = $cgmp{MLPmol2} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGmol2}\nSET,SPIN=$spin\ndummy";
#         foreach $at (@monoAreal) { print $handle ",$at"; }
#         print $handle "\n{df-hf,basis=jkfit,locorb=0; start,atdens; save,\$gb}\n";
#         print $handle "{df-ks,lhf,df_basis=dflhf,basis_coul=jkfitb,basis_exch=jkfitb; dftfac,1.0; start,\$gb; save,\$cb}\n";
#         print $handle "eehfb=energy; sapt; monomerB\n\n";
#
#         $spin = $cgmp{MLPsyst} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGsyst}\nSET,SPIN=$spin\n";
#         print $handle "{sapt,SAPT_LEVEL=3; intermol,ca=\$ca,cb=\$cb,icpks=0,fitlevel=3,nlexfac=0.0,cfac=0.0\n";
#         print $handle "dfit,basis_coul=jkfit,basis_exch=jkfit,cfit_scf=3}\n";
#         print $handle "eedisp=E2disp\n\n";
#
#    ],
}

#'dft-sapt-shift': [
#
#         # this is written in an inflexible way (fixed basis, functional) so that it is computed
#         #  only once, then used when writing DFT-SAPT inputs, which we'll be more flexible with
#
#         print $handle "basis={\n";
#         print $handle "set,orbital; default,aug-cc-pVQZ\n";
#         print $handle "set,jkfit;   default,avqz/jkfit\n";
#         print $handle "set,dflhf;   default,avqz/jkfit\n";
#         print $handle "}\n";
#
#         if    ($handle eq "M1OUT") { $charge = $cgmp{CHGmol1}; $spin = $cgmp{MLPmol1} - 1; }
#         elsif ($handle eq "M2OUT") { $charge = $cgmp{CHGmol2}; $spin = $cgmp{MLPmol2} - 1; }
#
#         print $handle "\ngdirect\n";
#         print $handle "{df-ks,pbex,pw91c,lhf; dftfac,0.75,1.0,0.25}\n";
#         print $handle "basis=tzvpp\n";
#         print $handle "{ks,pbe0; orbprint,0}\n";
#         print $handle "eeneut=energy\n";
#         $charge += 1;
#         $spin += 1;
#         print $handle "SET,CHARGE=$charge\nSET,SPIN=$spin\n";
#         print $handle "{ks,pbe0}\n";
#         print $handle "eecat=energy\n";
#         print $handle "eeie=eecat-eeneut\n";
#         print $handle "show[1,20f20.12],ee*,ce*,te*\n";
#         print $handle "show[1,60f20.12],_E*\n";
#    ]
#'dft-sapt': [
#
#         if ( ($asyA eq '') || ($asyB eq '') ) {
#            print "ERROR: asymptotic correction not defined for one or more monomers in index $system.\n";
#            close(DIOUT);
#            unlink("$pathDIOUT");
#         }
#
#         print $handle "gdirect\n";
#         print $handle "ca=2101.2; cb=2102.2\n\n";
#
#         $spin = $cgmp{MLPmol1} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGmol1}\nSET,SPIN=$spin\ndummy";
#         foreach $at (@monoBreal) { print $handle ",$at"; }
#         print $handle "\n{df-ks,pbex,pw91c,lhf,df_basis=dflhf,basis_coul=jkfitb,basis_exch=jkfitb; dftfac,0.75,1.0,0.25; asymp,$asyA; save,\$ca}\n";
#         print $handle "eehfa=energy; sapt; monomerA\n\n";
#
#         $spin = $cgmp{MLPmol2} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGmol2}\nSET,SPIN=$spin\ndummy";
#         foreach $at (@monoAreal) { print $handle ",$at"; }
#         print $handle "\n{df-ks,pbex,pw91c,lhf,df_basis=dflhf,basis_coul=jkfitb,basis_exch=jkfitb; dftfac,0.75,1.0,0.25; asymp,$asyB; save,\$cb}\n";
#         print $handle "eehfb=energy; sapt; monomerB\n\n";
#
#         $spin = $cgmp{MLPsyst} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGsyst}\nSET,SPIN=$spin\n";
#         print $handle "{sapt,sapt_level=3; intermol,ca=\$ca,cb=\$cb,icpks=0,fitlevel=3,nlexfac=0.0\n";
#         print $handle "dfit,basis_coul=jkfit,basis_exch=jkfit,basis_mp2=mp2fit,cfit_scf=3}\n";
#         print $handle "eeelst=E1pol\n";
#         print $handle "eeexch=E1ex\n";
#         print $handle "eeind=E2ind\n";
#         print $handle "eeexind=E2exind\n";
#         print $handle "eedisp=E2disp\n";
#         print $handle "eeexdisp=E2exdisp\n\n";
#
#    ]
#'dft-sapt-pbe0ac': [
#
#         if ( ($asyA eq '') || ($asyB eq '') ) {
#            print "ERROR: asymptotic correction not defined for one or more monomers in index $system.\n";
#            close(DIOUT);
#            unlink("$pathDIOUT");
#         }
#
#         print $handle "ca=2101.2; cb=2102.2\n\n";
#
#         $spin = $cgmp{MLPmol1} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGmol1}\nSET,SPIN=$spin\ndummy";
#         foreach $at (@monoBreal) { print $handle ",$at"; }
#         print $handle "\n{ks,pbe0; asymp,$asyA; save,\$ca}\n";
#         print $handle "eehfa=energy; sapt; monomerA\n\n";
#
#         $spin = $cgmp{MLPmol2} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGmol2}\nSET,SPIN=$spin\ndummy";
#         foreach $at (@monoAreal) { print $handle ",$at"; }
#         print $handle "\n{ks,pbe0; asymp,$asyB; save,\$cb}\n";
#         print $handle "eehfb=energy; sapt; monomerB\n\n";
#
#         $spin = $cgmp{MLPsyst} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGsyst}\nSET,SPIN=$spin\n";
#         print $handle "{sapt; intermol,ca=\$ca,cb=\$cb,icpks=0}\n";
#         print $handle "eeelst=E1pol\n";
#         print $handle "eeexch=E1ex\n";
#         print $handle "eeind=E2ind\n";
#         print $handle "eeexind=E2exind\n";
#         print $handle "eedisp=E2disp\n";
#         print $handle "eeexdisp=E2exdisp\n\n";
#    ]
#'dft-sapt-pbe0acalda': [
#
#         if ( ($asyA eq '') || ($asyB eq '') ) {
#            print "ERROR: asymptotic correction not defined for one or more monomers in index $system.\n";
#            close(DIOUT);
#            unlink("$pathDIOUT");
#         }
#
#         print $handle "ca=2101.2; cb=2102.2\n\n";
#
#         $spin = $cgmp{MLPmol1} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGmol1}\nSET,SPIN=$spin\ndummy";
#         foreach $at (@monoBreal) { print $handle ",$at"; }
#         print $handle "\n{ks,pbe0; asymp,$asyA; save,\$ca}\n";
#         print $handle "eehfa=energy; sapt; monomerA\n\n";
#
#         $spin = $cgmp{MLPmol2} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGmol2}\nSET,SPIN=$spin\ndummy";
#         foreach $at (@monoAreal) { print $handle ",$at"; }
#         print $handle "\n{ks,pbe0; asymp,$asyB; save,\$cb}\n";
#         print $handle "eehfb=energy; sapt; monomerB\n\n";
#
#         $spin = $cgmp{MLPsyst} - 1;
#         print $handle "SET,CHARGE=$cgmp{CHGsyst}\nSET,SPIN=$spin\n";
#         print $handle "{sapt,sapt_level=3; intermol,ca=\$ca,cb=\$cb,icpks=0,fitlevel=3,nlexfac=0.0\n";
#         print $handle "dfit,basis_coul=jkfit,basis_exch=jkfit,basis_mp2=mp2fit,cfit_scf=3}\n";
#         print $handle "eeelst=E1pol\n";
#         print $handle "eeexch=E1ex\n";
#         print $handle "eeind=E2ind\n";
#         print $handle "eeexind=E2exind\n";
#         print $handle "eedisp=E2disp\n";
#         print $handle "eeexdisp=E2exdisp\n\n";
#
#         print $handle "show[1,20f20.12],ee*,ce*,te*\n";
#         print $handle "show[1,60f20.12],_E*\n";
#      }
#
