#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import absolute_import

#
# @BEGIN LICENSE
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
# @END LICENSE
#

import sys
import re
import os
import math
import warnings

# CUimport psi4
# CUimport p4const
# CUimport p4util

from driver import *

# from extend_Molecule import *
# CUfrom molutil import *
# CUfrom p4regex import *
# never import aliases into this file

# _____________________________________________________________________________
#
# run_g3mp2 is a Psi4 implementation of Curtiss' G3(MP2) thermochemical method
# as described in:
#
#  Gaussian-3 theory using reduced Moeller-Plesset order
#  Curtiss, Redfern, Raghavachari,Rassolov, Pople
#  J.Chem.Phys., 110, 10(1999)
#
# The CCSD(T) and HLC revisions are in accord with the later article:
#
#  Gaussian-3 theory using coupled cluster energies
#  Curtiss, Raghavachari, Redfern, Baboul, Pople
#  Chemical Physics Letters 314(1999) 101-107
#
# Daniel R. Haney 3/28/2016 --------------------------
#
# ------------------------------------------------------------------------
#   TODO:
#       Add memory check on MP(2,fc)/G3MP2Large calc
#       to select for OUT_OF_CORE when it may exhaust memory.
#
#       Add QCISD(T) option for Gaussian compatibility.
#

ZPE_SCALE_FACTOR = 0.8929


# explicitly:  Boltzmann kT contribution @ 298.15K in Hartrees
#    = (k_boltzmann * T298.15 * Avogadro) / JoulePerkcal / kcalPerHartree / cal_per_kcal

kT = (psi_kb * 298.15 * psi_na) / psi_cal2J / psi_hartree2kcalmol / 1000.0

kT_kcal = (psi_kb * 298.15 * psi_na) / psi_cal2J



#
# Zero Point Energy scale factors vary with applications and vibration data.
#
# "Computational Thermochemistry: Scale Factor Databases and Scale Factors
# for Vibrational Frequencies Obtained from Electronic Model Chemistries,"
#  I. M. Alecu, Jingjing Zheng, Yan Zhao, and Donald G. Truhlar*
#  J. Chem. Theory Comput. 2010, 6, 2872.
#  DOI: 10.1021/ct100326h
#
# NIST's CCCBDB, using vibrations from ~270 compounds and recalculation
# excluding outliers, reports the scale factor as 0.899 +/- 0.025
#
# Curtiss' scale factor for HF/6-31G* ZPEs produced by Gaussian-94
# is overly precise at 4 sig fig, but is left unchanged since this module
# is a method replication and not a revision.
#

g3mp2_verbose_flag = False

def run_g3mp2(name, qcisd=False, **kwargs):

    def say(s=''):
        if g3mp2_verbose_flag:
            sys.stderr.write(s)  # print to console stderr without newline

    def log(line):
        psi4.print_out(line + '\n')
        if g3mp2_verbose_flag:
            print line

    def debug(s):
        return
        print s

    def dumpobj(obj, objname='obj'):
        print '%s is type %s' % (objname, type(obj).__name__)

        for attr in dir(obj):
            print '%s.%s = %s' % (objname, attr, getattr(obj, attr))

    #________________________________________________________________________
    def nfrozen_core ():
        _nfcz = 0
        _mol = psi4.get_active_molecule()
        natoms = _mol.natom()
        for A in range(0,natoms):
            if _mol.Z(A) > 2:
                _nfcz += 1
            if _mol.Z(A) > 10:
                _nfcz += 4
            if _mol.Z(A) > 18:
                _nfcz += 4
            if _mol.Z(A) > 36:
                _nfcz += 9
            if _mol.Z(A) > 54:
                _nfcz += 9
            if _mol.Z(A) > 86:
                _nfcz += 16
            if _mol.Z(A) > 108:
                raise ValidationError("run_g3mp2.nfrozen_core: Invalid atomic number")
        return _nfcz

    #________________________________________________________________________
    def E_spin_orbit(Z=0, charge=0):
        ''' spin orbit energy corrections for atomic species.
            See Curtiss ref. above.

            returns: values in Hartrees
        '''
        #       [neutral, ion+,     ion- ]
        _ESO = [[ 0.0,   0.0,    0.0 ],     # 00 index place holder
            [ 0.0,   0.0,    0.0 ],[ 0.0,   0.0,    0.0 ],  # 01 H   02 He
            [ 0.0,   0.0,    0.0 ],[ 0.0,   0.0,    0.0 ],  # 03 Li  04 Be
            [-0.05,  0.0,   -0.03],[-0.14, -0.2,    0.0 ],  # 05 B   06 C
            [ 0.0,  -0.43,   0.0 ],[-0.36,  0.0,   -0.26],  # 07 N   08 O
            [-0.61, -0.67,   0.0 ],[ 0.0,  -1.19,   0.0 ],  # 09 F   10 Ne
            [ 0.0,   0.0,    0.0 ],[ 0.0,   0.0,    0.0 ],  # 11 Na  12 Mg
            [-0.34,  0.0,   -0.28],[-0.68, -0.93,   0.0 ],  # 13 Al  14 Si
            [ 0.0,  -1.43,  -0.45],[-0.89,  0.0,   -0.88],  # 15 P   16 S
            [-1.34, -1.68,   0.0 ],[ 0.0,  -2.18,   0.0 ],  # 17 Cl  18 Ar
            [ 0.0,   0.0,    0.0 ],[ 0.0,   0.0,    0.0 ],  # 19 K   20 Ca
            [ 0.0,   0.0,    0.0 ],[ 0.0,   0.0,    0.0 ],  # 21 Sc  22 Ti
            [ 0.0,   0.0,    0.0 ],[ 0.0,   0.0,    0.0 ],  # 23 V   24 Cr
            [ 0.0,   0.0,    0.0 ],[ 0.0,   0.0,    0.0 ],  # 25 Mn  26 Fe
            [ 0.0,   0.0,    0.0 ],[ 0.0,   0.0,    0.0 ],  # 27 Co  28 Ni
            [ 0.0,   0.0,    0.0 ],[ 0.0,   0.0,    0.0 ],  # 29 Cu  30 Zn
            [-2.51,  0.0,    0.0 ],[-4.41, -5.37,   0.0 ],  # 31 Ga  32 Ge
            [ 0.0,  -8.04,   0.0 ],[-4.3,   0.0,    0.0 ],  # 33 As  34 Se
            [-5.6,  -6.71,   0.0 ],[ 0.0,  -8.16,   0.0 ]   # 35 Br  36 Kr
        ]

        if not Z in range(len(_ESO)):
            return 0.0

        if charge > 0:
            ion = 1
        elif charge < 0:
            ion = 2
        else:
            ion = 0

        debug("Spin orbit energy = %.2f mHa" % (_ESO[Z][ion]))

        milliHa_to_Ha = 0.001
        return (_ESO[Z][ion] * milliHa_to_Ha)
    #________________________________________________________________________
    def E0_ccsdt (Z=0):
        # G3(MP2) @0K energies of atomic species.
        #
        # returns:  atom E0 in Hartrees
        #
        # TODO: 3rd row
        debug("E0_ccsdt(%d)" % Z)

        _Z2E0 = [0.0,       # 0 index placeholder
                -0.501765, -2.902353, -7.366280, -14.629072,        # H  He Li Be
                -24.606789, -37.788989, -54.524780, -74.989201,     # B  C  N  O
                -99.640199, -128.827752, -161.847930, -199.650655,  # F  Ne Na Mg
                -241.936660, -288.939067, -340.826225, -397.663216, # Al Si P  S
                -459.686583, -527.060194, -599.160438, -676.790288] # Cl Ar K  Ca

        if Z < 1 and Z > len(_Z2E0):
            return 0.0
        else:
            return _Z2E0[Z]
    #________________________________________________________________________
    def E0_one_atom (Z=0):

        _e0 = E0_ccsdt(Z)
        _h298 = _e0 + (5.0/2 * kT)      # Ideal gas kinetic energy contribution

        debug("E0_one_atom: 5/2*kT = %.6f Ha" % (5.0/2 * kT))
        debug("E0_one_atom: e0 = %.6f h298 = %.6f" % (_e0, _h298))

        return (_e0, _h298)
    #________________________________________________________________________
    def dHf_one_atom(Z=0):
        '''return experimental atomic heats of formation
           returns:      tuple
           return type:     2 doubles
        '''
    #      atom [dHf(0), dHf(298)]  in kcal/mol

        atomDHF = [
            [0.0,0.0],            # 00 zero placeholder
            [51.63,   52.103  ],  # 01  Hydrogen
            [0.00,    0.00    ],  # 02  Helium
            [37.70,   38.07   ],  # 03  Lithium
            [76.40,   77.40   ],  # 04  Beryllium
            [135.10,  136.30  ],  # 05  Boron
            [169.98,  171.29  ],  # 06  Carbon
            [112.53,  112.97  ],  # 07  Nitrogen
            [58.99,   59.56   ],  # 08  Oxygen
            [18.47,   18.97   ],  # 09  Fluorine
            [0.00,    0.00    ],  # 10  Neon
            [25.76,   25.69   ],  # 11  Sodium
            [34.87,   35.16   ],  # 12  Magnesium
            [80.20,   80.80   ],  # 13  Aluminum
            [107.20,  108.20  ],  # 14  Silicon
            [75.45,   75.65   ],  # 15  Phosphorus
            [65.71,   66.25   ],  # 16  Sulfur
            [28.59,   28.99   ],  # 17  Chlorine
            [0.00,    0.00    ],  # 18  Argon
            [21.27,   21.49   ],  # 19  Potassium
            [42.50,   42.29   ],  # 20  Calcium
            [0.0,0.0],[0.0,0.0],  # transition elements 21-30 Sc-Zn
            [0.0,0.0],[0.0,0.0],  # transition elements
            [0.0,0.0],[0.0,0.0],  # transition elements
            [0.0,0.0],[0.0,0.0],  # transition elements
            [0.0,0.0],[0.0,0.0],  # transition elements
            [65.00,   65.00   ],  # 31  Gallium
            [88.91,   88.91   ],  # 32  Germanium
            [73.90,   72.42   ],  # 33  Arsenic
            [55.76,   54.27   ],  # 34  Selenium
            [26.74,   28.18   ],  # 35  Bromine
            [0.0,     0.0     ],  # 36  Krypton
        ]

        debug('dHf_one_atom(Z = %d)' % Z)
        if Z < len(atomDHF):
            debug('dHf_one_atom: E,H=%.2f,%.2f' %
                  (atomDHF[Z][0], atomDHF[Z][1]))
            return atomDHF[Z][0], atomDHF[Z][1]
        else:
            debug('dHf_one_atom: error: element %d not in table?' % Z)
            return 0.0, 0.0

    #_______________________________________________________
    # calculate heats of formation at 0K and 298K by atomization method
    #
    # input:    ab initio E0, H298 in Hartrees
    # returns: dHf@0K, dHf@298K -> 2 doubles

    def calc_deltaHf(_E0, _H298):

        debug("calc_deltaHf(%.6f, %.6f) entry" % (_E0, _H298))

        E0_sum_atoms = 0.0
        H298_sum_atoms = 0.0
        dhf0_sum_atoms = 0.0
        dhf298_sum_atoms = 0.0

        _natoms = mol.natom()

        for A in range(0,_natoms):
            Z = int(mol.Z(A))

            debug("\ncalc_deltaHf(Z=%d)" % (Z))

            e0_atom, h298_atom = E0_one_atom(Z)     # in Hartrees
            E0_sum_atoms += e0_atom
            H298_sum_atoms += h298_atom
            debug("E0_sum_atoms = %.6f Ha  H298_sum_atoms = %.6f Ha" % (E0_sum_atoms, H298_sum_atoms))

            d0, d298 = dHf_one_atom(Z)              # in kcal/mol
            dhf0_sum_atoms += d0
            dhf298_sum_atoms += d298
            debug('sumDHF0,sumDHF298 = %.2f kcal, %.2f kcal' % (dhf0_sum_atoms, dhf298_sum_atoms))

        debug("total E0 = %.6f Ha  E0_atoms = %.6f Ha" % (_E0, E0_sum_atoms))
        dhf0 = (_E0 - E0_sum_atoms) * psi_hartree2kcalmol + dhf0_sum_atoms

        debug("total H298 = %.6f Ha  H298_atoms = %.6f Ha" % (_H298, H298_sum_atoms))
        dhf298 = (_H298 - H298_sum_atoms) * psi_hartree2kcalmol + dhf298_sum_atoms

        #Eatomization = -(H298 - H298_sum_atoms) * psi_hartree2kcalmol

        debug('dhf0,dhf298 = %.2f,%.2f' % (dhf0, dhf298))

        return dhf0,dhf298
    #______________________________________________________


    use_QCISDT = False
    nAtoms = 0
    multiplicity = 1
    charge = 0

    kwargs = p4util.kwargs_lower(kwargs)
    g3mp2_verbose_flag = kwargs.pop('verbose', True)

    # debug('%s: %d args: %s' % (name, len(kwargs), ','.join([arg for arg in kwargs])))
    #
    # mol = kwargs.pop('molecule', psi4.get_active_molecule())

    mol = psi4.get_active_molecule()

    nAtoms = mol.natom()
    multiplicity = mol.multiplicity()
    charge = mol.molecular_charge()

    debug("charge = %d multiplicity=%d." % (charge, multiplicity))

    atomZlist = ",".join([str(int(mol.Z(i))) for i in range(nAtoms)])
    debug("There are %d atoms: %s" % (nAtoms, atomZlist))

    debug("frozen core count# %d" % (nfrozen_core()))

    #if 'qcisdt' in kwargs:
    #    use_QCISDT = kwargs.pop('qcisdt')
    #else:
    #   use_QCISDT = False

    # debug('use_QCISDT is %s' % str(use_QCISDT))

    if multiplicity > 1:
        psi4.set_global_option('REFERENCE', 'UHF')

    else:
        psi4.set_global_option('REFERENCE', 'RHF')

    # stash user options:
    optstash = p4util.OptionsState(['FREEZE_CORE'],
                                   ['MP2_TYPE'],
                                   ['SCF', 'SCF_TYPE'])

    # override default scf_type
    # OUT-OF-CORE, DIRECT, PK, DF available
    # DIRECT yields Eopt same as GAMESS, Nwchem to +/- 1.0E-06 Hartree for <= 6 big atoms
    # PK, a faster DIRECT
    # fastest: DF diverges beyond 1.0E-04 Ha ~= 1/15 kcal
    # OUT-OF-CORE: life is too short
    #psi4.set_local_option('SCF', 'SCF_TYPE', 'DIRECT')
    psi4.set_local_option('SCF', 'SCF_TYPE', 'PK')

    # optimize geometry at RHF/6-31G*
    say('HF optimize.')

    if nAtoms > 1:
        psi4.clean()
        psi4.set_global_option('BASIS', '6-31G(D)')
        optimize('scf')
        psi4.clean()

    else:
        pass

    # scf frequencies for ZPE
    say('ZPE.')

    if nAtoms == 1:
        # Set up Gas constant to Hartree/K units for obviousness.
        #temperature = 298.15        # TODO: fetch from molecule class instance
        #Rgas = psi_R / psi_cal2J / 1000.0 / psi_hartree2kcalmol
        #kT = Rgas * temperature
        Ezpe = 0.0
        Ethermal = (3.0/2.0) * kT
        Hthermal = Ethermal + kT
        Gthermal = 0.0              # TODO:  Sackur-Tetrode Strans+Selec. chYeah, right.

    else:
        frequency('scf')

        # thermodynamic properties
        Ethermal = psi4.get_variable('INTERNAL ENERGY CORRECTION')
        Hthermal = psi4.get_variable('ENTHALPY CORRECTION')
        Gthermal = psi4.get_variable('GIBBS FREE ENERGY CORRECTION')
        Ezpe = psi4.get_variable('ZPVE')

    Ezpe_scaled = Ezpe * ZPE_SCALE_FACTOR

    psi4.clean()

    # MP(2, FULL)/6-31G* geometry optimization
    say('MP2 optimize.')

    if nAtoms == 1:
        pass

    else:
        #psi4.set_local_option('SCF', 'SCF_TYPE', 'DIRECT')
        psi4.set_local_option('SCF', 'SCF_TYPE', 'PK')
        psi4.set_global_option('FREEZE_CORE', 'FALSE')
        psi4.set_global_option('MP2_TYPE', 'CONV')
        optimize('mp2')
        psi4.clean()

    # MP(2,frozen)/6-31G* energy
    # included to document original method
    #
    # psi4.set_global_option('FREEZE_CORE', 'TRUE')
    # psi4.set_global_option('MP2_TYPE', 'CONV')
    # optimize('mp2')
    # psi4.clean()

    if False:       # use_QCISDT is True:
        # QCISD(T,fc)/6-31G* energy for Gaussian compatibility
        say('QCISD(T).')
        psi4.set_global_option('BASIS', '6-31G*')
        psi4.set_global_option('FREEZE_CORE', 'TRUE')
        psi4.set_local_option('FNOCC', 'COMPUTE_MP4_TRIPLES', 'TRUE')
        ccref = run_fnocc('qcisd(t)', return_wfn=True, **kwargs)
        Ecc = psi4.get_variable('QCISD(T) TOTAL ENERGY')
    else:

        # CCSD(T,fc)/6-31G* energy
        say('CCSD(T).')
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DIRECT')
        psi4.set_global_option('FREEZE_CORE', 'TRUE')
        psi4.set_global_option('BASIS', '6-31G*')
        en, ccref = energy('ccsd(t)', return_wfn=True)
        Ecc = psi4.get_variable('CCSD(T) TOTAL ENERGY')

    debug('\nCC TOTAL ENERGY = %.6f Ha' % Ecc)

    Emp2_fc_631gd = psi4.get_variable('MP2 TOTAL ENERGY')
    debug('MP2 TOTAL ENERGY = %.6f Ha' % Emp2_fc_631gd)

    # ###  HLC: high-level correction based on number of valence electrons

    nfzc = nfrozen_core()
    nalpha = ccref.nalpha() - nfzc  # electron pairs
    nbeta = ccref.nbeta() - nfzc    # unpaired electrons

    debug("frozen:%d alpha:%d beta:%d" % (nfzc, nalpha, nbeta))

    # HLC constraint: nalpha >= nbeta
    if nalpha < nbeta:
        nalpha, nbeta = nbeta, nalpha

    # set QCISDT vs CCSD(T) empirical correction coefficients
    if use_QCISDT is True:
        A, B, C, D = 0.009279, 0.004471, 0.009345, 0.002021
    else:
        A, B, C, D = 0.009170, 0.004455, 0.009155, 0.001947

    # Calculate the Higher Level Correction AKA "Fudge Factor"
    if nAtoms > 1:
        Ehlc = -(A * nbeta) - B * (nalpha - nbeta)
    else:
        Ehlc = -(C * nbeta) - D * (nalpha - nbeta)

    debug("Higher Level Correction Ehlc = %.6f Ha" % (Ehlc))

    # ### Done HLC

    psi4.clean()

    # MP(2,fc)/G3MP2Large energy
    say('GMP2large.\n')

    psi4.set_global_option('BASIS', 'G3MP2Large')
    psi4.set_global_option('FREEZE_CORE', 'TRUE')
    psi4.set_local_option('SCF', 'SCF_TYPE', 'PK')
    psi4.set_global_option('MP2_TYPE', 'CONV')
    energy('mp2')

    Eg3mp2large = psi4.get_variable('MP2 TOTAL ENERGY')
    psi4.clean()

    # Approximate the MP2 correlation energy limit
    dMP2 = Eg3mp2large - Emp2_fc_631gd

    E0_g3mp2 = Ecc + dMP2 + Ezpe_scaled + Ehlc

    if nAtoms == 1:
        E0_g3mp2 += E_spin_orbit(int(mol.Z(0)), charge)

    E298_g3mp2 = E0_g3mp2 + Ethermal - Ezpe
    H298_g3mp2 = E0_g3mp2 + Hthermal - Ezpe
    G298_g3mp2 = H298_g3mp2 - (Hthermal - Gthermal)  # equiv. to "dG = dH - Tds"

    dHf0, dHf298 = calc_deltaHf(E0_g3mp2, H298_g3mp2)

    # LOG format was taken from the GAMESS-US G3MP2 implementation so that
    # data extraction tools for GAMESS G3(MP2) files may be reused here.
    #
    #    ____________________________________________________________Psi4
    #               SUMMARY OF G3(MP2, CCSD(T)) CALCULATIONS
    #    ________________________________________________________________
    #    MP2/6-31G(d)    =   -76.196848   CCSD(T)/6-31G(d) =   -76.207841
    #    MP2/G3MP2large  =   -76.314754   delta(MP2)       =    -0.117906
    #    ZPE(HF/6-31G(d))=     0.020516   ZPE Scale Factor =     0.892900
    #    HLC             =    -0.036680   Free Energy      =     0.004735
    #    Thermal Energy  =     0.025811   Thermal Enthalpy =     0.026755
    #    ________________________________________________________________
    #    E(G3(MP2)) @ 0K =   -76.341911   E(G3(MP2)) @298K =   -76.339077
    #    H(G3(MP2))      =   -76.338133   G(G3(MP2))       =   -76.360152
    #    ________________________________________________________________

    if use_QCISDT is True:
        ccstring = 'QCISDT '
    else:
        ccstring = 'CCSD(T)'

    log('    ____________________________________________________________Psi4')
    log('               SUMMARY OF G3(MP2, %s) CALCULATIONS'
        % ccstring)
    log('    ________________________________________________________________')
    log('    MP2/6-31G(d)    = % 12.6f   %s/6-31G(d) = % 12.6f'
        % (Emp2_fc_631gd, ccstring, Ecc))
    log('    MP2/G3MP2large  = % 12.6f   delta(MP2)       = % 12.6f'
        % (Eg3mp2large, dMP2))
    log('    ZPE(HF/6-31G(d))= % 12.6f   ZPE Scale Factor = % 12.6f'
        % (Ezpe_scaled, 0.8929))
    log('    HLC             = % 12.6f   Free Energy      = % 12.6f'
        % (Ehlc, Gthermal))
    log('    Thermal Energy  = % 12.6f   Thermal Enthalpy = % 12.6f'
        % (Ethermal, Hthermal))
    log('    ________________________________________________________________')
    log('    E(G3(MP2)) @ 0K = % 12.6f   E(G3(MP2)) @298K = % 12.6f'
        % (E0_g3mp2, E298_g3mp2))
    log('    H(G3(MP2))      = % 12.6f   G(G3(MP2))       = % 12.6f'
        % (H298_g3mp2, G298_g3mp2))
    log('    ________________________________________________________________')
    log("          HEAT OF FORMATION   (0K): % 10.2f kCal/mol" % dHf0)
    log("          HEAT OF FORMATION (298K): % 10.2f kCal/mol" % dHf298)
    log('    ________________________________________________________________')

    psi4.clean()
    optstash.restore()

    # return E @0K g3mp2 results

    return E0_g3mp2

# alias for g3mp2

procedures['energy']['g3mp2'] = run_g3mp2
