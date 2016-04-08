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
import os
import math
import warnings

from driver import *


# _____________________________________________________________________________
#
# run_g3mp2 is a Psi4 implementation of Curtiss' G3(MP2) thermochemical method.
#
# invoke:   en = energy('g3mp2')
#
# returns:  float en = g3mp2 total energy E0 in Hartrees,
#
#
# It is described in:
#
#  Gaussian-3 theory using reduced Moeller-Plesset order
#  Curtiss, Redfern, Raghavachari, Rassolov, Pople
#  J.Chem.Phys., 110, 10(1999)
#
# The CCSD(T) and HLC revisions are included and described in:
#
#  Gaussian-3 theory using coupled cluster energies
#  Curtiss, Raghavachari, Redfern, Baboul, Pople
#  Chem.Phys.Lett., 314(1999) 101-107
#
# Daniel R. Haney 3/28/2016 --------------------------
#
# ------------------------------------------------------------------------
#   TODO:
#       Add memory check on MP(2, fc)/G3MP2Large calc
#       to select for OUT_OF_CORE when it may exhaust memory.
#

# ________________________________________________________________________
# ________________________________________________________________________

# global

def run_g3mp2(name, **kwargs):

    g3mp2_verbose_flag = False
    use_QCISDT = False

    # Boltzmann kT contribution @ 298.15K in Hartrees
    psi_kT = (psi_kb * 298.15 * psi_na) / psi_cal2J / psi_hartree2kcalmol / 1000.0

    # Curtiss' HF/6-31G* Zero Point Energy scale factor from Gaussian-94.
    # Note that it varies by QC application and data set.

    ZPE_SCALE_FACTOR = 0.8929

    # ________________________________________________________________________
    # write to console stderr
    def say(s=''):
        if g3mp2_verbose_flag is True:
            sys.stderr.write(s)

    # write to log file
    def log(line):
        psi4.print_out(line + '\n')

    # write to log file and console
    def report(line):
        log(line)
        if g3mp2_verbose_flag is True:
            print line

    # ________________________________________________________________________
    def g3mp2_usage():

        log("\n")
        log("g3mp2 flags are:")
        log(" qcisdt  - calculate QCISD(T) energy. default: CCSD(T)")
        log(" verbose - write progress to console")
        log(" help    - this")
        log("\n")

    # ________________________________________________________________________
    def nfrozen_core(molecule):
        '''
        Calculate frozen core electron pair count for a molecule

        returns: frozen core count, integer
        '''

        _nfcz = 0
        natoms = molecule.natom()
        for A in range(0, natoms):
            Z = int(molecule.Z(A))
            if Z > 2:
                _nfcz += 1
            if Z > 10:
                _nfcz += 4
            if Z > 18:
                _nfcz += 4
            if Z > 36:
                _nfcz += 9
            if Z > 54:
                _nfcz += 9
            if Z > 86:
                _nfcz += 16
            if Z > 108:
                raise ValidationError("nfrozen_core: Invalid atomic number %d" % Z)
        return _nfcz

    # ________________________________________________________________________

    def empirical_correction(molecule, alpha, beta, qcisdt_f=False):
        '''
        Higher Level Correction of E0
        Applies Empirical coefficients derived from
        fit between experimental dHf values and
        valence electron counts.

        returns:  HLC in Hartrees,  float
        '''
        nAtoms = molecule.natom()
        nfzc = nfrozen_core(molecule)

        nalpha = alpha - nfzc  # electron pairs
        nbeta = beta - nfzc    # unpaired electrons

        # HLC constraint: nalpha >= nbeta
        if nalpha < nbeta:
            nalpha, nbeta = nbeta, nalpha

        # set QCISDT vs CCSD(T) empirical correction coefficients
        if qcisdt_f is True:
            A, B, C, D = 0.009279, 0.004471, 0.009345, 0.002021
        else:
            A, B, C, D = 0.009170, 0.004455, 0.009155, 0.001947

        if nAtoms > 1:
            hlc = - (A * nbeta) - B * (nalpha - nbeta)
        else:
            hlc = - (C * nbeta) - D * (nalpha - nbeta)

        return hlc

    # ________________________________________________________________________

    def spin_orbit_correction(Z=0, charge=0):
        '''
        spin orbit energy corrections for atomic species
        with respect to charge. See Curtiss ref. above.

        returns: value in Hartrees, float
        '''
        #       [neutral, Cation+, Anion- ]
        _ESO = [[0.0, 0.0, 0.00],       # 00 placeholder
                [0.0, 0.0, 0.00],       # 01  Hydrogen
                [0.0, 0.0, 0.00],       # 02  Helium
                [0.0, 0.0, 0.00],       # 03  Lithium
                [0.0, 0.0, 0.00],       # 04  Beryllium
                [-0.05, 0.0, -0.03],    # 05  Boron
                [-0.14, -0.2, 0.00],    # 06  Carbon
                [0.0, -0.43, 0.00],     # 07  Nitrogen
                [-0.36, 0.0, -0.26],    # 08  Oxygen
                [-0.61, -0.67, 0.00],   # 09  Fluorine
                [0.0, -1.19, 0.00],     # 10  Neon
                [0.0, 0.0, 0.00],       # 11  Sodium
                [0.0, 0.0, 0.00],       # 12  Magnesium
                [-0.34, 0.0, -0.28],    # 13  Aluminum
                [-0.68, -0.93, 0.00],   # 14  Silicon
                [0.0, -1.43, -0.45],    # 15  Phosphorus
                [-0.89, 0.0, -0.88],    # 16  Sulfur
                [-1.34, -1.68, 0.00],   # 17  Chlorine
                [0.0, -2.18, 0.00],     # 18  Argon
                [0.0, 0.0, 0.00],       # 19  Potassium
                [0.0, 0.0, 0.00],       # 20  Calcium
                [0.0, 0.0, 0.00],       # 21  Scandium
                [0.0, 0.0, 0.00],       # 22  Titanium
                [0.0, 0.0, 0.00],       # 23  Vanadium
                [0.0, 0.0, 0.00],       # 24  Chromium
                [0.0, 0.0, 0.00],       # 25  Manganese
                [0.0, 0.0, 0.00],       # 26  Iron
                [0.0, 0.0, 0.00],       # 27  Cobalt
                [0.0, 0.0, 0.00],       # 28  Nickel
                [0.0, 0.0, 0.00],       # 29  Copper
                [0.0, 0.0, 0.00],       # 30  Zinc
                [-2.51, 0.0, 0.00],     # 31  Gallium
                [-4.41, -5.37, 0.00],   # 32  Germanium
                [0.0, -8.04, 0.00],     # 33  Arsenic
                [-4.3, 0.0, 0.00],      # 34  Selenium
                [-5.6, -6.71, 0.00],    # 35  Bromine
                [0.0, -8.16, 0.00]]     # 36  Krypton

        if Z not in range(len(_ESO)):
            raise ValidationError("spin_orbit_correction: Invalid atomic number %d" % Z)

        if charge > 0:
            ion = 1
        elif charge < 0:
            ion = 2
        else:
            ion = 0

        milliHa_to_Ha = 0.001
        return(_ESO[Z][ion] * milliHa_to_Ha)

    # ________________________________________________________________________

    def E0_ccsdt(Z=0):
        '''
        G3(MP2, CCSD(T)) @0K energies of atomic species.

        returns:  atom E0 in Hartrees, float
        '''

        _ccsdtE0 = \
            [0.0,       # 0 index placeholder
             -0.501765, -2.902353, -7.366280, -14.629072,          # H  He Li Be
             -24.606789, -37.788989, -54.524780, -74.989201,       # B  C  N  O
             -99.640199, -128.827752, -161.847930, -199.650655,    # F  Ne Na Mg
             -241.936660, -288.939067, -340.826225, -397.663216,   # Al Si P  S
             -459.686583, -527.060194, -599.160438, -676.790288,   # Cl Ar K  Ca
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,     # Sc-Zn
             -1923.582394, -2075.683666, -2234.561806,             # Ga  Ge  As
             -2400.225135, -2572.828929, -2752.462850]             # Se  Br  Kr

        if Z < 1 and Z > len(_ccsdtE0):
            raise ValidationError("E0_ccsdt: Invalid atomic number %d" % Z)
        else:
            return _ccsdtE0[Z]
    # ________________________________________________________________________

    def E0_qcisdt(Z=0):
        '''
        G3(MP2, QCISDT) @0K energies of atomic species.
          Values from Nwchem until UHF-QCISD(T) is implemented.

        returns:  atom E0 in Hartrees
        '''

        _qcisdtE0 = \
            [0.0,
             -0.501839, -2.902543, -7.434048, -14.629262,           # H  He Li Be
             -24.607093, -37.789349, -54.525198, -74.989850,        # B  C  N  O
             -99.641120, -128.828970, -161.848004, -199.650845,     # F  Ne Na Mg
             -241.936973, -288.939460, -340.826670, -397.663794,    # Al Si P  S
             -459.687272, -527.060963, -599.160512, -676.789424,    # Cl Ar K  Ca
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,      # Sc-Zn
             -1923.585079, -2075.685727, -2234.563599,              # Ga  Ge  As
             -2400.226922, -2572.830747, -2752.464704]              # Se  Br  Kr

        if Z < 1 and Z > len(_qcisdtE0):
            raise ValidationError("E0_qcisdt: Invalid atomic number %d" % Z)

        else:
            return _qcisdtE0[Z]
    # ________________________________________________________________________

    def E0_one_atom(Z=0, qcisdt_f=False):

        if qcisdt_f is True:
            _e0 = E0_qcisdt(Z)
        else:
            _e0 = E0_ccsdt(Z)

        _h298 = _e0 + (5.0 / 2.0) * psi_kT  # Ideal gas kinetic energy contribution

        return(_e0, _h298)
    # ________________________________________________________________________

    def dHf_one_atom(Z=0):
        '''
        return experimental atomic heats of formation.

        Values harvested from the NIST Webbook at webbook.nist.gov/chemistry

        returns:    dHf(0K), dHf(298K)  tuple of 2 floats
        '''

        #    atom  [dHf(0), dHf(298)]  in kcal/mol
        atomDHF = [[0.0, 0.0],         # 00 place holder
                   [51.63, 52.103],    # 01 H     Hydrogen
                   [0.00, 0.00],       # 02 He    Helium
                   [37.70, 38.07],     # 03 Li    Lithium
                   [76.40, 77.40],     # 04 Be    Beryllium
                   [135.10, 136.30],   # 05 B     Boron
                   [169.98, 171.29],   # 06 C     Carbon
                   [112.53, 112.97],   # 07 N     Nitrogen
                   [58.99, 59.56],     # 08 O     Oxygen
                   [18.47, 18.97],     # 09 F     Fluorine
                   [0.00, 0.00],       # 10 Ne    Neon
                   [25.76, 25.69],     # 11 Na    Sodium
                   [34.87, 35.16],     # 12 Mg    Magnesium
                   [80.20, 80.80],     # 13 Al    Aluminum
                   [107.20, 108.20],   # 14 Si    Silicon
                   [75.45, 75.65],     # 15 P     Phosphorus
                   [65.71, 66.25],     # 16 S     Sulfur
                   [28.59, 28.99],     # 17 Cl    Chlorine
                   [0.00, 0.00],       # 18 Ar    Argon
                   [21.27, 21.49],     # 19 K     Potassium
                   [42.50, 42.29],     # 20 Ca    Calcium
                   [0.0, 0.0],         # 21 Sc    Scandium
                   [0.0, 0.0],         # 22 Ti    Titanium
                   [0.0, 0.0],         # 23 V     Vanadium
                   [0.0, 0.0],         # 24 Cr    Chromium
                   [0.0, 0.0],         # 25 Mn    Manganese
                   [0.0, 0.0],         # 26 Fe    Iron
                   [0.0, 0.0],         # 27 Co    Cobalt
                   [0.0, 0.0],         # 28 Ni    Nickel
                   [0.0, 0.0],         # 29 Cu    Copper
                   [0.0, 0.0],         # 30 Zn    Zinc
                   [65.00, 65.00],     # 31 Ga    Gallium
                   [88.91, 88.91],     # 32 Ge    Germanium
                   [73.90, 72.42],     # 33 As    Arsenic
                   [55.76, 54.27],     # 34 Se    Selenium
                   [26.74, 28.18],     # 35 Br    Bromine
                   [0.0, 0.0]]         # 36 Kr    Krypton

        if Z < 1 or Z >= len(atomDHF):
            raise ValidationError("dHf_one_atom: Invalid atomic number %d" % Z)
        else:
            return (atomDHF[Z][0], atomDHF[Z][1])

    # ________________________________________________________________________

    def heat_of_formation(molecule, _E0, _H298, qcisdt_f=False):
        '''
        calculate heats of formation at
        0K and 298K by atomization method

        input:    ab initio E0, H298 in Hartrees
        returns: dHf@0K, dHf@298K -> 2 doubles
        '''

        E0_sum_atoms = 0.0
        H298_sum_atoms = 0.0
        dhf0_sum_atoms = 0.0
        dhf298_sum_atoms = 0.0

        _natoms = molecule.natom()

        for A in range(0, _natoms):
            Z = int(molecule.Z(A))

            e0_atom, h298_atom = E0_one_atom(Z, qcisdt_f)     # in Hartrees
            E0_sum_atoms += e0_atom
            H298_sum_atoms += h298_atom

            d0, d298 = dHf_one_atom(Z)              # in kcal/mol
            dhf0_sum_atoms += d0
            dhf298_sum_atoms += d298

        dhf0 = (_E0 - E0_sum_atoms) * psi_hartree2kcalmol + dhf0_sum_atoms
        dhf298 = (_H298 - H298_sum_atoms) * psi_hartree2kcalmol + dhf298_sum_atoms

        # Eatomization = -(_H298 - H298_sum_atoms) * psi_hartree2kcalmol

        return dhf0, dhf298

    # ______________________________________________________
    # _____________________      ___________________________
    # _____________________ MAIN ___________________________
    # ______________________________________________________

    kwargs = p4util.kwargs_lower(kwargs)

    if "help" in kwargs:
        g3mp2_usage()
        return 0.0      # quit

    g3mp2_verbose_flag = kwargs.pop('verbose', False)
    use_QCISDT = kwargs.pop('qcisdt', False)

    g3mp2_mol = kwargs.pop('molecule')

    g3mp2_name = g3mp2_mol.name()
    nAtoms = g3mp2_mol.natom()
    charge = g3mp2_mol.molecular_charge()
    multiplicity = g3mp2_mol.multiplicity()

    if multiplicity > 1:
        if use_QCISDT is True:      # UHF-QCISD(T) unsupported
            report("\n")
            report("run_g3mp2: QCISD(T) doesn't support species with multiplicity > 1.")
            report("run_g3mp2: Use the default CCSD(T) if acceptable.\n")
            psi4.clean()
            raise ValidationError("UHF unsupported for QCISD(T)\n")

        else:
            psi4.set_global_option('REFERENCE', 'UHF')

    else:
        psi4.set_global_option('REFERENCE', 'RHF')

    # stash user options:
    optstash = p4util.OptionsState(
        ['FNOCC', 'COMPUTE_TRIPLES'],
        ['FNOCC', 'COMPUTE_MP4_TRIPLES'],
        ['FREEZE_CORE'],
        ['MP2_TYPE'],
        ['SCF', 'SCF_TYPE'])

    # frobnicate the scf_type if necessary
    # DF diverges from PK | DIRECT beyond 1.0E-04 Ha ~= 1/15 kcal
    # psi4.set_local_option('SCF', 'SCF_TYPE', 'DIRECT')
    psi4.set_local_option('SCF', 'SCF_TYPE', 'PK')

    # optimize geometry at RHF/6-31G*
    say('HF optimize. ')

    if nAtoms > 1:
        psi4.set_global_option('BASIS', '6-31G(D)')
        opt_e = optimize('scf', molecule=g3mp2_mol)
        g3mp2_mol = psi4.get_active_molecule()
        psi4.clean()

    else:
        pass

    # Harmonic Zero Point Energy, etc.
    say('ZPE. ')

    # until frequency() handles atoms properly.
    if nAtoms == 1:
        Ezpe = 0.0
        Ethermal = (3.0 / 2.0) * psi_kT     # 3/2 kT
        Hthermal = Ethermal + psi_kT
        Gthermal = 0.0              # TODO: Sackur-Tetrode Strans+Selec.

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
    say('MP2 optimize. ')

    if nAtoms == 1:
        pass

    else:

        psi4.set_local_option('SCF', 'SCF_TYPE', 'PK')
        psi4.set_global_option('FREEZE_CORE', 'FALSE')
        psi4.set_global_option('MP2_TYPE', 'CONV')
        opt_e, opt_wfn = optimize('mp2', return_wfn=True, molecule=g3mp2_mol)
        g3mp2_mol = opt_wfn.molecule()
        psi4.clean()

    # included to document original method
    #
    # MP(2, frozen)/6-31G* energy
    # psi4.set_local_option('SCF', 'SCF_TYPE', 'PK')
    # psi4.set_global_option('FREEZE_CORE', 'TRUE')
    # psi4.set_global_option('MP2_TYPE', 'CONV')
    # Emp2 = energy('mp2')
    # psi4.clean()

    if use_QCISDT is True:
        # QCISD(T, fc)/6-31G* energy for Gaussian compatibility
        say('QCISD(T). ')

        psi4.set_local_option('FNOCC', 'COMPUTE_MP4_TRIPLES', 'TRUE')
        # psi4.set_local_option('SCF', 'SCF_TYPE', 'DIRECT')
        psi4.set_global_option('FREEZE_CORE', 'TRUE')
        psi4.set_global_option('BASIS', '6-31G*')
        ccref = run_fnocc('qcisd(t)', return_wfn=True, molecule=g3mp2_mol)
        g3mp2_mol = ccref.molecule()
        Ecc = psi4.get_variable('QCISD(T) TOTAL ENERGY')

    else:
        # CCSD(T, fc)/6-31G* energy
        say('CCSD(T). ')
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DIRECT')
        psi4.set_global_option('FREEZE_CORE', 'TRUE')
        psi4.set_global_option('BASIS', '6-31G*')
        en, ccref = energy('ccsd(t)', return_wfn=True, molecule=g3mp2_mol)
        g3mp2_mol = ccref.molecule()

        Ecc = psi4.get_variable('CCSD(T) TOTAL ENERGY')

    Emp2_fc_631gd = psi4.get_variable('MP2 TOTAL ENERGY')

    Ehlc = empirical_correction(g3mp2_mol, ccref.nalpha(), ccref.nbeta(), use_QCISDT)

    psi4.clean()

    # MP(2, fc)/G3MP2Large energy
    say('GMP2large.\n')

    psi4.set_global_option('BASIS', 'G3MP2Large')
    psi4.set_global_option('FREEZE_CORE', 'TRUE')
    psi4.set_local_option('SCF', 'SCF_TYPE', 'PK')
    psi4.set_global_option('MP2_TYPE', 'CONV')
    en, wfn_g3mp2 = energy('mp2', molecule=g3mp2_mol, return_wfn=True)

    Eg3mp2large = psi4.get_variable('MP2 TOTAL ENERGY')
    psi4.clean()

    # Approximate the MP2 correlation energy limit
    dMP2 = Eg3mp2large - Emp2_fc_631gd

    E0_g3mp2 = Ecc + dMP2 + Ezpe_scaled + Ehlc

    if nAtoms == 1:
        Eso = spin_orbit_correction(int(g3mp2_mol.Z(0)), charge)
        E0_g3mp2 += Eso

    E298_g3mp2 = E0_g3mp2 + Ethermal - Ezpe
    H298_g3mp2 = E0_g3mp2 + Hthermal - Ezpe
    G298_g3mp2 = H298_g3mp2 - (Hthermal - Gthermal)
    psi4.set_variable('CURRENT ENERGY', E0_g3mp2)

    dHf0, dHf298 = heat_of_formation(g3mp2_mol, E0_g3mp2, H298_g3mp2, use_QCISDT)

    '''
    Summary format was taken from GAMESS-US G3MP2 implementation so that
    data extraction tools for GAMESS G3(MP2) files may be reused.

    ____________________________________________________________Psi4
               SUMMARY OF G3(MP2, CCSD(T)) CALCULATIONS
    ________________________________________________________________
    MP2/6-31G(d)    =   -76.196848   CCSD(T)/6-31G(d) =   -76.207840
    MP2/G3MP2large  =   -76.314759   delta(MP2)       =    -0.117911
    ZPE(HF/6-31G(d))=     0.020517   ZPE Scale Factor =     0.892900
    HLC             =    -0.036680   Free Energy      =     0.004726
    Thermal Energy  =     0.025811   Thermal Enthalpy =     0.026756
    ________________________________________________________________
    E(G3(MP2)) @ 0K =   -76.341914   E(G3(MP2)) @298K =   -76.339080
    H(G3(MP2))      =   -76.338136   G(G3(MP2))       =   -76.360165
    ________________________________________________________________
          HEAT OF FORMATION  (0K):     -56.87 kCal/mol
          HEAT OF FORMATION(298K):     -57.42 kCal/mol
    ________________________________________________________________
    '''

    if use_QCISDT is True:
        ccstring = 'QCISDT '
    else:
        ccstring = 'CCSD(T)'

    report('    ____________________________________________________________Psi4')
    report('               SUMMARY OF G3(MP2, %s) CALCULATIONS'
           % ccstring)
    report('    ________________________________________________________________')
    report('    MP2/6-31G(d)    = %12.6f   %s/6-31G(d) = %12.6f'
           % (Emp2_fc_631gd, ccstring, Ecc))
    report('    MP2/G3MP2large  = %12.6f   delta(MP2)       = %12.6f'
           % (Eg3mp2large, dMP2))
    report('    ZPE(HF/6-31G(d))= %12.6f   ZPE Scale Factor = %12.6f'
           % (Ezpe_scaled, 0.8929))
    report('    HLC             = %12.6f   Free Energy      = %12.6f'
           % (Ehlc, Gthermal))
    report('    Thermal Energy  = %12.6f   Thermal Enthalpy = %12.6f'
           % (Ethermal, Hthermal))
    report('    ________________________________________________________________')
    report('    E(G3(MP2)) @ 0K = %12.6f   E(G3(MP2)) @298K = %12.6f'
           % (E0_g3mp2, E298_g3mp2))
    report('    H(G3(MP2))      = %12.6f   G(G3(MP2))       = %12.6f'
           % (H298_g3mp2, G298_g3mp2))

    report('    ________________________________________________________________')
    report("          HEAT OF FORMATION  (0K): % 10.2f kCal/mol" % dHf0)
    report("          HEAT OF FORMATION(298K): % 10.2f kCal/mol" % dHf298)
    report('    ________________________________________________________________')

    # write final energies in CSV format for trivial extraction.
    # A database of several molecules may be roughly constructed by common
    # shell tools, e.g.,
    #
    #  grep -h -A3 "CSV format" *.out |grep -v CSV |sort -r | uniq >moldb.csv

    log("\n")
    log("## CSV format of G3(MP2) energies")

    log("#specie, E0, H298, dHf0, dHf298")
    log("#unit, Ha, Ha, kcal, kcal")
    log("\"%s\",%.6f,%.6f,%.2f,%.2f" %
        (g3mp2_name,
         E0_g3mp2, H298_g3mp2,
         dHf0, dHf298))

    log("\n")

    psi4.clean()
    optstash.restore()

    # return current energy E @0K g3mp2 and arbitrary wave function

    return E0_g3mp2, wfn_g3mp2

# alias for g3mp2

procedures['energy']['g3mp2'] = run_g3mp2
