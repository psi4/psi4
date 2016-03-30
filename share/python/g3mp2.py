from __future__ import absolute_import
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
import sys
import re
import os
import math
import warnings
#CUimport psi4
#CUimport p4const
#CUimport p4util
from driver import *
#from extend_Molecule import *
#CUfrom molutil import *
#CUfrom p4regex import *
# never import aliases into this file


#_____________________________________________________________________________
#
# run_g3mp2 is a Psi4 implementation of Curtiss' G3(MP2) thermochemical method
# as described in:
#
#  Gaussian-3 theory using reduced Moeller-Plesset order
#  Curtiss, Redfern, Raghavachari,Rassolov, Pople
#  J.Chem.Phys., 110, 10 (1999)
#
# The CCSD(T) and HLC revisions are in accord with the later article:
#
#  Gaussian-3 theory using coupled cluster energies
#  Curtiss, Raghavachari, Redfern, Baboul, Pople
#  Chemical Physics Letters 314 (1999) 101-107
#
# Daniel R. Haney 3/28/2016 --------------------------
#
#------------------------------------------------------------------------
#   TODO:   In This Order:
#
#           handle open shell species
#           add spin orbit correction for atomic species
#           generate E0 table for row 1,2,3 atoms (less D-group elements)
#           calculate enthalpy of formation at 0K, 298K
#
#       Other:  Add QCISD(T) option for Gaussian compatibility.
#

ZPE_SCALE_FACTOR = 0.8929

#
# Zero Point Energy scale factors vary with QC applications and vibrational data sets.
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


def run_g3mp2(name, **kwargs):

    def say (s=''):
        sys.stderr.write(s)             # print to console stderr without newline

    def log (s):
        psi4.print_out(line + '\n')     # print to log file
        print(line)                     # print to console with implicit newline

    def dumpobj (obj,objname='obj'):
        for attr in dir(obj):
            print "%s.%s = %s" % (objname,attr, getattr(obj, attr))

    # spin orbit energy stub
    def E_spin_orbit (Z=0,charge=0):
        return 0.0

    # Need nAtom count from Molecule object, an implicit argument.
    kwargs = p4util.kwargs_lower(kwargs)
    mol = kwargs.pop('molecule', psi4.get_active_molecule())

    # throw an exception for open-shells
    if (psi4.get_option('SCF','REFERENCE') != 'RHF' ):
        raise ValidationError("""g3mp2 computations require "reference rhf".""")

    # stash user options:
    optstash = p4util.OptionsState(
        ['FREEZE_CORE'],
        ['MP2_TYPE'],
        ['SCF','SCF_TYPE'])

    # override default scf_type
    psi4.set_local_option('SCF','SCF_TYPE','OUT_OF_CORE')

    # optimize geometry at RHF/6-31G* for rapid frequency analysis
    say('HF optimize.')
    psi4.clean()
    psi4.set_global_option('BASIS',"6-31G(D)")
    optimize('scf')
    psi4.clean()

    # scf frequencies for zpe
    say('ZPE.')
    scf_e, ref= frequency('scf', return_wfn=True)

    # thermodynamic properties
    Ethermal = psi4.get_variable('INTERNAL ENERGY CORRECTION')
    Hthermal = psi4.get_variable('ENTHALPY CORRECTION')
    Gthermal = psi4.get_variable('GIBBS FREE ENERGY CORRECTION')

    # Calculate Zero Point Energy
    freqs   = ref.frequencies()
    nfreq   = freqs.dim(0)
    freqsum = 0.0
    for i in range (0,nfreq):
        freqsum += freqs.get(i)

    Ezpe = freqsum / p4const.psi_hartree2wavenumbers * 0.5
    Ezpe_scaled = Ezpe * ZPE_SCALE_FACTOR
    psi4.clean()

    # MP(2,FULL)/6-31G* geometry optimization
    say('MP2 optimize.')
    psi4.set_global_option('FREEZE_CORE',"FALSE")
    psi4.set_global_option('MP2_TYPE', 'CONV')
    optimize('mp2')
    psi4.clean()

    # MP(2,fc)/6-31G* energy
    say('MP2(frozen).')
    psi4.set_global_option('FREEZE_CORE',"TRUE")
    psi4.set_global_option('MP2_TYPE', 'CONV')
    energy('mp2')
    Emp2frozen = psi4.get_variable('MP2 TOTAL ENERGY')
    psi4.clean()

    # CCSD(T,fc)/6-31G* energy
    say("CCSD(T).")
    psi4.set_global_option('FREEZE_CORE',"TRUE")
    psi4.set_global_option('BASIS',"6-31G*")
    en,ccref = energy('ccsd(t)', return_wfn=True)
    Eccsdt = psi4.get_variable('CCSD(T) TOTAL ENERGY')

    ####
    # HLC: high-level correction based on number of valence electrons
    # The CC module accurately reports enough frozen core info from
    # which to infer valence electron counts.
    nirrep = ccref.nirrep()
    frzcpi = ccref.frzcpi()
    nfzc = 0
    for i in range (0,nirrep):
        nfzc += frzcpi[i]

    nalpha = ccref.nalpha() - nfzc    # electron pairs
    nbeta  = ccref.nbeta() - nfzc     # unpaired electrons
    #log("nirrep:%d frozen:%d alpha:%d beta:%d" % (nirrep,nfzc,nalpha,nbeta))

    # HLC constraint: nalpha >= nbeta
    if nalpha < nbeta:
        nalpha,nbeta = nbeta,nalpha

    # Curtiss empirical correction factors for G3(MP(2),CCSD(T))
    A = 0.009170    # molecules
    B = 0.004455    # molecules
    C = 0.009155    # atoms
    D = 0.001947    # atoms

    if mol.natom() > 1:
        Ehlc = -(A * nbeta) - B * (nalpha - nbeta)
    else:
        Ehlc = -(C * nbeta) - D * (nalpha - nbeta)
    ####

    psi4.clean()

    # MP(2,fc)/G3MP2Large energy
    say('GMP2large.\n')
    psi4.set_global_option('BASIS',"G3MP2Large")
    psi4.set_global_option('FREEZE_CORE',"TRUE")
    psi4.set_global_option('MP2_TYPE', 'CONV')
    energy('mp2')
    Eg3mp2large = psi4.get_variable('MP2 TOTAL ENERGY')
    psi4.clean()

    dMP2 = Eg3mp2large - Emp2frozen

    E0   = Eccsdt + dMP2 + Ezpe_scaled + Ehlc
    E298 = E0 + (Ethermal - Ezpe)
    H298 = E0 + (Hthermal - Ezpe)
    G298 = H298 - (Hthermal - Gthermal)     # equiv. to "dG = dH - Tds"


# This format was taken from the GAMESS-US G3MP2 implementation.
#
#    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Psi4
#               SUMMARY OF G3(MP2,CCSD(T)) CALCULATIONS
#    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    MP2/6-31G(d)    =   -76.196848   CCSD(T)/6-31G(d) =   -76.207841
#    MP2/G3MP2large  =   -76.314754   delta(MP2)       =    -0.117906
#    ZPE(HF/6-31G(d))=     0.020516   ZPE Scale Factor =     0.892900
#    HLC             =    -0.036680   Free Energy      =     0.004735
#    Thermal Energy  =     0.025811   Thermal Enthalpy =     0.026755
#    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    E(G3(MP2)) @ 0K =   -76.341911   E(G3(MP2)) @298K =   -76.339077
#    H(G3(MP2))      =   -76.338133   G(G3(MP2))       =   -76.360152
#    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    summary = [
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Psi4"),
        ("               SUMMARY OF G3(MP2,CCSD(T)) CALCULATIONS              "),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        ("    MP2/6-31G(d)    = % 12.6f   CCSD(T)/6-31G(d) = % 12.6f" %
         (Emp2frozen, Eccsdt)),
        ("    MP2/G3MP2large  = % 12.6f   delta(MP2)       = % 12.6f" %
         (Eg3mp2large, dMP2)),
        ("    ZPE(HF/6-31G(d))= % 12.6f   ZPE Scale Factor = % 12.6f" %
         (Ezpe_scaled, 0.8929)),
        ("    HLC             = % 12.6f   Free Energy      = % 12.6f" %
         (Ehlc, Gthermal)),
        ("    Thermal Energy  = % 12.6f   Thermal Enthalpy = % 12.6f" %
         (Ethermal, Hthermal)),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"),
        ("    E(G3(MP2)) @ 0K = % 12.6f   E(G3(MP2)) @298K = % 12.6f" %
         (E0, E298)),
        ("    H(G3(MP2))      = % 12.6f   G(G3(MP2))       = % 12.6f" %
         (H298, G298)),
        ("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        ]

    for line in summary:
        log(line)                     # print to console stderr

    psi4.clean()
    optstash.restore()

    # return E @0K g3mp2 results
    return E0

# alias for g3mp2
procedures['energy']['g3mp2'] = run_g3mp2


