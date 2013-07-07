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

# Gn theory.

import re
import os
import math
import warnings
import psi4
import p4const
import p4util
from driver import *
#from extend_Molecule import *
from molutil import *
from p4regex import *
# never import aliases into this file

def run_gaussian_2(name, **kwargs):

    # throw an exception for open-shells
    if (psi4.get_option('SCF','REFERENCE') != 'RHF' ):
        raise ValidationError("""g2 computations require "reference rhf".""")

    # stash user options:
    optstash = p4util.OptionsState(
        ['FNOCC','COMPUTE_TRIPLES'],
        ['FNOCC','COMPUTE_MP4_TRIPLES'],
        ['FREEZE_CORE'],
        ['SCF','SCF_TYPE'])

    # override default scf_type
    psi4.set_local_option('SCF','SCF_TYPE','OUT_OF_CORE')

    # optimize geometry at scf level
    psi4.clean()
    psi4.set_global_option('BASIS',"6-31G(D)")
    optimize('scf')
    psi4.clean()

    # scf frequencies for zpe
    frequency('scf')

    # thermodynamic properties
    du = psi4.get_variable('INTERNAL ENERGY CORRECTION')
    dh = psi4.get_variable('ENTHALPY CORRECTION')
    dg = psi4.get_variable('GIBBS FREE ENERGY CORRECTION')

    ref     = psi4.wavefunction()
    freqs   = ref.frequencies()
    nfreq   = freqs.dim(0)
    freqsum = 0.0
    for i in range (0,nfreq):
        freqsum += freqs.get(i)
    zpe = freqsum / p4const.psi_hartree2wavenumbers * 0.8929 * 0.5
    psi4.clean()

    # optimize geometry at mp2 (no frozen core) level
    # note: freeze_core isn't an option in MP2
    psi4.set_global_option('FREEZE_CORE',"FALSE")
    optimize('conv-mp2')
    psi4.clean()

    # qcisd(t)
    psi4.set_local_option('FNOCC','COMPUTE_MP4_TRIPLES',"TRUE")
    psi4.set_global_option('FREEZE_CORE',"TRUE")
    psi4.set_global_option('BASIS',"6-311G(D_P)")
    run_fnocc('qcisd(t)',**kwargs)

    # HLC: high-level correction based on number of valence electrons
    ref    = psi4.wavefunction()
    nirrep = ref.nirrep()
    frzcpi = ref.frzcpi()
    nfzc = 0
    for i in range (0,nirrep):
        nfzc += frzcpi[i]
    nalpha = ref.nalpha() - nfzc
    nbeta  = ref.nbeta() - nfzc
    # hlc of gaussian-2
    hlc = -0.00481 * nalpha -0.00019 * nbeta
    # hlc of gaussian-1
    hlc1 = -0.00614 * nalpha

    eqci_6311gdp = psi4.get_variable("QCISD(T) TOTAL ENERGY")
    emp4_6311gd  = psi4.get_variable("MP4 TOTAL ENERGY")
    emp2_6311gd  = psi4.get_variable("MP2 TOTAL ENERGY")
    psi4.clean()

    # correction for diffuse functions
    psi4.set_global_option('BASIS',"6-311+G(D_P)")
    energy('mp4')
    emp4_6311pg_dp = psi4.get_variable("MP4 TOTAL ENERGY")
    emp2_6311pg_dp = psi4.get_variable("MP2 TOTAL ENERGY")
    psi4.clean()

    # correction for polarization functions
    psi4.set_global_option('BASIS',"6-311G(2DF_P)")
    energy('mp4')
    emp4_6311g2dfp = psi4.get_variable("MP4 TOTAL ENERGY")
    emp2_6311g2dfp = psi4.get_variable("MP2 TOTAL ENERGY")
    psi4.clean()

    # big basis mp2
    psi4.set_global_option('BASIS',"6-311+G(3DF_2P)")
    run_fnocc('_mp2',**kwargs)
    emp2_big = psi4.get_variable("MP2 TOTAL ENERGY")
    psi4.clean()

    eqci       = eqci_6311gdp
    e_delta_g2 = emp2_big + emp2_6311gd - emp2_6311g2dfp - emp2_6311pg_dp
    e_plus     = emp4_6311pg_dp - emp4_6311gd
    e_2df      = emp4_6311g2dfp - emp4_6311gd

    eg2 = eqci + e_delta_g2 + e_plus + e_2df
    eg2_mp2_0k = eqci + (emp2_big - emp2_6311gd) + hlc + zpe

    psi4.print_out('\n')
    psi4.print_out('  ==>  G1/G2 Energy Components  <==\n')
    psi4.print_out('\n')
    psi4.print_out('        QCISD(T):            %20.12lf\n' % eqci)
    psi4.print_out('        E(Delta):            %20.12lf\n' % e_delta_g2)
    psi4.print_out('        E(2DF):              %20.12lf\n' % e_2df)
    psi4.print_out('        E(+):                %20.12lf\n' % e_plus)
    psi4.print_out('        E(G1 HLC):           %20.12lf\n' % hlc1)
    psi4.print_out('        E(G2 HLC):           %20.12lf\n' % hlc)
    psi4.print_out('        E(ZPE):              %20.12lf\n' % zpe)
    psi4.print_out('\n')
    psi4.print_out('  ==>  0 Kelvin Results  <==\n')
    psi4.print_out('\n')
    eg2_0k = eg2 + zpe + hlc
    psi4.print_out('        G1:                  %20.12lf\n' % (eqci + e_plus + e_2df + hlc1 + zpe))
    psi4.print_out('        G2(MP2):             %20.12lf\n' % eg2_mp2_0k)
    psi4.print_out('        G2:                  %20.12lf\n' % eg2_0k)

    psi4.set_variable("G1 TOTAL ENERGY",eqci + e_plus + e_2df + hlc1 + zpe)
    psi4.set_variable("G2 TOTAL ENERGY",eg2_0k)
    psi4.set_variable("G2(MP2) TOTAL ENERGY",eg2_mp2_0k)

    psi4.print_out('\n')
    T = psi4.get_global_option('T')
    psi4.print_out('  ==>  %3.0lf Kelvin Results  <==\n'% T)
    psi4.print_out('\n')

    internal_energy = eg2_mp2_0k + du - zpe / 0.8929
    enthalpy        = eg2_mp2_0k + dh - zpe / 0.8929
    gibbs           = eg2_mp2_0k + dg - zpe / 0.8929

    psi4.print_out('        G2(MP2) energy:      %20.12lf\n' % internal_energy )
    psi4.print_out('        G2(MP2) enthalpy:    %20.12lf\n' % enthalpy)
    psi4.print_out('        G2(MP2) free energy: %20.12lf\n' % gibbs)
    psi4.print_out('\n')

    psi4.set_variable("G2(MP2) INTERNAL ENERGY",internal_energy)
    psi4.set_variable("G2(MP2) ENTHALPY",enthalpy)
    psi4.set_variable("G2(MP2) FREE ENERGY",gibbs)

    internal_energy = eg2_0k + du - zpe / 0.8929
    enthalpy        = eg2_0k + dh - zpe / 0.8929
    gibbs           = eg2_0k + dg - zpe / 0.8929

    psi4.print_out('        G2 energy:           %20.12lf\n' % internal_energy )
    psi4.print_out('        G2 enthalpy:         %20.12lf\n' % enthalpy)
    psi4.print_out('        G2 free energy:      %20.12lf\n' % gibbs)

    psi4.set_variable("CURRENT ENERGY",eg2_0k)

    psi4.set_variable("G2 INTERNAL ENERGY",internal_energy)
    psi4.set_variable("G2 ENTHALPY",enthalpy)
    psi4.set_variable("G2 FREE ENERGY",gibbs)

    psi4.clean()

    optstash.restore()

    # return 0K g2 results
    return eg2_0k

# aliases for g2
procedures['energy']['gaussian-2'] = run_gaussian_2
procedures['energy']['g2']         = run_gaussian_2
