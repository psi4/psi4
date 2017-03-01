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

# Gn theory.

import re
import os
import math
import warnings
from psi4 import core
from psi4.driver import driver
from psi4.driver import p4util
from psi4.driver import constants
#from driver import *
# never import aliases into this file

def run_gaussian_2(name, **kwargs):

    # throw an exception for open-shells
    if (core.get_option('SCF','REFERENCE') != 'RHF' ):
        raise ValidationError("""g2 computations require "reference rhf".""")

    # stash user options:
    optstash = p4util.OptionsState(
        ['FNOCC','COMPUTE_TRIPLES'],
        ['FNOCC','COMPUTE_MP4_TRIPLES'],
        ['FREEZE_CORE'],
        ['MP2_TYPE'],
        ['SCF','SCF_TYPE'])

    # override default scf_type
    core.set_local_option('SCF','SCF_TYPE','PK')

    # optimize geometry at scf level
    core.clean()
    core.set_global_option('BASIS',"6-31G(D)")
    driver.optimize('scf')
    core.clean()

    # scf frequencies for zpe
    # NOTE This line should not be needed, but without it there's a seg fault
    scf_e, ref = driver.frequency('scf', return_wfn=True)

    # thermodynamic properties
    du = core.get_variable('INTERNAL ENERGY CORRECTION')
    dh = core.get_variable('ENTHALPY CORRECTION')
    dg = core.get_variable('GIBBS FREE ENERGY CORRECTION')

    freqs   = ref.frequencies()
    nfreq   = freqs.dim(0)
    freqsum = 0.0
    for i in range(0, nfreq):
        freqsum += freqs.get(i)
    zpe = freqsum / constants.hartree2wavenumbers * 0.8929 * 0.5
    core.clean()

    # optimize geometry at mp2 (no frozen core) level
    # note: freeze_core isn't an option in MP2
    core.set_global_option('FREEZE_CORE',"FALSE")
    core.set_global_option('MP2_TYPE', 'CONV')
    driver.optimize('mp2')
    core.clean()

    # qcisd(t)
    core.set_local_option('FNOCC','COMPUTE_MP4_TRIPLES',"TRUE")
    core.set_global_option('FREEZE_CORE',"TRUE")
    core.set_global_option('BASIS',"6-311G(D_P)")
    ref = driver.proc.run_fnocc('qcisd(t)', return_wfn=True, **kwargs)

    # HLC: high-level correction based on number of valence electrons
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

    eqci_6311gdp = core.get_variable("QCISD(T) TOTAL ENERGY")
    emp4_6311gd  = core.get_variable("MP4 TOTAL ENERGY")
    emp2_6311gd  = core.get_variable("MP2 TOTAL ENERGY")
    core.clean()

    # correction for diffuse functions
    core.set_global_option('BASIS',"6-311+G(D_P)")
    driver.energy('mp4')
    emp4_6311pg_dp = core.get_variable("MP4 TOTAL ENERGY")
    emp2_6311pg_dp = core.get_variable("MP2 TOTAL ENERGY")
    core.clean()

    # correction for polarization functions
    core.set_global_option('BASIS',"6-311G(2DF_P)")
    driver.energy('mp4')
    emp4_6311g2dfp = core.get_variable("MP4 TOTAL ENERGY")
    emp2_6311g2dfp = core.get_variable("MP2 TOTAL ENERGY")
    core.clean()

    # big basis mp2
    core.set_global_option('BASIS',"6-311+G(3DF_2P)")
    #run_fnocc('_mp2',**kwargs)
    driver.energy('mp2')
    emp2_big = core.get_variable("MP2 TOTAL ENERGY")
    core.clean()
    eqci       = eqci_6311gdp
    e_delta_g2 = emp2_big + emp2_6311gd - emp2_6311g2dfp - emp2_6311pg_dp
    e_plus     = emp4_6311pg_dp - emp4_6311gd
    e_2df      = emp4_6311g2dfp - emp4_6311gd

    eg2 = eqci + e_delta_g2 + e_plus + e_2df
    eg2_mp2_0k = eqci + (emp2_big - emp2_6311gd) + hlc + zpe

    core.print_out('\n')
    core.print_out('  ==>  G1/G2 Energy Components  <==\n')
    core.print_out('\n')
    core.print_out('        QCISD(T):            %20.12lf\n' % eqci)
    core.print_out('        E(Delta):            %20.12lf\n' % e_delta_g2)
    core.print_out('        E(2DF):              %20.12lf\n' % e_2df)
    core.print_out('        E(+):                %20.12lf\n' % e_plus)
    core.print_out('        E(G1 HLC):           %20.12lf\n' % hlc1)
    core.print_out('        E(G2 HLC):           %20.12lf\n' % hlc)
    core.print_out('        E(ZPE):              %20.12lf\n' % zpe)
    core.print_out('\n')
    core.print_out('  ==>  0 Kelvin Results  <==\n')
    core.print_out('\n')
    eg2_0k = eg2 + zpe + hlc
    core.print_out('        G1:                  %20.12lf\n' % (eqci + e_plus + e_2df + hlc1 + zpe))
    core.print_out('        G2(MP2):             %20.12lf\n' % eg2_mp2_0k)
    core.print_out('        G2:                  %20.12lf\n' % eg2_0k)

    core.set_variable("G1 TOTAL ENERGY",eqci + e_plus + e_2df + hlc1 + zpe)
    core.set_variable("G2 TOTAL ENERGY",eg2_0k)
    core.set_variable("G2(MP2) TOTAL ENERGY",eg2_mp2_0k)

    core.print_out('\n')
    T = core.get_global_option('T')
    core.print_out('  ==>  %3.0lf Kelvin Results  <==\n'% T)
    core.print_out('\n')

    internal_energy = eg2_mp2_0k + du - zpe / 0.8929
    enthalpy        = eg2_mp2_0k + dh - zpe / 0.8929
    gibbs           = eg2_mp2_0k + dg - zpe / 0.8929

    core.print_out('        G2(MP2) energy:      %20.12lf\n' % internal_energy )
    core.print_out('        G2(MP2) enthalpy:    %20.12lf\n' % enthalpy)
    core.print_out('        G2(MP2) free energy: %20.12lf\n' % gibbs)
    core.print_out('\n')

    core.set_variable("G2(MP2) INTERNAL ENERGY",internal_energy)
    core.set_variable("G2(MP2) ENTHALPY",enthalpy)
    core.set_variable("G2(MP2) FREE ENERGY",gibbs)

    internal_energy = eg2_0k + du - zpe / 0.8929
    enthalpy        = eg2_0k + dh - zpe / 0.8929
    gibbs           = eg2_0k + dg - zpe / 0.8929

    core.print_out('        G2 energy:           %20.12lf\n' % internal_energy )
    core.print_out('        G2 enthalpy:         %20.12lf\n' % enthalpy)
    core.print_out('        G2 free energy:      %20.12lf\n' % gibbs)

    core.set_variable("CURRENT ENERGY",eg2_0k)

    core.set_variable("G2 INTERNAL ENERGY",internal_energy)
    core.set_variable("G2 ENTHALPY",enthalpy)
    core.set_variable("G2 FREE ENERGY",gibbs)

    core.clean()

    optstash.restore()

    # return 0K g2 results
    return eg2_0k

# aliases for g2
driver.procedures['energy']['gaussian-2'] = run_gaussian_2
driver.procedures['energy']['g2']         = run_gaussian_2
