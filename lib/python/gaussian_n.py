# Gn theory.  

import PsiMod
import re
import os
import math
import warnings
from driver import *
from molutil import *
from text import *
from procutil import *
from physconst import *
# never import aliases into this file

def run_gaussian_2(name, **kwargs):

    # throw an exception for open-shells
    if (PsiMod.get_option('SCF','REFERENCE') != 'RHF' ):
        raise ValidationError("""g2 computations require "reference rhf".""")

    # stash user options:
    optstash = OptionsState(
        ['FNOCC','COMPUTE_TRIPLES'],
        ['FNOCC','COMPUTE_MP4_TRIPLES'],
        ['FREEZE_CORE'],
        ['SCF','SCF_TYPE'])

    # override default scf_type
    PsiMod.set_local_option('SCF','SCF_TYPE','OUT_OF_CORE')

    # optimize geometry at scf level
    PsiMod.clean()
    PsiMod.set_global_option('BASIS',"6-31G(D)")
    optimize('scf')
    PsiMod.clean()

    # scf frequencies for zpe
    frequency('scf')

    # thermodynamic properties
    PsiMod.thermo()
    du = PsiMod.get_variable('DU')
    dh = PsiMod.get_variable('DH')
    dg = PsiMod.get_variable('DG')

    ref     = PsiMod.reference_wavefunction()
    freqs   = ref.frequencies()
    nfreq   = freqs.dim(0)
    freqsum = 0.0
    for i in range (0,nfreq):
        freqsum += freqs.get(i)
    zpe = freqsum / psi_hartree2wavenumbers * 0.8929 * 0.5
    PsiMod.clean()

    # optimize geometry at mp2 (no frozen core) level
    # note: freeze_core isn't an option in MP2
    PsiMod.set_global_option('FREEZE_CORE',"FALSE")
    optimize('conv-mp2')
    PsiMod.clean()

    # qcisd(t)
    PsiMod.set_local_option('FNOCC','COMPUTE_MP4_TRIPLES',"TRUE")
    PsiMod.set_global_option('FREEZE_CORE',"TRUE")
    PsiMod.set_global_option('BASIS',"6-311G(D_P)")
    run_fnocc('qcisd(t)',**kwargs)

    # HLC: high-level correction based on number of valence electrons
    ref    = PsiMod.reference_wavefunction()
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

    eqci_6311gdp = PsiMod.get_variable("QCISD(T) TOTAL ENERGY")
    emp4_6311gd  = PsiMod.get_variable("MP4 TOTAL ENERGY")
    emp2_6311gd  = PsiMod.get_variable("MP2 TOTAL ENERGY")
    PsiMod.clean()

    # correction for diffuse functions
    PsiMod.set_global_option('BASIS',"6-311+G(D_P)")
    energy('mp4')
    emp4_6311pg_dp = PsiMod.get_variable("MP4 TOTAL ENERGY")
    emp2_6311pg_dp = PsiMod.get_variable("MP2 TOTAL ENERGY")
    PsiMod.clean()

    # correction for polarization functions
    PsiMod.set_global_option('BASIS',"6-311G(2DF_P)")
    energy('mp4')
    emp4_6311g2dfp = PsiMod.get_variable("MP4 TOTAL ENERGY")
    emp2_6311g2dfp = PsiMod.get_variable("MP2 TOTAL ENERGY")
    PsiMod.clean()

    # big basis mp2
    PsiMod.set_global_option('BASIS',"6-311+G(3DF_2P)")
    run_fnocc('_mp2',**kwargs)
    emp2_big = PsiMod.get_variable("MP2 TOTAL ENERGY")
    PsiMod.clean()

    eqci       = eqci_6311gdp
    e_delta_g2 = emp2_big + emp2_6311gd - emp2_6311g2dfp - emp2_6311pg_dp
    e_plus     = emp4_6311pg_dp - emp4_6311gd
    e_2df      = emp4_6311g2dfp - emp4_6311gd

    eg2 = eqci + e_delta_g2 + e_plus + e_2df
    eg2_mp2_0k = eqci + (emp2_big - emp2_6311gd) + hlc + zpe

    PsiMod.print_out('\n')
    PsiMod.print_out('  ==>  G1/G2 Energy Components  <==\n')
    PsiMod.print_out('\n')
    PsiMod.print_out('        QCISD(T):            %20.12lf\n' % eqci)
    PsiMod.print_out('        E(Delta):            %20.12lf\n' % e_delta_g2)
    PsiMod.print_out('        E(2DF):              %20.12lf\n' % e_2df)
    PsiMod.print_out('        E(+):                %20.12lf\n' % e_plus)
    PsiMod.print_out('        E(G1 HLC):           %20.12lf\n' % hlc1)
    PsiMod.print_out('        E(G2 HLC):           %20.12lf\n' % hlc)
    PsiMod.print_out('        E(ZPE):              %20.12lf\n' % zpe)
    PsiMod.print_out('\n')
    PsiMod.print_out('  ==>  0 Kelvin Results  <==\n')
    PsiMod.print_out('\n')
    eg2_0k = eg2 + zpe + hlc
    PsiMod.print_out('        G1:                  %20.12lf\n' % (eqci + e_plus + e_2df + hlc1 + zpe))
    PsiMod.print_out('        G2(MP2):             %20.12lf\n' % eg2_mp2_0k)
    PsiMod.print_out('        G2:                  %20.12lf\n' % eg2_0k)

    PsiMod.set_variable("G1 TOTAL ENERGY",eqci + e_plus + e_2df + hlc1 + zpe)
    PsiMod.set_variable("G2 TOTAL ENERGY",eg2_0k)
    PsiMod.set_variable("G2(MP2) TOTAL ENERGY",eg2_mp2_0k)

    PsiMod.print_out('\n')
    T = PsiMod.get_global_option('T')
    PsiMod.print_out('  ==>  %3.0lf Kelvin Results  <==\n'% T)
    PsiMod.print_out('\n')

    internal_energy = eg2_mp2_0k + du - zpe / 0.8929
    enthalpy        = eg2_mp2_0k + dh - zpe / 0.8929
    gibbs           = eg2_mp2_0k + dg - zpe / 0.8929

    PsiMod.print_out('        G2(MP2) energy:      %20.12lf\n' % internal_energy )
    PsiMod.print_out('        G2(MP2) enthalpy:    %20.12lf\n' % enthalpy)
    PsiMod.print_out('        G2(MP2) free energy: %20.12lf\n' % gibbs)
    PsiMod.print_out('\n')

    PsiMod.set_variable("G2(MP2) INTERNAL ENERGY",internal_energy)
    PsiMod.set_variable("G2(MP2) ENTHALPY",enthalpy)
    PsiMod.set_variable("G2(MP2) FREE ENERGY",gibbs)

    internal_energy = eg2_0k + du - zpe / 0.8929
    enthalpy        = eg2_0k + dh - zpe / 0.8929
    gibbs           = eg2_0k + dg - zpe / 0.8929

    PsiMod.print_out('        G2 energy:           %20.12lf\n' % internal_energy )
    PsiMod.print_out('        G2 enthalpy:         %20.12lf\n' % enthalpy)
    PsiMod.print_out('        G2 free energy:      %20.12lf\n' % gibbs)

    PsiMod.set_variable("CURRENT ENERGY",eg2_0k)

    PsiMod.set_variable("G2 INTERNAL ENERGY",internal_energy)
    PsiMod.set_variable("G2 ENTHALPY",enthalpy)
    PsiMod.set_variable("G2 FREE ENERGY",gibbs)

    PsiMod.clean()

    optstash.restore()

    # return 0K g2 results
    return eg2_0k

# aliases for g2
procedures['energy']['gaussian-2'] = run_gaussian_2
procedures['energy']['g2']         = run_gaussian_2
