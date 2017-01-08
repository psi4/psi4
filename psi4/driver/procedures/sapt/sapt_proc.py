#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
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

import numpy as np

from psi4 import core
from psi4.driver import p4util
from psi4.driver.p4util.exceptions import *
from psi4.driver.molutil import *
from psi4.driver.procedures.proc import scf_helper

from . import sapt_jk_terms
from .sapt_util import print_sapt_summary
from . import sapt_mp2_terms

# Only export the run_ scripts
__all__ = ['run_sapt_dft']

def run_sapt_dft(name, **kwargs):
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'REFERENCE'],
        ['SCF', 'DFT_FUNCTIONAL'],
        ['SCF', 'DFT_GRAC_SHIFT'],
        ['SCF', 'SAVE_JK'])

    core.tstart()
    # Alter default algorithm
    if not core.has_option_changed('SCF', 'SCF_TYPE'):
        core.set_local_option('SCF', 'SCF_TYPE', 'DF')

    core.prepare_options_for_module("SAPT")

    # Get the molecule of interest
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        sapt_dimer = kwargs.pop('molecule', core.get_active_molecule())
    else:
        core.print_out('Warning! SAPT argument "ref_wfn" is only able to use molecule information.')
        sapt_dimer = ref_wfn.molecule()

    # Shifting to C1 so we need to copy the active molecule
    if sapt_dimer.schoenflies_symbol() != 'c1':
        core.print_out('  SAPT does not make use of molecular symmetry, further calculations in C1 point group.\n')

    # Make sure the geometry doesnt shift or rotate
    sapt_dimer = sapt_dimer.clone()
    sapt_dimer.reset_point_group('c1')
    sapt_dimer.fix_orientation(True)
    sapt_dimer.fix_com(True)
    sapt_dimer.update_geometry()

    # Grab overall settings
    mon_a_shift = core.get_option("SAPT", "SAPT_DFT_GRAC_SHIFT_A")
    mon_b_shift = core.get_option("SAPT", "SAPT_DFT_GRAC_SHIFT_B")
    do_delta_hf = core.get_option("SAPT", "SAPT_DFT_DO_DHF")
    sapt_dft_functional = core.get_option("SAPT", "SAPT_DFT_FUNCTIONAL")

    # Print out the title and some information
    core.print_out("\n");
    core.print_out("         ---------------------------------------------------------\n");
    core.print_out("         " + "SAPT(DFT) Procedure".center(58) + "\n");
    core.print_out("\n");
    core.print_out("         " + "by Daniel G. A. Smith".center(58) + "\n");
    core.print_out("         ---------------------------------------------------------\n");
    core.print_out("\n");

    core.print_out("  ==> Algorithm <==\n\n")
    core.print_out("   SAPT DFT Functional     %12s\n" % str(sapt_dft_functional))
    core.print_out("   Monomer A GRAC Shift    %12.6f\n" % mon_a_shift)
    core.print_out("   Monomer B GRAC Shift    %12.6f\n" % mon_b_shift)
    core.print_out("   Delta HF                %12s\n" % ("True" if do_delta_hf else "False"))
    core.print_out("   JK Algorithm            %12s\n" % core.get_option("SCF", "SCF_TYPE"))
    core.print_out("\n")
    core.print_out("   Required computations:\n")
    if (do_delta_hf):
        core.print_out("     HF (Dimer)\n")
        core.print_out("     HF (Monomer A)\n")
        core.print_out("     HF (Monomer B)\n")
    core.print_out("     DFT (Monomer A)\n")
    core.print_out("     DFT (Monomer B)\n")
    core.print_out("\n")



    if (mon_a_shift == 0.0) or (mon_b_shift == 0.0):
        raise ValidationError('SAPT(DFT): must set both "SAPT_DFT_GRAC_SHIFT_A" and "B".')

    if (core.get_option('SCF', 'REFERENCE') != 'RHF'):
        raise ValidationError('SAPT(DFT) currently only supports restricted references.')

    nfrag = sapt_dimer.nfragments()
    if nfrag != 2:
        raise ValidationError('SAPT requires active molecule to have 2 fragments, not %s.' % (nfrag))

    monomerA = sapt_dimer.extract_subsets(1, 2)
    monomerA.set_name('monomerA')
    monomerB = sapt_dimer.extract_subsets(2, 1)
    monomerB.set_name('monomerB')

    core.IO.set_default_namespace('dimer')
    data = {}

    # # Compute dimer wavefunction
    # core.print_out('\n')
    # p4util.banner('Dimer HF')
    # core.print_out('\n')
    # if (core.get_option('SCF', 'SCF_TYPE') == 'DF'):
    #     core.set_global_option('DF_INTS_IO', 'SAVE')

    # dimer_wfn = scf_helper("SCF", molecule=sapt_dimer, **kwargs)
    # data["SCF DIMER"] = core.get_variable("CURRENT ENERGY")

    # Set the primary functional
    core.set_global_option("DFT_FUNCTIONAL", core.get_option("SAPT", "SAPT_DFT_FUNCTIONAL"))
    core.set_local_option('SCF', 'REFERENCE', 'RKS')

    if (core.get_option('SCF', 'SCF_TYPE') == 'DF'):
        # core.set_global_option('DF_INTS_IO', 'LOAD')
        core.set_global_option('DF_INTS_IO', 'SAVE')

    # Compute Monomer A wavefunction
    if (core.get_option('SCF', 'SCF_TYPE') == 'DF'):
        core.IO.change_file_namespace(97, 'dimer', 'monomerA')

    if mon_a_shift:
        core.set_global_option("DFT_GRAC_SHIFT", mon_a_shift)

    # Save the JK object
    core.set_global_option("SAVE_JK", True)
    core.IO.set_default_namespace('monomerA')
    wfn_A = scf_helper("SCF", molecule=monomerA, banner="SAPT(DFT): DFT Monomer A", **kwargs)
    data["SCF MONOMERA"] = core.get_variable("CURRENT ENERGY")

    core.set_global_option("DFT_GRAC_SHIFT", 0.0)

    # Compute Monomer B wavefunction
    if (core.get_option('SCF', 'SCF_TYPE') == 'DF'):
        core.IO.change_file_namespace(97, 'monomerA', 'monomerB')

    if mon_b_shift:
        core.set_global_option("DFT_GRAC_SHIFT", mon_b_shift)


    core.IO.set_default_namespace('monomerB')
    wfn_B = scf_helper("SCF", molecule=monomerB, banner="SAPT(DFT): DFT Monomer B", **kwargs)
    data["SCF MONOMERB"] = core.get_variable("CURRENT ENERGY")

    core.set_global_option("DFT_GRAC_SHIFT", 0.0)

    # Print out the title and some information
    core.print_out("\n");
    core.print_out("         ---------------------------------------------------------\n");
    core.print_out("         " + "SAPT(DFT): Intermolecular Interaction Segment".center(58) + "\n");
    core.print_out("\n");
    core.print_out("         " + "by Daniel G. A. Smith and Rob Parrish".center(58) + "\n");
    core.print_out("         ---------------------------------------------------------\n");
    core.print_out("\n");

    core.print_out("  ==> Algorithm <==\n\n")
    core.print_out("   SAPT DFT Functional     %12s\n" % str(sapt_dft_functional))
    core.print_out("   Monomer A GRAC Shift    %12.6f\n" % mon_a_shift)
    core.print_out("   Monomer B GRAC Shift    %12.6f\n" % mon_b_shift)
    core.print_out("   Delta HF                %12s\n" % ("True" if do_delta_hf else "False"))
    core.print_out("   JK Algorithm            %12s\n" % core.get_option("SCF", "SCF_TYPE"))

    # Build cache and JK
    sapt_jk = wfn_B.jk()

    cache = sapt_jk_terms.build_sapt_jk_cache(wfn_A, wfn_B, sapt_jk, True)

    # Electostatics
    elst = sapt_jk_terms.electrostatics(cache, True)
    data.update(elst)

    # Exchange
    exch = sapt_jk_terms.exchange(cache, sapt_jk, True)
    data.update(exch)

    # Induction
    ind = sapt_jk_terms.induction(cache, sapt_jk, True,
                                  maxiter=core.get_option("SAPT", "MAXITER"),
                                  conv=core.get_option("SAPT", "D_CONVERGENCE"))
    data.update(ind)

    # Dispersion
    primary_basis = wfn_A.basisset()
    aux_basis = core.BasisSet.build(sapt_dimer, "DF_BASIS_MP2",
                                    core.get_option("DFMP2", "DF_BASIS_MP2"),
                                    "RIFIT", core.get_global_option('BASIS'))
    fdds_disp = sapt_mp2_terms.df_fdds_dispersion(primary_basis, aux_basis, cache)

    mp2_disp = sapt_mp2_terms.df_mp2_dispersion(wfn_A, primary_basis, aux_basis, cache, do_print=True)

    # Print out final data
    core.print_out("\n")
    core.print_out(print_sapt_summary(data, "SAPT(DFT)"))


    core.tstop()

    return None

