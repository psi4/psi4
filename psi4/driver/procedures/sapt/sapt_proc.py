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

import psi4
import numpy as np
from psi4 import core

from psi4 import extras
from psi4.driver import p4util
from psi4.driver import qcdb
from psi4.driver.p4util.exceptions import *
from psi4.driver.molutil import *
from psi4.driver.procedures.proc import scf_helper

from . import sapt_jk_terms

# Only export the run_ scripts
__all__ = ['run_sapt_dft']

def run_sapt_dft(name, **kwargs):
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'REFERENCE'],
        ['SCF', 'DFT_FUNCTIONAL'],
        ['SCF', 'DFT_GRAC_SHIFT'],
        ['SCF', 'SAVE_JK'])

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

    # Grab the overall data
    mon_a_shift = core.get_option("SAPT", "SAPT_DFT_GRAC_SHIFT_A")
    mon_b_shift = core.get_option("SAPT", "SAPT_DFT_GRAC_SHIFT_B")
    do_delta_hf = core.get_option("SAPT", "SAPT_DFT_DO_DHF")

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
    # data["SCF DIMER"] = psi4.core.get_variable("CURRENT ENERGY")

    # Set the primary functional
    core.set_global_option("DFT_FUNCTIONAL", core.get_option("SAPT", "SAPT_DFT_FUNCTIONAL"))
    core.set_local_option('SCF', 'REFERENCE', 'RKS')

    if (core.get_option('SCF', 'SCF_TYPE') == 'DF'):
        core.set_global_option('DF_INTS_IO', 'LOAD')
        core.set_global_option('DF_INTS_IO', 'SAVE')

    # Compute Monomer A wavefunction
    if (core.get_option('SCF', 'SCF_TYPE') == 'DF'):
        core.IO.change_file_namespace(97, 'dimer', 'monomerA')

    if mon_a_shift:
        core.set_global_option("DFT_GRAC_SHIFT", mon_a_shift)

    # Save the JK object
    core.set_global_option("SAVE_JK", True)

    core.IO.set_default_namespace('monomerA')
    core.print_out('\n')
    p4util.banner('Monomer A HF')
    core.print_out('\n')
    wfn_A = scf_helper("SCF", molecule=monomerA, **kwargs)
    data["SCF MONOMERA"] = psi4.core.get_variable("CURRENT ENERGY")

    core.set_global_option("DFT_GRAC_SHIFT", 0.0)

    # Compute Monomer B wavefunction
    if (core.get_option('SCF', 'SCF_TYPE') == 'DF'):
        core.IO.change_file_namespace(97, 'monomerA', 'monomerB')

    if mon_b_shift:
        core.set_global_option("DFT_GRAC_SHIFT", mon_b_shift)


    core.IO.set_default_namespace('monomerB')
    core.print_out('\n')
    p4util.banner('Monomer B HF')
    core.print_out('\n')
    wfn_B = scf_helper("SCF", molecule=monomerB, **kwargs)
    data["SCF MONOMERB"] = psi4.core.get_variable("CURRENT ENERGY")

    core.set_global_option("DFT_GRAC_SHIFT", 0.0)

    # Compute SAPT interaction

    # Build cache and JK
    sapt_jk = wfn_B.jk()

    cache = sapt_jk_terms.build_sapt_jk_cache(wfn_A, wfn_B, sapt_jk, True)


    elst = sapt_jk_terms.electrostatics(cache, True)
    exch = sapt_jk_terms.exchange(cache, sapt_jk, True)
    ind = sapt_jk_terms.induction(cache, sapt_jk, True)


