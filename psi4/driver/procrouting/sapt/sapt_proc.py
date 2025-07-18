#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import numpy as np

from psi4 import core

from ... import p4util
from ...constants import constants
from ...p4util.exceptions import *
from .. import proc_util
from ..proc import scf_helper, run_scf
from . import sapt_jk_terms, sapt_mp2_terms, sapt_sf_terms
from .sapt_util import print_sapt_dft_summary, print_sapt_hf_summary, print_sapt_var
import qcelemental as qcel
from ...p4util.exceptions import ConvergenceError

# Only export the run_ scripts
__all__ = ["run_sapt_dft", "sapt_dft", "run_sf_sapt"]


def run_sapt_dft(name, **kwargs):
    optstash = p4util.OptionsState(
        ["SCF_TYPE"],
        ["SCF", "REFERENCE"],
        ["SCF", "DFT_GRAC_SHIFT"],
        ["SCF", "SAVE_JK"],
    )

    core.tstart()
    # Alter default algorithm
    if not core.has_global_option_changed("SCF_TYPE"):
        core.set_global_option("SCF_TYPE", "DF")

    core.prepare_options_for_module("SAPT")

    # Grab overall settings
    do_mon_grac_shift_A = False
    do_mon_grac_shift_B = False
    mon_a_shift = core.get_option("SAPT", "SAPT_DFT_GRAC_SHIFT_A")
    mon_b_shift = core.get_option("SAPT", "SAPT_DFT_GRAC_SHIFT_B")
    grac_compute = core.get_option("SAPT", "SAPT_DFT_GRAC_COMPUTE")

    if not core.has_option_changed("SAPT", "SAPT_DFT_GRAC_SHIFT_A") and grac_compute.upper() != "NONE":
        do_mon_grac_shift_A = True
    if not core.has_option_changed("SAPT", "SAPT_DFT_GRAC_SHIFT_B") and grac_compute.upper() != "NONE":
        do_mon_grac_shift_B = True

    do_delta_hf = core.get_option("SAPT", "SAPT_DFT_DO_DHF")
    sapt_dft_functional = core.get_option("SAPT", "SAPT_DFT_FUNCTIONAL")
    do_dft = sapt_dft_functional != "HF"

    # Get the molecule of interest
    ref_wfn = kwargs.get("ref_wfn", None)
    if ref_wfn is None:
        sapt_dimer = kwargs.pop("molecule", core.get_active_molecule())
    else:
        core.print_out(
            'Warning! SAPT argument "ref_wfn" is only able to use molecule information.'
        )
        sapt_dimer = ref_wfn.molecule()

    sapt_dimer, monomerA, monomerB = proc_util.prepare_sapt_molecule(
        sapt_dimer, "dimer"
    )

    if do_mon_grac_shift_A or do_mon_grac_shift_B:
        monomerA_mon_only_bf = sapt_dimer.extract_subsets(1)
        monomerB_mon_only_bf = sapt_dimer.extract_subsets(2)

    # Print out the title and some information
    core.print_out("\n")
    core.print_out(
        "         ---------------------------------------------------------\n"
    )
    core.print_out("         " + "SAPT(DFT) Procedure".center(58) + "\n")
    core.print_out("\n")
    core.print_out("         " + "by Daniel G. A. Smith".center(58) + "\n")
    core.print_out(
        "         ---------------------------------------------------------\n"
    )
    core.print_out("\n")

    # core.print_out("  !!!  WARNING:  SAPT(DFT) capability is in beta. Please use with caution. !!!\n\n")
    core.print_out(
        "Warning! The default value of SAPT_DFT_EXCH_DISP_SCALE_SCHEME has changed from DISP to FIXED. Please be careful comparing results with earlier versions. \n\n"
    )

    core.print_out("  ==> Algorithm <==\n\n")
    core.print_out("   SAPT DFT Functional     %12s\n" %
                   str(sapt_dft_functional))
    core.print_out("   Monomer A GRAC Shift    %12.6f\n" % mon_a_shift)
    core.print_out("   Monomer B GRAC Shift    %12.6f\n" % mon_b_shift)
    # fmt: off
    core.print_out("   Delta HF                %12s\n" % ("True" if do_delta_hf else "False"))
    core.print_out("   JK Algorithm            %12s\n" % core.get_global_option("SCF_TYPE"))
    # fmt: on
    core.print_out("\n")
    core.print_out("   Required computations:\n")
    if do_delta_hf:
        core.print_out("     HF   (Dimer)\n")
        core.print_out("     HF   (Monomer A)\n")
        core.print_out("     HF   (Monomer B)\n")
    if do_dft:
        core.print_out("     DFT  (Monomer A)\n")
        core.print_out("     DFT  (Monomer B)\n")
    if do_mon_grac_shift_A:
        core.print_out("     GRAC (Monomer A)\n")
        mon_a_shift = compute_GRAC_shift(
            monomerA_mon_only_bf,
            grac_compute,
            "Monomer A",
        )
    if do_mon_grac_shift_B:
        core.print_out("     GRAC (Monomer B)\n")
        mon_b_shift = compute_GRAC_shift(
            monomerB_mon_only_bf,
            grac_compute,
            "Monomer B",
        )

    core.set_variable("SAPT DFT GRAC SHIFT A", mon_a_shift)  # P::e SAPT
    core.set_variable("SAPT DFT GRAC SHIFT B", mon_b_shift)  # P::e SAPT
    core.print_out("\n")

    if do_dft and ((not core.has_option_changed("SAPT", "SAPT_DFT_GRAC_SHIFT_A")) or (not core.has_option_changed("SAPT", "SAPT_DFT_GRAC_SHIFT_B"))) and grac_compute == "NONE":
        raise ValidationError(
            'SAPT(DFT): User must set both "SAPT_DFT_GRAC_SHIFT_A" and "_B".  Or, to automatically compute the GRAC shift, set SAPT_DFT_GRAC_COMPUTE to "ITERATIVE" or "SINGLE".'
        )

    if core.get_option("SCF", "REFERENCE") != "RHF":
        raise ValidationError(
            "SAPT(DFT) currently only supports restricted references."
        )

    core.IO.set_default_namespace("dimer")
    data = {}

    # Save integrals
    # We want to try to re-use itegrals for the dimer and monomer SCF's. If we
    # are using Disk based DF (DISK_DF) then we can use the DF_INTS_IO option.
    # MemDF does not know about this option but setting it will be harmless
    # there.
    core.set_global_option('DF_INTS_IO', 'SAVE')

    # Compute dimer wavefunction
    hf_wfn_dimer = None
    if do_delta_hf:
        core.set_global_option('DF_INTS_IO', 'SAVE')
        core.timer_on("SAPT(DFT):Dimer SCF")
        hf_data = {}

        core.set_local_option("SCF", "SAVE_JK", True)
        hf_wfn_dimer = scf_helper("SCF", molecule=sapt_dimer, banner="SAPT(DFT): delta HF Dimer", **kwargs)
        hf_data["HF DIMER"] = core.variable("CURRENT ENERGY")
        core.timer_off("SAPT(DFT):Dimer SCF")

        core.timer_on("SAPT(DFT):Monomer A SCF")
        core.IO.change_file_namespace(97, 'dimer', 'monomerA')

        jk_obj = hf_wfn_dimer.jk()
        hf_wfn_A = scf_helper("SCF", molecule=monomerA, banner="SAPT(DFT): delta HF Monomer A", jk=jk_obj, **kwargs)
        hf_data["HF MONOMER A"] = core.variable("CURRENT ENERGY")
        core.timer_off("SAPT(DFT):Monomer A SCF")

        core.timer_on("SAPT(DFT):Monomer B SCF")
        core.IO.change_file_namespace(97, 'monomerA', 'monomerB')

        hf_wfn_B = scf_helper("SCF", molecule=monomerB, banner="SAPT(DFT): delta HF Monomer B", jk=jk_obj, **kwargs)
        hf_data["HF MONOMER B"] = core.variable("CURRENT ENERGY")
        core.set_global_option("SAVE_JK", False)
        core.timer_off("SAPT(DFT):Monomer B SCF")

        if do_dft:  # For SAPT(HF) do the JK terms in sapt_dft()
            # Grab JK object and set to A (so we do not save many JK objects)
            sapt_jk = hf_wfn_B.jk()
            hf_wfn_A.set_jk(sapt_jk)
            core.set_global_option("SAVE_JK", False)

            # Move it back to monomer A
            core.IO.change_file_namespace(97, 'monomerB', 'dimer')

            core.print_out("\n")
            core.print_out(
                "         ---------------------------------------------------------\n"
            )
            core.print_out(
                "         " + "SAPT(DFT): delta HF Segment".center(58) + "\n"
            )
            core.print_out("\n")
            core.print_out(
                "         " +
                "by Daniel G. A. Smith and Rob Parrish".center(58) + "\n"
            )
            core.print_out(
                "         ---------------------------------------------------------\n"
            )
            core.print_out("\n")

            # Build cache
            hf_cache = sapt_jk_terms.build_sapt_jk_cache(
                hf_wfn_A, hf_wfn_B, sapt_jk, True
            )

            # Electrostatics
            core.timer_on("SAPT(HF):elst")
            elst = sapt_jk_terms.electrostatics(hf_cache, True)
            hf_data.update(elst)
            core.timer_off("SAPT(HF):elst")

            # Exchange
            core.timer_on("SAPT(HF):exch")
            exch = sapt_jk_terms.exchange(hf_cache, sapt_jk, True)
            hf_data.update(exch)
            core.timer_off("SAPT(HF):exch")

            # Induction
            core.timer_on("SAPT(HF):ind")
            ind = sapt_jk_terms.induction(
                hf_cache,
                sapt_jk,
                True,
                maxiter=core.get_option("SAPT", "MAXITER"),
                conv=core.get_option("SAPT", "CPHF_R_CONVERGENCE"),
                Sinf=core.get_option("SAPT", "DO_IND_EXCH_SINF"),
            )
            hf_data.update(ind)
            core.timer_off("SAPT(HF):ind")

            dhf_value = (
                hf_data["HF DIMER"] - hf_data["HF MONOMER A"] -
                hf_data["HF MONOMER B"]
            )

            core.print_out("\n")
            core.print_out(
                print_sapt_hf_summary(hf_data, "SAPT(HF)", delta_hf=dhf_value)
            )

            data["Delta HF Correction"] = core.variable("SAPT(DFT) Delta HF")
            sapt_jk.finalize()

            del hf_wfn_A, hf_wfn_B, sapt_jk

        else:
            wfn_A = hf_wfn_A
            wfn_B = hf_wfn_B
            data["DFT MONOMER A"] = hf_data["HF MONOMER A"]
            data["DFT MONOMER B"] = hf_data["HF MONOMER B"]
            dhf_value = (
                hf_data["HF DIMER"] - hf_data["HF MONOMER A"] -
                hf_data["HF MONOMER B"]
            )
            data["DHF VALUE"] = dhf_value

    if hf_wfn_dimer is None:
        dimer_wfn = core.Wavefunction.build(
            sapt_dimer, core.get_global_option("BASIS"))
    else:
        dimer_wfn = hf_wfn_dimer

    if do_dft or not do_delta_hf:

        # Set the primary functional
        core.set_local_option("SCF", "REFERENCE", "RKS")

        # Compute Monomer A wavefunction
        core.timer_on("SAPT(DFT): Monomer A DFT")
        core.IO.change_file_namespace(97, 'dimer', 'monomerA')

        if mon_a_shift:
            core.set_global_option("DFT_GRAC_SHIFT", mon_a_shift)

        core.IO.set_default_namespace('monomerA')
        core.set_global_option("SAVE_JK", True)
        wfn_A = scf_helper(sapt_dft_functional,
                           post_scf=False,
                           molecule=monomerA,
                           banner="SAPT(DFT): DFT Monomer A",
                           **kwargs)
        data["DFT MONOMERA"] = core.variable("CURRENT ENERGY")

        core.set_global_option("DFT_GRAC_SHIFT", 0.0)
        core.timer_off("SAPT(DFT): Monomer A DFT")

        # Compute Monomer B wavefunction
        core.timer_on("SAPT(DFT): Monomer B DFT")
        core.IO.change_file_namespace(97, 'monomerA', 'monomerB')

        if mon_b_shift:
            core.set_global_option("DFT_GRAC_SHIFT", mon_b_shift)

        core.set_global_option("SAVE_JK", True)
        core.IO.set_default_namespace('monomerB')
        wfn_B = scf_helper(sapt_dft_functional,
                           post_scf=False,
                           molecule=monomerB,
                           banner="SAPT(DFT): DFT Monomer B",
                           jk=wfn_A.jk(),
                           **kwargs)
        data["DFT MONOMERB"] = core.variable("CURRENT ENERGY")
        core.timer_off("SAPT(DFT): Monomer B DFT")

    # Save JK object
    sapt_jk = wfn_B.jk()
    wfn_A.set_jk(sapt_jk)
    core.set_global_option("SAVE_JK", False)

    core.set_global_option("DFT_GRAC_SHIFT", 0.0)

    # Write out header
    scf_alg = core.get_global_option("SCF_TYPE")
    sapt_dft_header(
        sapt_dft_functional, mon_a_shift, mon_b_shift, bool(
            do_delta_hf), scf_alg
    )

    # Compute Delta HF for SAPT(HF)?
    delta_hf = do_delta_hf and not do_dft

    # Call SAPT(DFT)
    sapt_jk = wfn_B.jk()
    sapt_dft(
        dimer_wfn,
        wfn_A,
        wfn_B,
        do_dft=do_dft,
        sapt_jk=sapt_jk,
        data=data,
        print_header=False,
        delta_hf=delta_hf,
    )

    # Copy data back into globals
    for k, v in data.items():
        core.set_variable(k, v)

    core.tstop()

    optstash.restore()

    return dimer_wfn


sapt_dft_grac_convergence_tier_options = {
    "SINGLE": [
        {
            "SCF_INITIAL_ACCELERATOR": "ADIIS",
        }
    ],
    "ITERATIVE": [
        {
            "SCF_INITIAL_ACCELERATOR": "ADIIS",
        },
        {
            "LEVEL_SHIFT": 0.1,
            "LEVEL_SHIFT_CUTOFF": 1e-05,
            "SCF_INITIAL_ACCELERATOR": "ADIIS",
            "MAXITER": 200,
        },
        {
            "LEVEL_SHIFT": 0.5,
            "LEVEL_SHIFT_CUTOFF": 1e-3,
            "SCF_INITIAL_ACCELERATOR": "ADIIS",
            "MAXITER": 200,
        },
    ],
}


def compute_GRAC_shift(molecule, sapt_dft_grac_convergence_tier, label):
    optstash = p4util.OptionsState(
        ["SCF_TYPE"],
        ["SCF", "REFERENCE"],
        ["SCF", "DFT_GRAC_SHIFT"],
        ["SCF", "SAVE_JK"],
        ["SCF", "MAXITER"],
        ["SCF", "LEVEL_SHIFT"],
        ["SCF", "LEVEL_SHIFT_CUTOFF"],
        ["SCF", "SCF_INITIAL_ACCELERATOR"],
        ["BASIS"],
    )

    dft_functional = core.get_option("SAPT", "SAPT_DFT_FUNCTIONAL")
    grac_basis = core.get_option("SAPT", "SAPT_DFT_GRAC_BASIS")
    if grac_basis != "AUTO":
        core.set_global_option("BASIS", grac_basis)

    core.print_out(f"Computing GRAC shift for {label} using {sapt_dft_grac_convergence_tier}...")
    grac_options = sapt_dft_grac_convergence_tier_options[
        sapt_dft_grac_convergence_tier
    ]
    for options in grac_options:
        for key, val in options.items():
            core.set_local_option("SCF", key, val)
        # Need to get the initial and cation to estimate ionization energy for
        # GRAC shift
        mol_qcel_dict = molecule.to_schema(dtype=2)
        del mol_qcel_dict["fragment_charges"]
        del mol_qcel_dict["fragment_multiplicities"]
        del mol_qcel_dict["molecular_multiplicity"]
        given_charge = mol_qcel_dict["molecular_charge"]

        mol_qcel = qcel.models.Molecule(**mol_qcel_dict)
        mol_given = core.Molecule.from_schema(mol_qcel.dict())

        mol_qcel_dict["molecular_charge"] += 1
        mol_qcel = qcel.models.Molecule(**mol_qcel_dict)
        mol_cation = core.Molecule.from_schema(mol_qcel.dict())

        core.print_out(f"\n\n  ==> GRAC {label} Given Molecule: charge={mol_given.molecular_charge()} mult={mol_given.multiplicity()} <==\n\n")
        try:
            if mol_given.multiplicity() != 1:
                core.set_local_option("SCF", "REFERENCE", "UHF")
            else:
                core.set_local_option("SCF", "REFERENCE", "RHF")
            wfn_given = run_scf(
                dft_functional.lower(),
                molecule=mol_given,
            )
            if mol_cation.multiplicity() != 1:
                core.set_local_option("SCF", "REFERENCE", "UHF")
            else:
                core.set_local_option("SCF", "REFERENCE", "RHF")
            core.print_out(f"\n\n  ==> GRAC {label} Electron Removed Molecule: charge={mol_cation.molecular_charge()} mult={mol_cation.multiplicity()} <==\n\n")
            wfn_cation = run_scf(
                dft_functional.lower(),
                molecule=mol_cation,
            )
        except ConvergenceError:
            if len(options) == 1:
                raise Exception(
                    "Convergence error in GRAC shift calculation, please try a different convergence tier."
                )
            else:
                core.print_out("Convergence error, trying next GRAC iteration...")
            continue
        occ_given = wfn_given.epsilon_a_subset(basis="SO", subset="OCC").to_array(
            dense=True
        )
        HOMO = np.amax(occ_given)

        E_given = wfn_given.energy()
        E_cation = wfn_cation.energy()
        grac = E_cation - E_given + HOMO
        if grac >= 1 or grac <= -1:
            raise Exception(f"The computed GRAC shift ({grac} [E_h]) for {label} exceeds the bounds of -1 < x < 1 and should not be used to approximate the ionization potential.")
    optstash.restore()
    return grac


def sapt_dft_header(
    sapt_dft_functional="unknown",
    mon_a_shift=None,
    mon_b_shift=None,
    do_delta_hf="N/A",
    jk_alg="N/A",
):
    # Print out the title and some information
    core.print_out("\n")
    core.print_out(
        "         ---------------------------------------------------------\n"
    )
    core.print_out(
        "         " +
        "SAPT(DFT): Intermolecular Interaction Segment".center(58) + "\n"
    )
    core.print_out("\n")
    core.print_out(
        "         " + "by Daniel G. A. Smith and Rob Parrish".center(58) + "\n"
    )
    core.print_out(
        "         ---------------------------------------------------------\n"
    )
    core.print_out("\n")

    core.print_out("  ==> Algorithm <==\n\n")
    core.print_out("   SAPT DFT Functional     %12s\n" %
                   str(sapt_dft_functional))
    if mon_a_shift:
        core.print_out("   Monomer A GRAC Shift    %12.6f\n" % mon_a_shift)
    if mon_b_shift:
        core.print_out("   Monomer B GRAC Shift    %12.6f\n" % mon_b_shift)
    core.print_out("   Delta HF                %12s\n" % do_delta_hf)
    core.print_out("   JK Algorithm            %12s\n" % jk_alg)


def sapt_dft(
    dimer_wfn,
    wfn_A,
    wfn_B,
    do_dft=True,
    sapt_jk=None,
    sapt_jk_B=None,
    data=None,
    print_header=True,
    cleanup_jk=True,
    delta_hf=False,
):
    """
    The primary SAPT(DFT) algorithm to compute the interaction energy once the wavefunctions have been built.

    Example
    -------

    dimer = psi4.geometry('''
      Ne
      --
      Ar 1 6.5
      units bohr
    ''')

    psi4.set_options({"BASIS": "aug-cc-pVDZ"})

    # Prepare the fragments
    sapt_dimer, monomerA, monomerB = psi4.proc_util.prepare_sapt_molecule(sapt_dimer, "dimer")

    # Run the first monomer
    set DFT_GRAC_SHIFT 0.203293
    wfnA, energyA = psi4.energy("PBE0", monomer=monomerA, return_wfn=True)

    # Run the second monomer
    set DFT_GRAC_SHIFT 0.138264
    wfnB, energyB = psi4.energy("PBE0", monomer=monomerB, return_wfn=True)

    # Build the dimer wavefunction
    wfnD = psi4.core.Wavefunction.build(sapt_dimer)

    # Compute SAPT(DFT) from the provided wavefunctions
    data = psi4.procrouting.sapt.sapt_dft(wfnD, wfnA, wfnB)
    """

    # Handle the input options
    core.timer_on("SAPT(DFT):Build JK")
    if print_header:
        sapt_dft_header()

    if sapt_jk is None:

        core.print_out("\n   => Building SAPT JK object <= \n\n")
        sapt_jk = core.JK.build(dimer_wfn.basisset())
        sapt_jk.set_do_J(True)
        sapt_jk.set_do_K(True)
        if wfn_A.functional().is_x_lrc():
            sapt_jk.set_do_wK(True)
            sapt_jk.set_omega(wfn_A.functional().x_omega())
        sapt_jk.initialize()
        sapt_jk.print_header()
        if wfn_B.functional().is_x_lrc() and (wfn_A.functional().x_omega() != wfn_B.functional().x_omega()):
            core.print_out("   => Monomer B: Building SAPT JK object <= \n\n")
            core.print_out(
                "      Reason: MonomerA Omega != MonomerB Omega\n\n")
            sapt_jk_B = core.JK.build(dimer_wfn.basisset())
            sapt_jk_B.set_do_J(True)
            sapt_jk_B.set_do_K(True)
            sapt_jk_B.set_do_wK(True)
            sapt_jk_B.set_omega(wfn_B.functional().x_omega())
            sapt_jk_B.initialize()
            sapt_jk_B.print_header()

    else:
        sapt_jk.set_do_K(True)

    sapt_jk.set_do_J(True)
    sapt_jk.set_do_K(True)

    if wfn_A.functional().is_x_lrc():
        sapt_jk.set_do_wK(True)
        sapt_jk.set_omega(wfn_A.functional().x_omega())


    if data is None:
        data = {}

    # Build SAPT cache
    cache = sapt_jk_terms.build_sapt_jk_cache(wfn_A, wfn_B, sapt_jk, True)
    core.timer_off("SAPT(DFT):Build JK")

    # Electrostatics
    core.timer_on("SAPT(DFT):elst")
    elst = sapt_jk_terms.electrostatics(cache, True)
    data.update(elst)
    core.timer_off("SAPT(DFT):elst")

    # Exchange
    core.timer_on("SAPT(DFT):exch")
    exch = sapt_jk_terms.exchange(cache, sapt_jk, True)
    data.update(exch)
    core.timer_off("SAPT(DFT):exch")

    # Induction
    core.timer_on("SAPT(DFT):ind")
    ind = sapt_jk_terms.induction(
        cache,
        sapt_jk,
        True,
        sapt_jk_B=sapt_jk_B,
        maxiter=core.get_option("SAPT", "MAXITER"),
        conv=core.get_option("SAPT", "CPHF_R_CONVERGENCE"),
        Sinf=core.get_option("SAPT", "DO_IND_EXCH_SINF"),
    )
    data.update(ind)

    # Set Delta HF for SAPT(HF)
    if delta_hf:
        total_sapt = (
            data["Elst10,r"] + data["Exch10"] +
            data["Ind20,r"] + data["Exch-Ind20,r"]
        )
        sapt_hf_delta = data["DHF VALUE"] - total_sapt
        core.set_variable("SAPT(DFT) Delta HF", sapt_hf_delta)
        data["Delta HF Correction"] = core.variable("SAPT(DFT) Delta HF")

    core.timer_off("SAPT(DFT):ind")

    # Blow away JK object before doing MP2 for memory considerations
    if cleanup_jk:
        sapt_jk.finalize()

    # Hybrid xc kernel check
    do_hybrid = core.get_option("SAPT", "SAPT_DFT_DO_HYBRID")
    is_x_hybrid = wfn_B.functional().is_x_hybrid()
    is_x_lrc = wfn_B.functional().is_x_lrc()
    hybrid_specified = core.has_option_changed("SAPT", "SAPT_DFT_DO_HYBRID")
    if is_x_lrc:
        if do_hybrid:
            if hybrid_specified:
                raise ValidationError(
                    "SAPT(DFT): Hybrid xc kernel not yet implemented for range-separated funtionals."
                )
            else:
                core.print_out(
                    "Warning: Hybrid xc kernel not yet implemented for range-separated funtionals; hybrid kernel capability is turned off.\n"
                )
        is_hybrid = False
    else:
        if do_hybrid:
            is_hybrid = is_x_hybrid
        else:
            is_hybrid = False

    # Dispersion
    core.timer_on("SAPT(DFT):disp")

    primary_basis = wfn_A.basisset()
    aux_basis = core.BasisSet.build(
        dimer_wfn.molecule(),
        "DF_BASIS_MP2",
        core.get_option("DFMP2", "DF_BASIS_MP2"),
        "RIFIT",
        core.get_global_option("BASIS"),
    )

    if do_dft:
        core.timer_on("FDDS disp")
        core.print_out("\n")
        x_alpha = wfn_B.functional().x_alpha()
        if not is_hybrid:
            x_alpha = 0.0
        fdds_disp = sapt_mp2_terms.df_fdds_dispersion(
            primary_basis, aux_basis, cache, is_hybrid, x_alpha
        )
        data.update(fdds_disp)
        nfrozen_A = 0
        nfrozen_B = 0
        core.timer_off("FDDS disp")
    else:
        # this is where we actually need to figure out the number of frozen-core orbitals
        # if SAPT_DFT_MP2_DISP_ALG == FISAPT, the code will not figure it out on its own
        nfrozen_A = wfn_A.basisset().n_frozen_core(
            core.get_global_option("FREEZE_CORE"), wfn_A.molecule()
        )
        nfrozen_B = wfn_B.basisset().n_frozen_core(
            core.get_global_option("FREEZE_CORE"), wfn_B.molecule()
        )

    core.timer_on("MP2 disp")
    if core.get_option("SAPT", "SAPT_DFT_MP2_DISP_ALG") == "FISAPT":
        mp2_disp = sapt_mp2_terms.df_mp2_fisapt_dispersion(
            wfn_A, primary_basis, aux_basis, cache, nfrozen_A, nfrozen_B, do_print=True
        )
    else:
        mp2_disp = sapt_mp2_terms.df_mp2_sapt_dispersion(
            dimer_wfn, wfn_A, wfn_B, primary_basis, aux_basis, cache, do_print=True
        )
    data.update(mp2_disp)

    # Exchange-dispersion scaling
    if do_dft:
        exch_disp_scheme = core.get_option(
            "SAPT", "SAPT_DFT_EXCH_DISP_SCALE_SCHEME")
        core.print_out("    %-33s % s\n" %
                       ("Scaling Scheme", exch_disp_scheme))
        if exch_disp_scheme == "NONE":
            data["Exch-Disp20,r"] = data["Exch-Disp20,u"]
        elif exch_disp_scheme == "FIXED":
            exch_disp_scale = core.get_option(
                "SAPT", "SAPT_DFT_EXCH_DISP_FIXED_SCALE")
            core.print_out("    %-28s % 10.3f\n" %
                           ("Scaling Factor", exch_disp_scale))
            data["Exch-Disp20,r"] = exch_disp_scale * data["Exch-Disp20,u"]
        elif exch_disp_scheme == "DISP":
            exch_disp_scale = data["Disp20"] / data["Disp20,u"]
            data["Exch-Disp20,r"] = exch_disp_scale * data["Exch-Disp20,u"]
        if exch_disp_scheme != "NONE":
            core.print_out(
                print_sapt_var("Est. Exch-Disp20,r",
                               data["Exch-Disp20,r"], short=True)
                + "\n"
            )

    core.timer_off("MP2 disp")
    core.timer_off("SAPT(DFT):disp")

    # Print out final data
    core.print_out("\n")
    core.print_out(print_sapt_dft_summary(data, "SAPT(DFT)", do_dft=do_dft))

    return data


def run_sf_sapt(name, **kwargs):
    optstash = p4util.OptionsState(
        ["SCF_TYPE"],
        ["SCF", "REFERENCE"],
        ["SCF", "DFT_GRAC_SHIFT"],
        ["SCF", "SAVE_JK"],
    )

    core.tstart()

    # Alter default algorithm
    if not core.has_global_option_changed("SCF_TYPE"):
        core.set_global_option("SCF_TYPE", "DF")

    core.prepare_options_for_module("SAPT")

    # Get the molecule of interest
    ref_wfn = kwargs.get("ref_wfn", None)
    if ref_wfn is None:
        sapt_dimer = kwargs.pop("molecule", core.get_active_molecule())
    else:
        core.print_out(
            'Warning! SAPT argument "ref_wfn" is only able to use molecule information.'
        )
        sapt_dimer = ref_wfn.molecule()

    sapt_dimer, monomerA, monomerB = proc_util.prepare_sapt_molecule(
        sapt_dimer, "dimer"
    )

    # Print out the title and some information
    core.print_out("\n")
    core.print_out(
        "         ---------------------------------------------------------\n"
    )
    core.print_out("         " + "Spin-Flip SAPT Procedure".center(58) + "\n")
    core.print_out("\n")
    core.print_out(
        "         " +
        "by Daniel G. A. Smith and Konrad Patkowski".center(58) + "\n"
    )
    core.print_out(
        "         ---------------------------------------------------------\n"
    )
    core.print_out("\n")

    core.print_out("  ==> Algorithm <==\n\n")
    core.print_out(
        "   JK Algorithm            %12s\n" % core.get_option(
            "SCF", "SCF_TYPE")
    )
    core.print_out("\n")
    core.print_out("   Required computations:\n")
    core.print_out("     HF  (Monomer A)\n")
    core.print_out("     HF  (Monomer B)\n")
    core.print_out("\n")

    if core.get_option("SCF", "REFERENCE") != "ROHF":
        raise ValidationError(
            "Spin-Flip SAPT currently only supports restricted open-shell references."
        )

    # Run the two monomer computations
    core.IO.set_default_namespace('dimer')

    core.set_global_option('DF_INTS_IO', 'SAVE')

    # Compute dimer wavefunction
    wfn_A = scf_helper(
        "SCF", molecule=monomerA, banner="SF-SAPT: HF Monomer A", **kwargs
    )

    core.set_global_option("SAVE_JK", True)
    wfn_B = scf_helper(
        "SCF", molecule=monomerB, banner="SF-SAPT: HF Monomer B", **kwargs
    )
    sapt_jk = wfn_B.jk()
    core.set_global_option("SAVE_JK", False)
    core.print_out("\n")
    core.print_out(
        "         ---------------------------------------------------------\n"
    )
    core.print_out(
        "         " +
        "Spin-Flip SAPT Exchange and Electrostatics".center(58) + "\n"
    )
    core.print_out("\n")
    core.print_out(
        "         " +
        "by Daniel G. A. Smith and Konrad Patkowski".center(58) + "\n"
    )
    core.print_out(
        "         ---------------------------------------------------------\n"
    )
    core.print_out("\n")

    sf_data = sapt_sf_terms.compute_sapt_sf(sapt_dimer, sapt_jk, wfn_A, wfn_B)

    # Print the results
    core.print_out("   Spin-Flip SAPT Results\n")
    core.print_out("  " + "-" * 103 + "\n")

    for key, value in sf_data.items():
        value = sf_data[key]
        print_vals = (
            key,
            value * 1000,
            value * constants.hartree2kcalmol,
            value * constants.hartree2kJmol,
        )
        string = (
            "    %-26s % 15.8f [mEh] % 15.8f [kcal/mol] % 15.8f [kJ/mol]\n" % print_vals
        )
        core.print_out(string)
    core.print_out("  " + "-" * 103 + "\n\n")

    dimer_wfn = core.Wavefunction.build(sapt_dimer, wfn_A.basisset())

    # Set variables
    psivar_tanslator = {
        "Elst10": "SAPT ELST ENERGY",
        "Exch10(S^2) [diagonal]": "SAPT EXCH10(S^2),DIAGONAL ENERGY",
        "Exch10(S^2) [off-diagonal]": "SAPT EXCH10(S^2),OFF-DIAGONAL ENERGY",
        "Exch10(S^2) [highspin]": "SAPT EXCH10(S^2),HIGHSPIN ENERGY",
    }

    for k, v in sf_data.items():
        psi_k = psivar_tanslator[k]

        dimer_wfn.set_variable(psi_k, v)
        core.set_variable(psi_k, v)

    # Copy over highspin
    core.set_variable("SAPT EXCH ENERGY", sf_data["Exch10(S^2) [highspin]"])

    core.tstop()

    return dimer_wfn
