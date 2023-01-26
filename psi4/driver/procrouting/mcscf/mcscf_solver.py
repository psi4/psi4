#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
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

import os

import numpy as np

from psi4 import core
from psi4.driver import p4util
from psi4.driver import qcdb
from psi4.driver.p4util import solvers
from .augmented_hessian import ah_iteration
from .. import proc_util


def print_iteration(mtype, niter, energy, de, orb_rms, ci_rms, nci, norb, stype):
    core.print_out("%s %2d:  % 18.12f   % 1.4e  %1.2e  %1.2e  %3d  %3d  %s\n" %
                    (mtype, niter, energy, de, orb_rms, ci_rms, nci, norb, stype))

def mcscf_solver(ref_wfn):

    # Build CIWavefunction
    core.prepare_options_for_module("DETCI")
    ciwfn = core.CIWavefunction(ref_wfn)
    ciwfn.set_module("detci")

    # Hush a lot of CI output
    ciwfn.set_print(0)

    # Begin with a normal two-step
    step_type = 'Initial CI'
    total_step = core.Matrix("Total step", ciwfn.get_dimension('OA'), ciwfn.get_dimension('AV'))
    start_orbs = ciwfn.get_orbitals("ROT").clone()
    ciwfn.set_orbitals("ROT", start_orbs)

    # Grab da options
    mcscf_orb_grad_conv = core.get_option("DETCI", "MCSCF_R_CONVERGENCE")
    mcscf_e_conv = core.get_option("DETCI", "MCSCF_E_CONVERGENCE")
    mcscf_max_macroiteration = core.get_option("DETCI", "MCSCF_MAXITER")
    mcscf_type = core.get_option("DETCI", "MCSCF_TYPE")
    mcscf_d_file = core.get_option("DETCI", "CI_FILE_START") + 3
    mcscf_nroots = core.get_option("DETCI", "NUM_ROOTS")
    mcscf_wavefunction_type = core.get_option("DETCI", "WFN")
    mcscf_ndet = ciwfn.ndet()
    mcscf_nuclear_energy = ciwfn.molecule().nuclear_repulsion_energy()
    mcscf_steplimit = core.get_option("DETCI", "MCSCF_MAX_ROT")
    mcscf_rotate = core.get_option("DETCI", "MCSCF_ROTATE")

    # DIIS info
    mcscf_diis_start = core.get_option("DETCI", "MCSCF_DIIS_START")
    mcscf_diis_freq = core.get_option("DETCI", "MCSCF_DIIS_FREQ")
    mcscf_diis_error_type = core.get_option("DETCI", "MCSCF_DIIS_ERROR_TYPE")
    mcscf_diis_max_vecs = core.get_option("DETCI", "MCSCF_DIIS_MAX_VECS")

    # One-step info
    mcscf_target_conv_type = core.get_option("DETCI", "MCSCF_ALGORITHM")
    mcscf_so_start_grad = core.get_option("DETCI", "MCSCF_SO_START_GRAD")
    mcscf_so_start_e = core.get_option("DETCI", "MCSCF_SO_START_E")
    mcscf_current_step_type = 'Initial CI'

    # Start with SCF energy and other params
    scf_energy = ciwfn.variable("HF TOTAL ENERGY")
    eold = scf_energy
    norb_iter = 1
    converged = False
    ah_step = False
    qc_step = False
    approx_integrals_only = True

    # Fake info to start with the initial diagonalization
    ediff = 1.e-4
    orb_grad_rms = 1.e-3

    # Grab needed objects
    diis_obj = solvers.DIIS(mcscf_diis_max_vecs)
    mcscf_obj = ciwfn.mcscf_object()

    # Execute the rotate command
    for rot in mcscf_rotate:
        if len(rot) != 4:
            raise p4util.PsiException("Each element of the MCSCF rotate command requires 4 arguements (irrep, orb1, orb2, theta).")

        irrep, orb1, orb2, theta = rot
        if irrep > ciwfn.Ca().nirrep():
            raise p4util.PsiException("MCSCF_ROTATE: Expression %s irrep number is larger than the number of irreps" %
                                    (str(rot)))

        if max(orb1, orb2) > ciwfn.Ca().coldim()[irrep]:
            raise p4util.PsiException("MCSCF_ROTATE: Expression %s orbital number exceeds number of orbitals in irrep" %
                                    (str(rot)))

        theta = np.deg2rad(theta)

        x = ciwfn.Ca().nph[irrep][:, orb1].copy()
        y = ciwfn.Ca().nph[irrep][:, orb2].copy()

        xp = np.cos(theta) * x - np.sin(theta) * y
        yp = np.sin(theta) * x + np.cos(theta) * y

        ciwfn.Ca().nph[irrep][:, orb1] = xp
        ciwfn.Ca().nph[irrep][:, orb2] = yp


    # Limited RAS functionality
    if core.get_local_option("DETCI", "WFN") == "RASSCF" and mcscf_target_conv_type != "TS":
        core.print_out("\n  Warning! Only the TS algorithm for RASSCF wavefunction is currently supported.\n")
        core.print_out("             Switching to the TS algorithm.\n\n")
        mcscf_target_conv_type = "TS"

    # Print out headers
    if mcscf_type == "CONV":
        mtype = "   @MCSCF"
        core.print_out("\n   ==> Starting MCSCF iterations <==\n\n")
        core.print_out("        Iter         Total Energy       Delta E   Orb RMS    CI RMS  NCI NORB\n")
    elif mcscf_type == "DF":
        mtype = "   @DF-MCSCF"
        core.print_out("\n   ==> Starting DF-MCSCF iterations <==\n\n")
        core.print_out("           Iter         Total Energy       Delta E   Orb RMS    CI RMS  NCI NORB\n")
    else:
        mtype = "   @AO-MCSCF"
        core.print_out("\n   ==> Starting AO-MCSCF iterations <==\n\n")
        core.print_out("           Iter         Total Energy       Delta E   Orb RMS    CI RMS  NCI NORB\n")

    # Iterate !
    for mcscf_iter in range(1, mcscf_max_macroiteration + 1):

        # Transform integrals, diagonalize H
        ciwfn.transform_mcscf_integrals(approx_integrals_only)
        nci_iter = ciwfn.diag_h(abs(ediff) * 1.e-2, orb_grad_rms * 1.e-3)

        # After the first diag we need to switch to READ
        ciwfn.set_ci_guess("DFILE")

        ciwfn.form_opdm()
        ciwfn.form_tpdm()
        ci_grad_rms = ciwfn.variable("DETCI AVG DVEC NORM")

        # Update MCSCF object
        Cocc = ciwfn.get_orbitals("DOCC")
        Cact = ciwfn.get_orbitals("ACT")
        Cvir = ciwfn.get_orbitals("VIR")
        opdm = ciwfn.get_opdm(-1, -1, "SUM", False)
        tpdm = ciwfn.get_tpdm("SUM", True)
        mcscf_obj.update(Cocc, Cact, Cvir, opdm, tpdm)

        current_energy = ciwfn.variable("MCSCF TOTAL ENERGY")

        ciwfn.reset_ci_H0block()

        orb_grad_rms = mcscf_obj.gradient_rms()
        ediff = current_energy - eold

        # Print iterations
        print_iteration(mtype, mcscf_iter, current_energy, ediff, orb_grad_rms, ci_grad_rms,
                        nci_iter, norb_iter, mcscf_current_step_type)
        eold = current_energy

        if mcscf_current_step_type == 'Initial CI':
            mcscf_current_step_type = 'TS'

        # Check convergence
        if (orb_grad_rms < mcscf_orb_grad_conv) and (abs(ediff) < abs(mcscf_e_conv)) and\
            (mcscf_iter > 3) and not qc_step:

            core.print_out("\n       %s has converged!\n\n" % mtype)
            converged = True
            break


        # Which orbital convergence are we doing?
        if ah_step:
            converged, norb_iter, step = ah_iteration(mcscf_obj, print_micro=False)
            norb_iter += 1

            if converged:
                mcscf_current_step_type = 'AH'
            else:
                core.print_out("      !Warning. Augmented Hessian did not converge. Taking an approx step.\n")
                step = mcscf_obj.approx_solve()
                mcscf_current_step_type = 'TS, AH failure'

        else:
            step = mcscf_obj.approx_solve()
            step_type = 'TS'

        maxstep = step.absmax()
        if maxstep > mcscf_steplimit:
            core.print_out('      Warning! Maxstep = %4.2f, scaling to %4.2f\n' % (maxstep, mcscf_steplimit))
            step.scale(mcscf_steplimit / maxstep)

        xstep = total_step.clone()
        total_step.add(step)

        # Do or add DIIS
        if (mcscf_iter >= mcscf_diis_start) and ("TS" in mcscf_current_step_type):

            # Figure out DIIS error vector
            if mcscf_diis_error_type == "GRAD":
                error = core.triplet(ciwfn.get_orbitals("OA"),
                                            mcscf_obj.gradient(),
                                            ciwfn.get_orbitals("AV"),
                                            False, False, True)
            else:
                error = step

            diis_obj.add(total_step, error)

            if not (mcscf_iter % mcscf_diis_freq):
                total_step = diis_obj.extrapolate()
                mcscf_current_step_type = 'TS, DIIS'

        # Build the rotation by continuous updates
        if mcscf_iter == 1:
            totalU = mcscf_obj.form_rotation_matrix(total_step)
        else:
            xstep.axpy(-1.0, total_step)
            xstep.scale(-1.0)
            Ustep = mcscf_obj.form_rotation_matrix(xstep)
            totalU = core.doublet(totalU, Ustep, False, False)

        # Build the rotation directly (not recommended)
        # orbs_mat = mcscf_obj.Ck(start_orbs, total_step)

        # Finally rotate and set orbitals
        orbs_mat = core.doublet(start_orbs, totalU, False, False)
        ciwfn.set_orbitals("ROT", orbs_mat)

        # Figure out what the next step should be
        if (orb_grad_rms < mcscf_so_start_grad) and (abs(ediff) < abs(mcscf_so_start_e)) and\
                (mcscf_iter >= 2):

            if mcscf_target_conv_type == 'AH':
                approx_integrals_only = False
                ah_step = True
            else:
                continue
        #raise p4util.PsiException("")

    core.print_out(mtype + " Final Energy: %20.15f\n" % current_energy)

    # Die if we did not converge
    if (not converged):
        if core.get_global_option("DIE_IF_NOT_CONVERGED"):
            raise p4util.PsiException("MCSCF: Iterations did not converge!")
        else:
            core.print_out("\nWarning! MCSCF iterations did not converge!\n\n")

    # For orbital invariant methods we transform the orbitals to the natural or
    # semicanonical basis. Frozen doubly occupied and virtual orbitals are not
    # modified.
    if core.get_option("DETCI", "WFN") == "CASSCF":
        # Do we diagonalize the opdm?
        if core.get_option("DETCI", "NAT_ORBS"):
            ciwfn.ci_nat_orbs()
        else:
            ciwfn.semicanonical_orbs()

        # Retransform intragrals and update CI coeffs., OPDM, and TPDM
        ciwfn.transform_mcscf_integrals(approx_integrals_only)
        ciwfn.set_print(1)
        ciwfn.set_ci_guess("H0_BLOCK")
        nci_iter = ciwfn.diag_h(mcscf_e_conv, mcscf_e_conv ** 0.5)

        ciwfn.form_opdm()
        ciwfn.form_tpdm()

    proc_util.print_ci_results(ciwfn, "MCSCF", scf_energy, current_energy, print_opdm_no=True)

    # Set final energy
    ciwfn.set_variable("CURRENT ENERGY", ciwfn.variable("MCSCF TOTAL ENERGY"))
    ciwfn.set_energy(ciwfn.variable("MCSCF TOTAL ENERGY"))

    # What do we need to cleanup?
    if core.get_option("DETCI", "MCSCF_CI_CLEANUP"):
        ciwfn.cleanup_ci()
    if core.get_option("DETCI", "MCSCF_DPD_CLEANUP"):
        ciwfn.cleanup_dpd()

    del diis_obj
    del mcscf_obj
    return ciwfn

