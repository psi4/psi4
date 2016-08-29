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


import os
import psi4
import numpy as np
import diis_helper
from augmented_hessian import ah_iteration
from quadratic_step import qc_iteration

# Relative hack for now
import sys, inspect
path_dir = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../../")))
sys.path.append(path_dir)
import p4util
from p4util.exceptions import *
np.set_printoptions(precision=5, linewidth=200, threshold=2000, suppress=True)
import scipy.linalg

def print_iteration(mtype, niter, energy, de, orb_rms, ci_rms, nci, norb, stype):
    psi4.print_out("%s %2d:  % 18.12f   % 1.4e  %1.2e  %1.2e  %3d  %3d  %s\n" %
                    (mtype, niter, energy, de, orb_rms, ci_rms, nci, norb, stype))

def mcscf_solver(ref_wfn):

    # Build CIWavefunction
    psi4.prepare_options_for_module("DETCI")
    ciwfn = psi4.CIWavefunction(ref_wfn)

    # Hush a lot of CI output
    ciwfn.set_print(0)

    # Begin with a normal two-step
    step_type = 'Initial CI'
    total_step = psi4.Matrix("Total step", ciwfn.get_dimension('OA'), ciwfn.get_dimension('AV'))
    start_orbs = ciwfn.get_orbitals("ROT").clone()
    ciwfn.set_orbitals("ROT", start_orbs)

    # Grab da options
    mcscf_orb_grad_conv = psi4.get_option("DETCI", "MCSCF_R_CONVERGENCE")
    mcscf_e_conv = psi4.get_option("DETCI", "MCSCF_E_CONVERGENCE")
    mcscf_max_macroiteration = psi4.get_option("DETCI", "MCSCF_MAXITER")
    mcscf_type = psi4.get_option("DETCI", "MCSCF_TYPE")
    mcscf_d_file = psi4.get_option("DETCI", "CI_FILE_START") + 3
    mcscf_nroots = psi4.get_option("DETCI", "NUM_ROOTS")
    mcscf_wavefunction_type = psi4.get_option("DETCI", "WFN")
    mcscf_ndet = ciwfn.ndet()
    mcscf_nuclear_energy = ciwfn.molecule().nuclear_repulsion_energy()
    mcscf_steplimit = psi4.get_option("DETCI", "MCSCF_MAX_ROT")
    mcscf_rotate = psi4.get_option("DETCI", "MCSCF_ROTATE")

    # DIIS info
    mcscf_diis_start = psi4.get_option("DETCI", "MCSCF_DIIS_START")
    mcscf_diis_freq = psi4.get_option("DETCI", "MCSCF_DIIS_FREQ")
    mcscf_diis_error_type = psi4.get_option("DETCI", "MCSCF_DIIS_ERROR_TYPE")
    mcscf_diis_max_vecs = psi4.get_option("DETCI", "MCSCF_DIIS_MAX_VECS")

    # One-step info
    mcscf_target_conv_type = psi4.get_option("DETCI", "MCSCF_ALGORITHM")
    mcscf_so_start_grad = psi4.get_option("DETCI", "MCSCF_SO_START_GRAD")
    mcscf_so_start_e = psi4.get_option("DETCI", "MCSCF_SO_START_E")
    mcscf_current_step_type = 'Initial CI'

    # Start with SCF energy and other params
    scf_energy = psi4.get_variable("HF TOTAL ENERGY")
    eold = scf_energy
    norb_iter = 1
    converged = False
    ah_step = False
    qc_step = False

    # Fake info to start with the inital diagonalization
    ediff = 1.e-4
    orb_grad_rms = 1.e-3

    # Grab needed objects
    diis_obj = diis_helper.DIIS_helper(mcscf_diis_max_vecs)
    mcscf_obj = ciwfn.mcscf_object()

    # Execute the rotate command
    for rot in mcscf_rotate:
        if len(rot) != 4:
            raise PsiException("Each element of the MCSCF rotate command requires 4 arguements (irrep, orb1, orb2, theta).")

        irrep, orb1, orb2, theta = rot
        if irrep > ciwfn.Ca().nirrep():
            raise PsiException("MCSCF_ROTATE: Expression %s irrep number is larger than the number of irreps" %
                                    (str(rot)))

        if max(orb1, orb2) > ciwfn.Ca().coldim()[irrep]:
            raise PsiException("MCSCF_ROTATE: Expression %s orbital number exceeds number of orbitals in irrep" %
                                    (str(rot)))

        theta = np.deg2rad(theta)

        x = ciwfn.Ca().nph[irrep][:, orb1].copy()
        y = ciwfn.Ca().nph[irrep][:, orb2].copy()

        xp = np.cos(theta) * x - np.sin(theta) * y
        yp = np.sin(theta) * x + np.cos(theta) * y
        
        ciwfn.Ca().nph[irrep][:, orb1] = xp
        ciwfn.Ca().nph[irrep][:, orb2] = yp


    # Limited RAS functionality
    if psi4.get_local_option("DETCI", "WFN") == "RASSCF" and mcscf_target_conv_type != "TS":
        psi4.print_out("\n  Warning! Only the TS algorithm for RASSCF wavefunction is currently supported.\n")
        psi4.print_out("             Switching to the TS algorithm.\n\n")
        mcscf_target_conv_type = "TS"

    # Print out headers
    if mcscf_type == "CONV":
        mtype = "   @MCSCF"
        psi4.print_out("\n   ==> Starting MCSCF iterations <==\n\n")
        psi4.print_out("        Iter         Total Energy       Delta E   Orb RMS    CI RMS  NCI NORB\n")
    else:
        mtype = "   @DF-MCSCF"
        psi4.print_out("\n   ==> Starting DF-MCSCF iterations <==\n\n")
        psi4.print_out("           Iter         Total Energy       Delta E   Orb RMS    CI RMS  NCI NORB\n")

    # Iterate !
    for mcscf_iter in range(1, mcscf_max_macroiteration + 1):

        # Transform integrals, diagonalize H
        ciwfn.transform_mcscf_integrals(mcscf_current_step_type == 'TS')
        nci_iter = ciwfn.diag_h(abs(ediff) * 1.e-2, orb_grad_rms * 1.e-3)

        ciwfn.form_opdm()
        ciwfn.form_tpdm()
        ci_grad_rms = psi4.get_variable("DETCI AVG DVEC NORM")

        # Update MCSCF object
        mcscf_obj.update(ciwfn.get_orbitals("DOCC"), ciwfn.get_orbitals("ACT"),
                         ciwfn.get_orbitals("VIR"), ciwfn.get_opdm(-1, -1, "SUM", False),
                         ciwfn.get_tpdm("SUM", True))
        current_energy = psi4.get_variable('CURRENT ENERGY')

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

            psi4.print_out("\n       %s has converged!\n\n" % mtype);
            converged = True
            break


        # Which orbital convergence are we doing?
        if ah_step:
            converged, norb_iter, step = ah_iteration(mcscf_obj, print_micro=False)
            norb_iter += 1

            if converged:
                mcscf_current_step_type = 'AH'
            else:
                psi4.print_out("      !Warning. Augmented Hessian did not converge. Taking an approx step.\n")
                step = mcscf_obj.approx_solve()
                mcscf_current_step_type = 'TS, AH failure'

        else:
            step = mcscf_obj.approx_solve()
            step_type = 'TS'

        maxstep = step.absmax()
        if maxstep > mcscf_steplimit:
            psi4.print_out('      Warning! Maxstep = %4.2f, scaling to %4.2f\n' % (maxstep, mcscf_steplimit))
            step.scale(mcscf_steplimit / maxstep)

        total_step.add(step)

        # Do or add DIIS
        if (mcscf_iter >= mcscf_diis_start) and ("TS" in mcscf_current_step_type):

            # Figure out DIIS error vector
            if mcscf_diis_error_type == "GRAD":
                error = psi4.Matrix.triplet(ciwfn.get_orbitals("OA"),
                                            mcscf_obj.gradient(),
                                            ciwfn.get_orbitals("AV"),
                                            False, False, True)
            else:
                error = step

            diis_obj.add(total_step, error)

            if not (mcscf_iter % mcscf_diis_freq):
                total_step = diis_obj.extrapolate()
                mcscf_current_step_type = 'TS, DIIS'

        # Finally rotate and set orbitals
        orbs_mat = mcscf_obj.Ck(start_orbs, total_step)
        ciwfn.set_orbitals("ROT", orbs_mat)

        # Figure out what the next step should be
        if (orb_grad_rms < mcscf_so_start_grad) and (abs(ediff) < abs(mcscf_so_start_e)) and\
                (mcscf_iter >= 2):

            if mcscf_target_conv_type == 'AH':
                ah_step = True
            elif mcscf_target_conv_type == 'OS':
                mcscf_current_step_type = 'OS, Prep'
                break
            else:
                continue
        #raise PsiException("")

    # If we converged do not do onestep
    if converged or (mcscf_target_conv_type != 'OS'):
        one_step_iters = []

    # If we are not converged load in Dvec and build iters array
    else:
        one_step_iters = range(mcscf_iter + 1, mcscf_max_macroiteration + 1)
        dvec = ciwfn.new_civector(1, mcscf_d_file, True, True)
        dvec.set_nvec(1)
        dvec.init_io_files(True)
        dvec.read(0, 0)
        dvec.symnormalize(1.0, 0)

        ci_grad = ciwfn.new_civector(1, mcscf_d_file + 1, True, True)
        ci_grad.set_nvec(1)
        ci_grad.init_io_files(True)

    # Loop for onestep
    for mcscf_iter in one_step_iters:

        # Transform integrals and update the MCSCF object
        ciwfn.transform_mcscf_integrals(False)
        ciwfn.form_opdm()
        ciwfn.form_tpdm()

        mcscf_obj.update(ciwfn.get_orbitals("DOCC"), ciwfn.get_orbitals("ACT"),
                         ciwfn.get_orbitals("VIR"), ciwfn.get_opdm(-1, -1, "SUM", False),
                         ciwfn.get_tpdm("SUM", True))

        orb_grad_rms = mcscf_obj.gradient_rms()

        # Warning! Does not work for SA-MCSCF
        current_energy = mcscf_obj.current_total_energy()
        current_energy += mcscf_nuclear_energy

        psi4.set_variable("CI ROOT %d TOTAL ENERGY" % 1, current_energy)
        psi4.set_variable("CURRENT ENERGY", current_energy)

        docc_energy = mcscf_obj.current_docc_energy()
        ci_energy = mcscf_obj.current_ci_energy()

        # Compute CI gradient
        ciwfn.sigma(dvec, ci_grad, 0, 0)
        ci_grad.scale(2.0, 0)
        ci_grad.axpy(-2.0 * ci_energy, dvec, 0, 0)

        ci_grad_rms = ci_grad.norm(0)
        orb_grad_rms = mcscf_obj.gradient().rms()

        ediff = current_energy - eold

        print_iteration(mtype, mcscf_iter, current_energy, ediff, orb_grad_rms, ci_grad_rms,
                        nci_iter, norb_iter, mcscf_current_step_type)
        mcscf_current_step_type = 'OS'

        eold = current_energy

        if (orb_grad_rms < mcscf_orb_grad_conv) and (abs(ediff) < abs(mcscf_e_conv)):

            psi4.print_out("\n       %s has converged!\n\n" % mtype);
            converged = True
            break

        # Take a step
        converged, norb_iter, nci_iter, step = qc_iteration(dvec, ci_grad, ciwfn, mcscf_obj)

        # Rotate integrals to new frame
        total_step.add(step)
        orbs_mat = mcscf_obj.Ck(ciwfn.get_orbitals("ROT"), step)
        ciwfn.set_orbitals("ROT", orbs_mat)


    psi4.print_out(mtype + " Final Energy: %20.15f\n" % current_energy)

    # Die if we did not converge
    if (not converged):
        if psi4.get_global_option("DIE_IF_NOT_CONVERGED"):
            raise PsiException("MCSCF: Iterations did not converge!")
        else:
            psi4.print_out("\nWarning! MCSCF iterations did not converge!\n\n")

    # Print out energetics
    psi4.print_out("\n   ==> Energetics <==\n\n")
    psi4.print_out("    SCF energy =         %20.15f\n" % scf_energy)
    psi4.print_out("    Total CI energy =    %20.15f\n\n" % current_energy)

    # Print out CI vector information
    if mcscf_target_conv_type != 'SO':
        dvec = ciwfn.new_civector(mcscf_nroots, mcscf_d_file, True, True)
        dvec.init_io_files(True)

    for root in range(mcscf_nroots):
        psi4.print_out("\n   ==> CI root %d information <==\n\n" % (root + 1))

        # Print total energy
        root_e = psi4.get_variable("CI ROOT %d TOTAL ENERGY" % (root + 1))
        psi4.print_out("    CI Root %2d energy =  %20.15f\n" % (root + 1, root_e))

        # Print natural occupations
        psi4.print_out("\n   Natural occupation numbers:\n\n")
        ropdm = ciwfn.get_opdm(root, root, "SUM", False).to_array(dense=True)
        nocc, rot = np.linalg.eigh(ropdm)
        nocc = nocc[::-1]

        irrep_info = ['??' for x in range(nocc.shape[0])]

        cnt = 0
        for ln in range(nocc.shape[0] // 3 + 1):
            for sp in range(3):
                if cnt >= nocc.shape[0]: break

                psi4.print_out("      %4s  % 8.6f" % (irrep_info[cnt], nocc[cnt]))
                cnt += 1
            psi4.print_out("\n")

        # Print CIVector information
        ciwfn.print_vector(dvec, root)


    if psi4.get_option("DETCI", "MCSCF_CLEANUP"):
        ciwfn.cleanup()

    del diis_obj
    del mcscf_obj
    return ciwfn

