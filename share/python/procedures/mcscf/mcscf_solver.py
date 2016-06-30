import psi4
import numpy as np
import diis_helper
from augmented_hessian import ah_iteration

def print_iteration(mtype, niter, energy, de, orb_rms, ci_rms, nci, norb, stype):
    psi4.print_out("%s %2d:  % 18.12f   % 1.4e  %1.2e  %1.2e  %3d  %3d  %s\n" % 
                    (mtype, niter, energy, de, orb_rms, ci_rms, nci, norb, stype))

def mcscf_solver(ref_wfn):

    # Build CIWavefunction
    psi4.prepare_options_for_module("DETCI")
    ciwfn = psi4.CIWavefunction(ref_wfn)
    ciwfn.set_print(0)

    # Build up some dimension information    
    docc_dim = ciwfn.get_dimension("DOCC")
    act_dim = ciwfn.get_dimension("ACT")
    vir_dim = ciwfn.get_dimension("VIR")
    oa_dim = docc_dim + act_dim
    av_dim = act_dim + vir_dim

    # Begin with a normal two-step    
    step_type = 'Initial CI'
    total_step = psi4.Matrix("Total step", oa_dim, av_dim)
    start_orbs = ciwfn.get_orbitals("ROT").clone()
    ciwfn.set_orbitals("ROT", start_orbs)

    # Start with SCF energy and other params 
    eold = psi4.get_variable("HF TOTAL ENERGY") 
    scf_energy = eold
    converged = False
    ediff = 1.e-4
    orb_grad_rms = 1.e-3
    one_step = False
    norb_iter = 1

    # Grab da options
    mcscf_orb_grad_conv = psi4.get_option("DETCI", "MCSCF_R_CONVERGENCE")
    mcscf_e_conv = psi4.get_option("DETCI", "MCSCF_E_CONVERGENCE")
    mcscf_max_macroiteration = psi4.get_option("DETCI", "MCSCF_MAXITER")
    mcscf_type = psi4.get_option("DETCI", 'MCSCF_TYPE')
    mcscf_steplimit = psi4.get_option("DETCI", 'MCSCF_MAX_ROT')
    mcscf_diis_start = psi4.get_option("DETCI", 'MCSCF_DIIS_START')
    mcscf_diis_freq = psi4.get_option("DETCI", 'MCSCF_DIIS_FREQ')
    mcscf_diis_max_vecs = psi4.get_option("DETCI", 'MCSCF_DIIS_MAX_VECS')
    mcscf_d_file = psi4.get_option("DETCI", 'CI_FILE_START') + 3
    mcscf_nroots = psi4.get_option("DETCI", 'NUM_ROOTS')

    # Grab needed objects 
    diis_obj = diis_helper.DIIS_helper(mcscf_diis_max_vecs)
    mcscf_obj = ciwfn.new_mcscf_object()

    if mcscf_type == "CONV":
        mtype = '   @MCSCF'
        psi4.print_out("\n   ==> Starting MCSCF iterations <==\n\n")
        psi4.print_out("        Iter         Total Energy       Delta E   Orb RMS    CI RMS  NCI NORB\n")
    else:
        mtype = '   @DF-MCSCF'
        psi4.print_out("\n   ==> Starting DF-MCSCF iterations <==\n\n") 
        psi4.print_out("           Iter         Total Energy       Delta E   Orb RMS    CI RMS  NCI NORB\n")

    # Iterate !
    for mcscf_iter in range(1, mcscf_max_macroiteration + 1):

        # Transform integrals, diagonalize H
        ciwfn.transform_mcscf_integrals(False)
        nci_iter = ciwfn.diag_h(abs(ediff) * 1.e-2, orb_grad_rms * 1.e-2)

        ciwfn.form_opdm()
        ciwfn.form_tpdm()
        ci_grad_norm = psi4.get_variable("DETCI AVG DVEC NORM")

        # Update MCSCF object    
        mcscf_obj.update(ciwfn.get_orbitals("DOCC"), ciwfn.get_orbitals("ACT"),
                         ciwfn.get_orbitals("VIR"), ciwfn.get_opdm(-1, -1, "SUM", False),
                         ciwfn.get_tpdm("SUM", True))
    
        CI_energy = psi4.get_variable('CURRENT ENERGY')
    
        orb_grad_rms = mcscf_obj.gradient_rms()
        ediff = CI_energy - eold

        # Print iterations
        print_iteration(mtype, mcscf_iter, CI_energy, ediff, orb_grad_rms, ci_grad_norm,
                        nci_iter, norb_iter, step_type) 

        eold = CI_energy

        # Check convergence
        if (orb_grad_rms < mcscf_orb_grad_conv) and (abs(ediff) < abs(mcscf_e_conv)) and\
            (mcscf_iter > 3):

            psi4.print_out("\n       %s has converged!\n\n" % mtype);
            converged = True
            break

        if orb_grad_rms > 1.e-3:
            step = mcscf_obj.approx_solve()
            step_type = 'TS'
        else:
            converged, norb_iter, step = ah_iteration(mcscf_obj, print_micro=False)
            norb_iter += 1     # Initial gradient
                
            step_type = 'AH'
            one_step = True

            if not converged:
                step = mcscf_obj.approx_solve()
                step_type = 'TS'
                one_step = False

    
        maxstep = step.absmax()
        if maxstep > mcscf_steplimit:
            psi4.print_out('    Warning! Maxstep = %4.2f, scaling to %4.2f\n' % (maxstep, mcscf_steplimit))
            step.scale(mcscf_steplimit / maxstep)
    
        total_step.add(step)

        # Do or add DIIS
        if (mcscf_iter >= mcscf_diis_start) and not one_step:
            #error = np.dot(orbs[:, :noa], np.array(grad)).dot(orbs[:, docc:].T)
            error = psi4.Matrix.triplet(ciwfn.get_orbitals("OA"),
                                        mcscf_obj.gradient(),
                                        ciwfn.get_orbitals("AV"),
                                        False, False, True)

            #print error
            diis_obj.add(total_step, error)

            if not (mcscf_iter % mcscf_diis_freq):
                total_step = diis_obj.extrapolate()
                step_type = 'TS, DIIS'
   
        # Finall rotate and set orbitals 
        orbs_mat = mcscf_obj.Ck(start_orbs, total_step)
        ciwfn.set_orbitals("ROT", orbs_mat)

    psi4.print_out(mtype + " Final Energy: %20.15f\n" % CI_energy)

    # Die if we did not converge
    if (not converged) and psi4.get_global_option("DIE_IF_NOT_CONVERGED"):
        raise Exception("MCSCF: Iterations did not converge!")

    # Print out energetics#
    psi4.print_out("\n   => Energetics <=\n\n")
    psi4.print_out("   SCF energy =         %20.15f\n" % scf_energy)
    psi4.print_out("   Total CI energy =    %20.15f\n\n" % CI_energy)

    dvec = ciwfn.new_civector(1, 53, True, True)
    dvec.set_nvec(mcscf_nroots)
    dvec.init_io_files(True)  
    for root in range(mcscf_nroots):
        root_e = psi4.get_variable("CI ROOT %d TOTAL ENERGY" % (root + 1))
        psi4.print_out("   CI Root %2d energy =  %20.15f\n" % (root + 1, root_e))
        ciwfn.print_vector(dvec, root)

 
    ciwfn.cleanup() 
    del diis_obj
    del mcscf_obj
    return ciwfn
    
