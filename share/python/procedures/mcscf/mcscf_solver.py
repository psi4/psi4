import os
import psi4
import numpy as np
import diis_helper
from augmented_hessian import ah_iteration
from quadratic_step import qc_iteration
import MCSCF

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
    qc_step = False

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
    mcscf_wavefunction_type = psi4.get_option('DETCI', 'WFN')
    mcscf_conv_type = 'OS'
    mcscf_ndet = ciwfn.ndet() 
    mcscf_nuclear_energy = ciwfn.molecule().nuclear_repulsion_energy()


    ########## DELETE ME
    docc = sum(ciwfn.get_dimension("DOCC").to_tuple())
    act = sum(ciwfn.get_dimension("ACT").to_tuple())
    vir = sum(ciwfn.get_dimension("VIR").to_tuple())

    noa = docc + act
    nav = act + vir
    mints = psi4.MintsHelper(ref_wfn.basisset())
    mints.integrals()
    Hcore = np.array(mints.ao_potential()) + np.array(mints.ao_kinetic())
    
    nbf = Hcore.shape[0]
    vir = Hcore.shape[0] - act - docc
    noa = docc + act
    nav = act + vir
    
    #I = np.array(mints.ao_eri())
    
    aux = psi4.BasisSet.pyconstruct_auxiliary(ciwfn.molecule(), 'DF_BASIS_SCF', psi4.get_global_option('DF_BASIS_MCSCF'),
                                              'JKFIT', psi4.get_global_option('BASIS'))
    zeroC = psi4.Matrix(nbf, nbf)
    dfobj = psi4.DFTensor(ciwfn.basisset(), aux, zeroC, noa, nbf - noa)
    
    Qpq = np.asarray(dfobj.Qso())
    nQ = Qpq.shape[0]
    
    
    Qpq = np.asarray(dfobj.Qso())
    I = np.einsum('Qpq,Qrs->pqrs', Qpq, Qpq)
    
    
    mcscf = MCSCF.MCSCF(docc, act, vir, Hcore, I, None, True)
    ########## DELETE ME


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

        nci_iter = ciwfn.diag_h(1.e-14, 1.e-14)
        #nci_iter = ciwfn.diag_h(abs(ediff) * 1.e-2, orb_grad_rms * 1.e-3)

        ciwfn.form_opdm()
        ciwfn.form_tpdm()
        ci_grad_rms = psi4.get_variable("DETCI AVG DVEC NORM")

        # Update MCSCF object    
        mcscf_obj.update(ciwfn.get_orbitals("DOCC"), ciwfn.get_orbitals("ACT"),
                         ciwfn.get_orbitals("VIR"), ciwfn.get_opdm(-1, -1, "SUM", False),
                         ciwfn.get_tpdm("SUM", True))
        mcscf.update(ciwfn.get_orbitals("ROT").np,
                     ciwfn.get_opdm(-1, -1, "SUM", False).np,
                     ciwfn.get_tpdm("SUM", True).np)
        current_energy = mcscf_obj.current_total_energy()
        current_energy += mcscf_nuclear_energy
    
        if abs(current_energy - psi4.get_variable('CURRENT ENERGY')) > 1.e-10:
            raise ValidationError("MCSCF: CI and Orbital energies are diverging!")
    
        orb_grad_rms = mcscf_obj.gradient_rms()
        ediff = current_energy - eold

        # Print iterations
        print_iteration(mtype, mcscf_iter, current_energy, ediff, orb_grad_rms, ci_grad_rms,
                        nci_iter, norb_iter, step_type) 

        #print np.allclose(mcscf_obj.gradient(), mcscf.gradient())
        #print np.allclose(mcscf_obj.H_approx_diag(), mcscf.diag_hessian())
        eold = current_energy

        if (mcscf_conv_type == 'OS') and (orb_grad_rms < 1.e-3) and (mcscf_iter >= 3):
            qc_step = True

        # Check convergence
        if (orb_grad_rms < mcscf_orb_grad_conv) and (abs(ediff) < abs(mcscf_e_conv)) and\
            (mcscf_iter > 3) and not qc_step:

            psi4.print_out("\n       %s has converged!\n\n" % mtype);
            converged = True
            break


        # Which orbital convergence are we doing?
        if (orb_grad_rms < 5.e-3) and mcscf_wavefunction_type == 'CASSCF' and False:
            converged, norb_iter, step = ah_iteration(mcscf_obj, print_micro=False)
            norb_iter += 1     # Initial gradient
                
            if converged:
                step_type = 'AH'
                one_step = True
            else:
                psi4.print_out("      !Warning. Augmented Hessian did not converge. Taking an approx step.\n")
                step = mcscf_obj.approx_solve()
                step_type = 'TS'
                one_step = False

        else:
            step = mcscf_obj.approx_solve()
            step_type = 'TS'

        maxstep = step.absmax()
        if maxstep > mcscf_steplimit:
            psi4.print_out('    Warning! Maxstep = %4.2f, scaling to %4.2f\n' % (maxstep, mcscf_steplimit))
            step.scale(mcscf_steplimit / maxstep)
    
        total_step.add(step)

        # Do or add DIIS
#        if (mcscf_iter >= mcscf_diis_start) and not one_step:
        if False:
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
        #orbs_mat.np[:] = mcscf.rotate_orbs(start_orbs.np, -total_step.np)
        ciwfn.set_orbitals("ROT", orbs_mat)

        # If we are taking a QC step do the following:
        if qc_step:
            break

    # If we converged do not do onestep
    if converged or (mcscf_conv_type != 'OS'):
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
        step_type = 'OS, Prep'
   
    # Loop for onstep 
    for mcscf_iter in one_step_iters:
#    for mcscf_iter in range(3):

        # Check convergence
        #if (orb_grad_rms < mcscf_orb_grad_conv) and (abs(ediff) < abs(mcscf_e_conv)) and\
        #    (mcscf_iter > 3):

        #    psi4.print_out("\n       %s has converged!\n\n" % mtype);
        #    converged = True
        #    break

        # Transform integrals and update the MCSCF object
        ciwfn.transform_mcscf_integrals(False)
        ciwfn.form_opdm()
        ciwfn.form_tpdm()

        mcscf_obj.update(ciwfn.get_orbitals("DOCC"), ciwfn.get_orbitals("ACT"),
                         ciwfn.get_orbitals("VIR"), ciwfn.get_opdm(-1, -1, "SUM", False),
                         ciwfn.get_tpdm("SUM", True))

        orb_grad_rms = mcscf_obj.gradient_rms()
        current_energy = mcscf_obj.current_total_energy()
        current_energy += mcscf_nuclear_energy

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
                        nci_iter, norb_iter, step_type) 
        step_type = 'OS'

        eold = current_energy

        if (orb_grad_rms < mcscf_orb_grad_conv) and (abs(ediff) < abs(mcscf_e_conv)):

            psi4.print_out("\n       %s has converged!\n\n" % mtype);
            converged = True
            break
        if mcscf_iter > 7:
            break

        # Take a step
        converged, norb_iter, nci_iter, step = qc_iteration(dvec, ci_grad, ciwfn, mcscf_obj)

        # Rotate integrals to new frame
        total_step.add(step)
        #orbs_mat = mcscf_obj.Ck(start_orbs, total_step)
        orbs_mat = mcscf_obj.Ck(ciwfn.get_orbitals("ROT"), step)
        ciwfn.set_orbitals("ROT", orbs_mat)


    psi4.print_out(mtype + " Final Energy: %20.15f\n" % current_energy)

    # Die if we did not converge
    if (not converged) and psi4.get_global_option("DIE_IF_NOT_CONVERGED"):
        raise Exception("MCSCF: Iterations did not converge!")

    # Print out energetics#
    psi4.print_out("\n   => Energetics <=\n\n")
    psi4.print_out("   SCF energy =         %20.15f\n" % scf_energy)
    psi4.print_out("   Total CI energy =    %20.15f\n\n" % current_energy)

    #dvec = ciwfn.new_civector(1, 53, True, True)
    #dvec.set_nvec(mcscf_nroots)
    #dvec.init_io_files(True)  
    #for root in range(mcscf_nroots):
    #    root_e = psi4.get_variable("CI ROOT %d TOTAL ENERGY" % (root + 1))
    #    psi4.print_out("   CI Root %2d energy =  %20.15f\n" % (root + 1, root_e))
    #    ciwfn.print_vector(dvec, root)

 
    ciwfn.cleanup() 
    del diis_obj
    del mcscf_obj
    return ciwfn
    
