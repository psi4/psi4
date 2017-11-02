"""
The SCF iteration functions
"""

import psi4
import numpy

from psi4 import core


def scf_compute_energy(self):
    self.initialize()
    self.py_iterate()
    scf_energy = self.finalize_energy()
    return scf_energy


def scf_iterate(self):

    reference = core.get_option('SCF', "REFERENCE")

    MOM_performed = False
    diis_performed = False
    frac_performed = False  # revisit
    reset_occ = True if (core.get_option('SCF', 'GUESS') == 'SAD') else False
    # todo reset_occ was ripped out of nice logic so revisit
    print_lvl = core.get_option('SCF', "PRINT")

    # First, did the user request a different number of diis vectors?
    min_diis_vectors = core.get_option('SCF', "DIIS_MIN_VECS")
    max_diis_vectors = core.get_option('SCF', "DIIS_MAX_VECS")
    diis_start = core.get_option('SCF', "DIIS_START")
    diis_enabled = core.get_option('SCF', "DIIS")

    # Thresholds
    energy_threshold = core.get_option("SCF", "E_CONVERGENCE");
    density_threshold = core.get_option("SCF", "D_CONVERGENCE");

    # Don't perform DIIS if less than 2 vectors requested, or user requested a negative number
    if min_diis_vectors < 2:
        diis_enabled = False

    # damping?
    damping_enabled = core.has_option_changed('SCF', 'DAMPING_PERCENTAGE')
    if damping_enabled:
        damping_percentage = core.get_option('SCF', "DAMPING_PERCENTAGE") / 100.0
        if damping_percentage < 0.0 or damping_percentage > 1.0:
            raise ValidationError("DAMPING_PERCENTAGE must be between 0 and 100.")
        damping_convergence = core.get_option('SCF', "DAMPING_CONVERGENCE")

    df = core.get_option('SCF', "SCF_TYPE") == "DF"
    is_dfjk = core.get_global_option('SCF_TYPE').endswith('DF')

    core.print_out("  ==> Iterations <==\n\n")
    core.print_out("%s                        Total Energy        Delta E     RMS |[F,P]|\n\n" % ("   " if is_dfjk else ""))

    # Iterate !
    # SCF iterations
    SCFE_old = 0.0
    SCFE = 0.0
    Drms = 0.0
    for iteration in range(core.get_option('SCF', 'MAXITER')):

        self.save_density_and_energy()

        # #ifdef USING_libefp
        #Horig = self.H().clone()
        # self.H().copy(Horig)
        # self.H().axpy(1.0, Vefp)

        #         # add efp contribution to Fock matrix
        #         if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
        #             H_->copy(Horig_)
        #             std::shared_ptr<Matrix> Vefp = Process::environment.get_efp()->modify_Fock_induced()
        #             H_->add(Vefp)
        #         }
        # #endif

        SCFE = 0.0

        core.timer_on("HF: Form G")
        self.form_G()
        core.timer_off("HF: Form G")

        # Reset fractional SAD occupation
        if (iteration == 0) and reset_occ:
            self.reset_occupation()

        core.timer_on("HF: Form F")
        self.form_F()
        core.timer_off("HF: Form F")

        if print_lvl > 3:
            self.Fa().print_out()
            self.Fb().print_out()

        SCFE += self.compute_E()

        ##ifdef USING_libefp
        #        # add efp contribution to energy
        #        if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
        #            double efp_wfn_dependent_energy = Process::environment.get_efp()->scf_energy_update()
        #            E_ += efp_wfn_dependent_energy
        #        }
        ##endif
        #
        ##ifdef USING_PCMSolver
        #        # The PCM potential must be added to the Fock operator *after* the
        #        # energy computation, not in form_F()
        #        if(pcm_enabled_) {
        #          # Prepare the density
        #          SharedMatrix D_pcm(Da_->clone())
        #          if(same_a_b_orbs()) {
        #            D_pcm->scale(2.0) # PSI4's density doesn't include the occupation
        #          } else {
        #            D_pcm->add(Db_)
        #          }
        #
        #          # Compute the PCM charges and polarization energy
        #          double epcm = 0.0
        #          if (options_.get_str("PCM_SCF_TYPE") == "TOTAL") {
        #            epcm = hf_pcm_->compute_E(D_pcm, PCM::Total)
        #          } else {
        #            epcm = hf_pcm_->compute_E(D_pcm, PCM::NucAndEle)
        #          }
        #          energies_["PCM Polarization"] = epcm
        #          variables_["PCM POLARIZATION ENERGY"] = energies_["PCM Polarization"]
        #          Process::environment.globals["PCM POLARIZATION ENERGY"] = energies_["PCM Polarization"]
        #          E_ += epcm
        #
        #          # Add the PCM potential to the Fock matrix
        #          SharedMatrix V_pcm
        #          V_pcm = hf_pcm_->compute_V()
        #          if (same_a_b_orbs()) {
        #            Fa_->add(V_pcm)
        #          } else {
        #            Fa_->add(V_pcm)
        #            Fb_->add(V_pcm)
        #          }
        #        } else {
        #          energies_["PCM Polarization"] = 0.0
        #          variables_["PCM POLARIZATION ENERGY"] = energies_["PCM Polarization"]
        #          Process::environment.globals["PCM POLARIZATION ENERGY"] = energies_["PCM Polarization"]
        #        }
        ##endif
        # std::string status = ""
        self.set_energies("Total Energy", SCFE)
        Ediff = SCFE - SCFE_old
        SCFE_old = SCFE

        status = ""

        # We either do SOSCF or DIIS
        did_soscf = False
        soscf_enabled = core.get_option('SCF', 'SOSCF')
        soscf_r_start = core.get_option('SCF', 'SOSCF_START_CONVERGENCE')
        if (soscf_enabled and (Drms < soscf_r_start) and (iteration > 3)):

            Drms = self.compute_orbital_gradient(False)
            diis_performed = False
            base_name = ""
            if self.functional().needs_xc():
                base_name = "SOKS, nmicro = "
            else:
                base_name = "SOSCF, nmicro = "

            if (abs(Ediff) > energy_threshold) or (Drms > density_threshold):
                nmicro = self.soscf_update()
                if nmicro > 0:
                    # If zero the soscf call bounced for some reason
                    self.find_occupation()
                    status += base_name + str(nmicro)
                    did_soscf = True  # Stops DIIS
                else:
                    if print_lvl:
                        core.print_out("Did not take a SOSCF step, using normal convergence methods\n")
                    did_soscf = False  # Back to DIIS

            else:
                # We need to ensure orthogonal orbitals and set epsilon
                status += base_name + "conv"
                core.timer_on("HF: Form C")
                self.form_C()
                core.timer_off("HF: Form C")
                did_soscf = True  # Stops DIIS

        if not did_soscf:
            # Normal convergence procedures if we do not do SOSCF

            core.timer_on("HF: DIIS")
            add_to_diis_subspace = False
            if (diis_enabled and (iteration > 0) and (iteration >= diis_start)):
                add_to_diis_subspace = True

            Drms = self.compute_orbital_gradient(add_to_diis_subspace)

            if (diis_enabled and (iteration >= diis_start + min_diis_vectors - 1)):
                diis_performed = self.diis()
            else:
                diis_performed = False

            core.timer_off("HF: DIIS")

            if (print_lvl > 4) and diis_performed:
                core.print_out("  After DIIS:\n")
                self.Fa().print_out()
                self.Fb().print_out()

            core.timer_on("HF: Form C")
            self.form_C()
            core.timer_off("HF: Form C")

        # If we're too well converged, or damping wasn't enabled, do DIIS
        damping_performed = (damping_enabled and (iteration > 1) and (Drms > damping_convergence))

        if diis_performed:
            if (status != ""):
                status += "/"
            status += "DIIS"

        if MOM_performed:
            if (status != ""):
                status += "/"
            status += "MOM"

        if damping_performed:
            if (status != ""):
                status += "/"
            status += "DAMP"

        if frac_performed:
            if (status != ""):
                status += "/"
            status += "FRAC"

        core.timer_on("HF: Form D")
        self.form_D()
        core.timer_off("HF: Form D")

        core.set_variable("SCF ITERATION ENERGY", SCFE)

        # After we've built the new D, damp the update if
        if damping_performed:
            self.damp_update()

        if print_lvl > 3:
            self.Ca().print_out()
            self.Cb().print_out()
            self.Da().print_out()
            self.Db().print_out()

        converged = (abs(Ediff) < energy_threshold) and (Drms < density_threshold)

        # Print out the iteration
        df = core.get_option('SCF', "SCF_TYPE") == "DF"
        core.print_out("   @%s%s iter %3d: %20.14f   %12.5e   %-11.5e %s\n" %
                       ("DF-" if is_dfjk else "", reference, iteration + 1, SCFE, Ediff, Drms,
                        status))

        # If a an excited MOM is requested but not started, don't stop yet
        #if (MOM_excited_ && !MOM_started_) converged_ = false
        # revisit TODO

        # If a fractional occupation is requested but not started, don't stop yet
        #if (frac_enabled_ && !frac_performed_) converged_ = false

        ## If a DF Guess environment, reset the JK object, and keep running
        #if (converged_ && options_.get_bool("DF_SCF_GUESS") && (old_scf_type_ == "DIRECT")) {
        #    core.print_out( "\n  DF guess converged.\n\n") # Be cool dude.
        #    converged_ = false
        #    if(initialized_diis_manager_)
        #        diis_manager_->reset_subspace()
        #    scf_type_ = old_scf_type_
        #    options_.set_global_str("SCF_TYPE",old_scf_type_)
        #    old_scf_type_ = "DF"
        #    integrals()
        #    is_dfjk = False
        #}

        # Call any postiteration callbacks
        if converged:
            break

    return True


def finalize_scf_energy(self):
    # Perform wavefunction stability analysis before doing
    # anything on a wavefunction that may not be truly converged.
    # if core.get_option("STABILITY_ANALYSIS") != None:
    #     # We need the integral file, make sure it is written and
    #     # compute it if needed

    #     if core.get_option("STABILITY_ANALYSIS") != "UHF":

    # if(options_.get_str("STABILITY_ANALYSIS") != "NONE") {
    #     if(options_.get_str("REFERENCE") != "UHF") {
    #         psio_->open(PSIF_SO_TEI, PSIO_OPEN_OLD)
    #         if (psio_->tocscan(PSIF_SO_TEI, IWL_KEY_BUF) == NULL) {
    #             psio_->close(PSIF_SO_TEI,1)
    #             core.print_out("    SO Integrals not on disk, computing...")
    #             std::shared_ptr<MintsHelper> mints(new MintsHelper(basisset_, options_, 0))
    #             mints->integrals()
    #             core.print_out("done.\n")
    #         } else {
    #             psio_->close(PSIF_SO_TEI,1)
    #         }

    #     }
    #     bool follow = stability_analysis()

    #     while ( follow && !(attempt_number_ > max_attempts_) ) {

    #       attempt_number_++
    #       core.print_out("    Running SCF again with the rotated orbitals.\n")

    #       if(initialized_diis_manager_) diis_manager_->reset_subspace()
    #       # Reading the rotated orbitals in before starting iterations
    #       form_D()
    #       E_ = compute_initial_E()
    #       iterations()
    #       follow = stability_analysis()
    #     }
    #     if ( follow && (attempt_number_ > max_attempts_) ) {
    #       core.print_out( "    There's still a negative eigenvalue. Try modifying FOLLOW_STEP_SCALE\n")
    #       core.print_out("    or increasing MAX_ATTEMPTS (not available for PK integrals).\n")
    #     }
    # }

    # # At this point, we are not doing any more SCF cycles
    # # and we can compute and print final quantities.
    # #ifdef USING_libefp
    #     if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
    #         Process::environment.get_efp()->compute()

    #         double efp_wfn_independent_energy = Process::environment.globals["EFP TOTAL ENERGY"] -
    #                                             Process::environment.globals["EFP IND ENERGY"]
    #         energies_["EFP"] = Process::environment.globals["EFP TOTAL ENERGY"]

    #         core.print_out("    EFP excluding EFP Induction   %20.12f [Eh]\n", efp_wfn_independent_energy)
    #         core.print_out("    SCF including EFP Induction   %20.12f [Eh]\n", E_)

    #         E_ += efp_wfn_independent_energy

    #         core.print_out("    Total SCF including Total EFP %20.12f [Eh]\n", E_)
    #     }
    # #endif

    core.print_out("\n  ==> Post-Iterations <==\n\n")

    self.check_phases()
    self.compute_spin_contamination()
    self.frac_renormalize()
    reference = core.get_option("SCF", "REFERENCE")
    converged = True
    print_lvl = 0

    energy = self.get_energies("Total Energy")
    fail_on_maxiter = core.get_option("SCF", "FAIL_ON_MAXITER")
    if converged or not fail_on_maxiter:

        # Print the orbitals
        if print_lvl:
            self.print_orbitals()

        if converged:
            core.print_out("  Energy converged.\n\n")
        else:
            core.print_out("  Energy did not converge, but proceeding anyway.\n\n")

        prefix = "DF-" if is_dfjk else ''

        core.print_out("  @%s%s Final Energy: %20.14f" % (prefix, reference, energy))
        # if (perturb_h_) {
        #     core.print_out( " with %f %f %f perturbation" % (perturb_dipoles_[0], perturb_dipoles_[1], perturb_dipoles_[2]))
        # }
        core.print_out("\n\n")
        self.print_energies()

        # Need to recompute the Fock matrices, as they are modified during the SCF iteration
        # and might need to be dumped to checkpoint later
        self.form_F()
        #ifdef USING_PCMSolver
        # if(pcm_enabled_) {
        #     # Prepare the density
        #     SharedMatrix D_pcm(Da_->clone())
        #     if(same_a_b_orbs()) {
        #       D_pcm->scale(2.0) # PSI4's density doesn't include the occupation
        #     }
        #     else {
        #       D_pcm->add(Db_)
        #     }

        #     # Add the PCM potential to the Fock matrix
        #     SharedMatrix V_pcm
        #     V_pcm = hf_pcm_->compute_V()
        #     if(same_a_b_orbs()) Fa_->add(V_pcm)
        #     else {
        #       Fa_->add(V_pcm)
        #       Fb_->add(V_pcm)
        #     }
        # }
        #endif

        # Properties
        #  Comments so that autodoc utility will find these PSI variables
        #  Process::environment.globals["SCF DIPOLE X"] =
        #  Process::environment.globals["SCF DIPOLE Y"] =
        #  Process::environment.globals["SCF DIPOLE Z"] =
        #  Process::environment.globals["SCF QUADRUPOLE XX"] =
        #  Process::environment.globals["SCF QUADRUPOLE XY"] =
        #  Process::environment.globals["SCF QUADRUPOLE XZ"] =
        #  Process::environment.globals["SCF QUADRUPOLE YY"] =
        #  Process::environment.globals["SCF QUADRUPOLE YZ"] =
        #  Process::environment.globals["SCF QUADRUPOLE ZZ"] =

    else:
        core.print_out("  Failed to converge.\n")
        energy = 0.0

        # Throw if we didn't converge?
        self.die_if_not_converged()

    # Orbitals are always saved, in case an MO guess is requested later
    # save_orbitals()

    self.finalize()

    core.print_out("\nComputation Completed\n")

    return energy


core.HF.py_iterate = scf_iterate
core.HF.compute_energy = scf_compute_energy
core.HF.finalize_energy = finalize_scf_energy
# core.RHF.py_iterate = scf_iterate
# core.UHF.py_iterate = scf_iterate
# core.ROHF.py_iterate = scf_iterate
# core.CUHF.py_iterate = scf_iterate
