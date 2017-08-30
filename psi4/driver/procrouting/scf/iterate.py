"""
The SCF iteration functions
"""

import psi4
import numpy

@classmethod
def scf_iterate(self):

    reference = psi4.core.get_option("REFERENCE")

    MOM_performed_ = False
    diis_performed_ = False

    df = psi4.core.get_option("SCF_TYPE") == "DF"

    psi4.core.print_out( "  ==> Iterations <==\n\n");
    psi4.core.print_out( "%s                        Total Energy        Delta E     RMS |[F,P]|\n\n", df  "   " : "");


    // SCF iterations
    iteration = 0
    self.save_density_and_energy()

    Horig = self.H().clone()
    # self.H().copy(Horig)
    # self.H().axpy(1.0, Vefp)

# #ifdef USING_libefp
#         // add efp contribution to Fock matrix
#         if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
#             H_->copy(Horig_);
#             std::shared_ptr<Matrix> Vefp = Process::environment.get_efp()->modify_Fock_induced();
#             H_->add(Vefp);
#         }
# #endif

        E_ = 0.0;

        psi4.core.timer_on("HF: Form G");
        self.form_G();
        psi4.core.timer_off("HF: Form G");

        // Reset fractional SAD occupation
        if ((iteration_ == 0) && reset_occ_){
            reset_occupation();
        }

        timer_on("HF: Form F");
        form_F();
        timer_off("HF: Form F");

        if (print_>3) {
            Fa_->print("outfile");
            Fb_->print("outfile");
        }

        E_ += compute_E();

#ifdef USING_libefp
        // add efp contribution to energy
        if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
            double efp_wfn_dependent_energy = Process::environment.get_efp()->scf_energy_update();
            E_ += efp_wfn_dependent_energy;
        }
#endif

#ifdef USING_PCMSolver
        // The PCM potential must be added to the Fock operator *after* the
        // energy computation, not in form_F()
        if(pcm_enabled_) {
          // Prepare the density
          SharedMatrix D_pcm(Da_->clone());
          if(same_a_b_orbs()) {
            D_pcm->scale(2.0); // PSI4's density doesn't include the occupation
          } else {
            D_pcm->add(Db_);
          }

          // Compute the PCM charges and polarization energy
          double epcm = 0.0;
          if (options_.get_str("PCM_SCF_TYPE") == "TOTAL") {
            epcm = hf_pcm_->compute_E(D_pcm, PCM::Total);
          } else {
            epcm = hf_pcm_->compute_E(D_pcm, PCM::NucAndEle);
          }
          energies_["PCM Polarization"] = epcm;
          variables_["PCM POLARIZATION ENERGY"] = energies_["PCM Polarization"];
          Process::environment.globals["PCM POLARIZATION ENERGY"] = energies_["PCM Polarization"];
          E_ += epcm;

          // Add the PCM potential to the Fock matrix
          SharedMatrix V_pcm;
          V_pcm = hf_pcm_->compute_V();
          if (same_a_b_orbs()) {
            Fa_->add(V_pcm);
          } else {
            Fa_->add(V_pcm);
            Fb_->add(V_pcm);
          }
        } else {
          energies_["PCM Polarization"] = 0.0;
          variables_["PCM POLARIZATION ENERGY"] = energies_["PCM Polarization"];
          Process::environment.globals["PCM POLARIZATION ENERGY"] = energies_["PCM Polarization"];
        }
#endif
        std::string status = "";

        // We either do SOSCF or DIIS
        bool did_soscf = false;
        if (soscf_enabled_ && (Drms_ < soscf_r_start_) && (iteration_ > 3)) {
            compute_orbital_gradient(false);
            diis_performed_ = false;
            std::string base_name;
            if (functional_->needs_xc()) {
                base_name = "SOKS, nmicro = ";
            } else {
                base_name = "SOSCF, nmicro = ";
            }

            if (!test_convergency()) {
                int nmicro = soscf_update();
                if (nmicro > 0) {  // If zero the soscf call bounced for some reason
                    find_occupation();
                    status += base_name + psi::to_string(nmicro);
                    did_soscf = true;  // Stops DIIS
                } else {
                    if (print_) {
                        outfile->Printf(
                            "Did not take a SOSCF step, using normal convergence methods\n");
                    }
                    did_soscf = false;  // Back to DIIS
                }
            } else {
                // We need to ensure orthogonal orbitals and set epsilon
                status += base_name + "conv";
                timer_on("HF: Form C");
                form_C();
                timer_off("HF: Form C");
                did_soscf = true;  // Stops DIIS
            }
        }  // End SOSCF block

        if (!did_soscf){ // Normal convergence procedures if we do not do SOSCF

            timer_on("HF: DIIS");
            bool add_to_diis_subspace = false;
            if (diis_enabled_ && iteration_ > 0 && iteration_ >= diis_start_ )
                add_to_diis_subspace = true;

            compute_orbital_gradient(add_to_diis_subspace);

            if (diis_enabled_ == true && iteration_ >= diis_start_ + min_diis_vectors_ - 1) {
                diis_performed_ = diis();
            } else {
                diis_performed_ = false;
            }
            timer_off("HF: DIIS");

            if (print_>4 && diis_performed_) {
                outfile->Printf("  After DIIS:\n");
                Fa_->print("outfile");
                Fb_->print("outfile");
            }

            timer_on("HF: Form C");
            form_C();
            timer_off("HF: Form C");
        }

        // If we're too well converged, or damping wasn't enabled, do DIIS
        damping_performed_ = (damping_enabled_ && iteration_ > 1 && Drms_ > damping_convergence_);

        if(diis_performed_){
            if(status != "") status += "/";
            status += "DIIS";
        }
        if(MOM_performed_){
            if(status != "") status += "/";
            status += "MOM";
        }
        if(damping_performed_){
            if(status != "") status += "/";
            status += "DAMP";
        }
        if(frac_performed_){
            if(status != "") status += "/";
            status += "FRAC";
        }

        timer_on("HF: Form D");
        form_D();
        timer_off("HF: Form D");

        Process::environment.globals["SCF ITERATION ENERGY"] = E_;

        // After we've built the new D, damp the update if
        if(damping_performed_) damp_update();

        if (print_ > 3){
            Ca_->print("outfile");
            Cb_->print("outfile");
            Da_->print("outfile");
            Db_->print("outfile");
        }

        converged_ = test_convergency();

        df = (options_.get_str("SCF_TYPE") == "DF");


        outfile->Printf( "   @%s%s iter %3d: %20.14f   %12.5e   %-11.5e %s\n", df ? "DF-" : "",
                          reference.c_str(), iteration_, E_, E_ - Eold_, Drms_, status.c_str());


        // If a an excited MOM is requested but not started, don't stop yet
        if (MOM_excited_ && !MOM_started_) converged_ = false;

        // If a fractional occupation is requested but not started, don't stop yet
        if (frac_enabled_ && !frac_performed_) converged_ = false;

        // If a DF Guess environment, reset the JK object, and keep running
        if (converged_ && options_.get_bool("DF_SCF_GUESS") && (old_scf_type_ == "DIRECT")) {
            outfile->Printf( "\n  DF guess converged.\n\n"); // Be cool dude.
            converged_ = false;
            if(initialized_diis_manager_)
                diis_manager_->reset_subspace();
            scf_type_ = old_scf_type_;
            options_.set_str("SCF","SCF_TYPE",old_scf_type_);
            old_scf_type_ = "DF";
            integrals();
        }

        // Call any postiteration callbacks
//        call_postiteration_callbacks();

    } while (!converged_ && iteration_ < maxiter_ );

    return True



psi4.core.RHF.py_iterate = scf_iterate
psi4.core.UHF.py_iterate = scf_iterate
psi4.core.ROHF.py_iterate = scf_iterate
psi4.core.CUHF.py_iterate = scf_iterate






