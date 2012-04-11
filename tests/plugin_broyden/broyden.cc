#include "broyden.h" 
#include <psifiles.h>
#include <libmints/view.h>
#include <libmints/mints.h>
#include <libfock/apps.h>
#include <libfock/jk.h>
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

#include <sstream>

using namespace psi;

namespace psi{ namespace scf{

BroydenRHF::BroydenRHF(Options& options) 
    : RHF(options, _default_psio_lib_)
{
    broyden_iteration_ = 0;
    broyden_status_ = "ROOTHAAN";
}
BroydenRHF::~BroydenRHF()
{
}

double BroydenRHF::compute_energy()
{
    std::string reference = options_.get_str("REFERENCE");

    bool converged = false;
    
    if (options_.get_str("GUESS") == "SAD" || options_.get_str("GUESS") == "READ")
        iteration_ = -1;
    else
        iteration_ = 0;

    if (print_ && (Communicator::world->me() == 0))
        fprintf(outfile, "  ==> Pre-Iterations <==\n\n");

    timer_on("Form H");
    form_H(); //Core Hamiltonian
    timer_off("Form H");

    timer_on("Form S/X");
    form_Shalf(); //S and X Matrix
    timer_off("Form S/X");

    timer_on("Guess");
    guess(); // Guess
    timer_off("Guess");

    if (print_)
        print_preiterations();

    integrals();

    if (Communicator::world->me() == 0) {
        fprintf(outfile, "  ==> Iterations <==\n\n");
        fprintf(outfile, "                        Total Energy        Delta E     Density RMS\n\n");
    }
    fflush(outfile);

    // SCF iterations
    do {
        iteration_++;

        save_density_and_energy();

        // Call any preiteration callbacks
        call_preiteration_callbacks();

        timer_on("Form G");
        form_G();
        timer_off("Form G");

        // Reset fractional SAD occupation
        if (iteration_ == 0 && options_.get_str("GUESS") == "SAD")
            reset_SAD_occupation();

        timer_on("Form F");
        form_F();
        timer_off("Form F");

        if (print_>3) {
            Fa_->print(outfile);
            Fb_->print(outfile);
        }

        E_ = compute_E();

        if (Communicator::world->me() == 0) {
            fprintf(outfile, "   @%s iter %3d: %20.14f   %12.5e   %-11.5e %s\n",
                              reference.c_str(), iteration_, E_, E_ - Eold_, Drms_, broyden_status_.c_str());
            fflush(outfile);
        }

        Process::environment.globals["SCF ITERATION ENERGY"] = E_;

        if (options_.get_str("SCF_STEP") == "BROYDEN" && iteration_ >= options_.get_int("BROYDEN_START")) {
            // Broyden
            broyden_step(); 
        } else if (options_.get_str("SCF_STEP") == "BFGS" && iteration_ >= options_.get_int("BROYDEN_START")) {
            // BFGS
            bfgs_step();
        } else if (options_.get_str("SCF_STEP") == "DFP" && iteration_ >= options_.get_int("BROYDEN_START")) {
            // DFP 
            dfp_step();
        } else if (options_.get_str("SCF_STEP") == "SR1" && iteration_ >= options_.get_int("BROYDEN_START")) {
            // SR1 
            sr1_step();
        } else {
            // Roothaan
            timer_on("Form C");
            form_C(); // Calls find_occupation
            timer_off("Form C");
        }

        timer_on("Form D");
        form_D();
        timer_off("Form D");

        if (print_ > 3){
            Ca_->print(outfile);
            Cb_->print(outfile);
            Da_->print(outfile);
            Db_->print(outfile);
        }

        converged = test_convergency();

        // Call any postiteration callbacks
        call_postiteration_callbacks();

    } while (!converged && iteration_ < maxiter_ );

    if (Communicator::world->me() == 0)
        fprintf(outfile, "\n  ==> Post-Iterations <==\n\n");

    compute_spin_contamination();

    if (converged) {
        // Need to recompute the Fock matrices, as they are modified during the SCF interation
        // and might need to be dumped to checkpoint later
        form_F();

        // Orbitals need to be updated to canonicalize
        form_C();    
    
        // Print the orbitals
        if(print_)
            print_orbitals();

        if (Communicator::world->me() == 0) {
            fprintf(outfile, "\n  Energy converged.\n");
            fprintf(outfile, "\n  @%s Final Energy: %20.14f",reference.c_str(), E_);
            if (perturb_h_) {
                fprintf(outfile, " with %f perturbation", lambda_);
            }
            fprintf(outfile, "\n");
        }

        // Properties
        //if (print_) {
        //    boost::shared_ptr<OEProp> oe(new OEProp());
        //    oe->set_title("SCF");
        //    oe->add("DIPOLE");

        //    if (print_ >= 2) {
        //        oe->add("QUADRUPOLE");
        //        oe->add("MULLIKEN_CHARGES");
        //    }

        //    if (print_ >= 3) {
        //        oe->add("LOWDIN_CHARGES");
        //        oe->add("MAYER_INDICES");
        //        oe->add("WIBERG_LOWDIN_INDICES");
        //    }

        //    if (Communicator::world->me() == 0)
        //        fprintf(outfile, "\n  ==> Properties <==\n\n");
        //    oe->compute();
        //}

        save_information();
    } else {
        if (Communicator::world->me() == 0) {
            fprintf(outfile, "\n  Failed to converged.\n");
            fprintf(outfile, "    NOTE: MO Coefficients will not be saved to Checkpoint.\n");
        }
        E_ = 0.0;
        if(psio_->open_check(PSIF_CHKPT))
            psio_->close(PSIF_CHKPT, 1);
    }

    // Orbitals are always saved, in case a dual basis is required later
    save_orbitals();
    if (options_.get_str("SAPT") != "FALSE") //not a bool because it has types
        save_sapt_info();

    Communicator::world->sync();

    // Clean memory off, handle diis closeout, etc
    finalize();

    //fprintf(outfile,"\nComputation Completed\n");
    fflush(outfile);
    return E_;
}
void BroydenRHF::sr1_seed()
{
    // Form orthonormal guess 
    form_C();
    
    // Sizing may have changed
    nocc_ = nalphapi_;
    nvir_ = nmopi_;
    nvir_ -= nocc_;

    // Grab new preconditioner
    eps_occ_ = epsilon_a_subset("SO","OCC");
    eps_vir_ = epsilon_a_subset("SO","VIR");

    // Reshape F/X registers
    Fia_ = SharedMatrix(new Matrix("Fia", nocc_, nvir_));
    Xia_ = SharedMatrix(new Matrix("Xia", nocc_, nvir_));

    // New orthonormal guess 
    C0_ = SharedMatrix(Ca_->clone());

    // Random stuff 
    broyden_iteration_ = 0;
    V_ = SharedMatrix(new Matrix("V", nmopi_, nmopi_));

    // L-BFGS state vectors
    y_.clear();
    s_.clear();
    p_.clear();

    Fia_old_.reset();
    Xia_old_.reset();
}
void BroydenRHF::sr1_step()
{
    if (broyden_iteration_ == 0) {
        sr1_seed();
        broyden_iteration_++;
    } else {
        // Gradient
        build_Fia();

        // Build new S/Y vectors
        if (Fia_old_) {
            y_.push_back(SharedMatrix(Fia_->clone()));
            y_[y_.size() - 1]->subtract(Fia_old_);
            s_.push_back(SharedMatrix(Xia_->clone()));
            s_[s_.size() - 1]->subtract(Xia_old_);
            p_.push_back(1.0 / y_[y_.size() - 1]->vector_dot(s_[s_.size() - 1]));
            Fia_old_->copy(Fia_);
            Xia_old_->copy(Xia_);
        } else {
            Fia_old_ = SharedMatrix(Fia_->clone());
            Xia_old_ = SharedMatrix(Xia_->clone());
        }        

        std::deque<double> alpha;
        std::deque<double> beta;

        // Produce the effective gradient
        for (int i = ((int)s_.size()) - 1; i >= 0; i--) {
            alpha.push_front(p_[i] * s_[i]->vector_dot(Fia_));
            Xia_->copy(y_[i]);
            Xia_->scale(alpha[0]);
            Fia_->subtract(Xia_);
        }

        // Apply the preconditioner (can be dynamic)
        precondition_Fia();

        for (int i = 0; i < s_.size(); i++) {
            beta.push_back(p_[i] * s_[i]->vector_dot(Fia_));
            Xia_->copy(s_[i]);
            Xia_->scale(alpha[i] - beta[i]);
            Fia_->add(Xia_);
        }

        // And here's the step
        Xia_->copy(Fia_);
        Xia_->scale(-1.0);
        Xia_->add(Xia_old_);

        rotate_orbitals();
        broyden_iteration_++;

        if (s_.size() >= options_.get_int("BROYDEN_MAXITER")) {
            s_.pop_front();
            y_.pop_front();
            p_.pop_front();
        }
    }

    // Set up status message
    std::stringstream ss;
    ss << "SR1 " << (broyden_iteration_);
    broyden_status_ = ss.str();
}
void BroydenRHF::dfp_seed()
{
    // Form orthonormal guess 
    form_C();
    
    // Sizing may have changed
    nocc_ = nalphapi_;
    nvir_ = nmopi_;
    nvir_ -= nocc_;

    // Grab new preconditioner
    eps_occ_ = epsilon_a_subset("SO","OCC");
    eps_vir_ = epsilon_a_subset("SO","VIR");

    // Reshape F/X registers
    Fia_ = SharedMatrix(new Matrix("Fia", nocc_, nvir_));
    Xia_ = SharedMatrix(new Matrix("Xia", nocc_, nvir_));

    // New orthonormal guess 
    C0_ = SharedMatrix(Ca_->clone());

    // Random stuff 
    broyden_iteration_ = 0;
    V_ = SharedMatrix(new Matrix("V", nmopi_, nmopi_));

    // L-BFGS state vectors
    y_.clear();
    s_.clear();
    p_.clear();

    Fia_old_.reset();
    Xia_old_.reset();
}
void BroydenRHF::dfp_step()
{
    if (broyden_iteration_ == 0) {
        dfp_seed();
        broyden_iteration_++;
    } else {
        // Gradient
        build_Fia();

        // Build new S/Y vectors
        if (Fia_old_) {
            y_.push_back(SharedMatrix(Fia_->clone()));
            y_[y_.size() - 1]->subtract(Fia_old_);
            s_.push_back(SharedMatrix(Xia_->clone()));
            s_[s_.size() - 1]->subtract(Xia_old_);
            p_.push_back(1.0 / y_[y_.size() - 1]->vector_dot(s_[s_.size() - 1]));
            Fia_old_->copy(Fia_);
            Xia_old_->copy(Xia_);
        } else {
            Fia_old_ = SharedMatrix(Fia_->clone());
            Xia_old_ = SharedMatrix(Xia_->clone());
        }        

        std::deque<double> alpha;
        std::deque<double> beta;

        // Produce the effective gradient
        for (int i = ((int)s_.size()) - 1; i >= 0; i--) {
            alpha.push_front(p_[i] * s_[i]->vector_dot(Fia_));
            Xia_->copy(y_[i]);
            Xia_->scale(alpha[0]);
            Fia_->subtract(Xia_);
        }

        // Apply the preconditioner (can be dynamic)
        precondition_Fia();

        for (int i = 0; i < s_.size(); i++) {
            beta.push_back(p_[i] * s_[i]->vector_dot(Fia_));
            Xia_->copy(s_[i]);
            Xia_->scale(alpha[i] - beta[i]);
            Fia_->add(Xia_);
        }

        // And here's the step
        Xia_->copy(Fia_);
        Xia_->scale(-1.0);
        Xia_->add(Xia_old_);

        rotate_orbitals();
        broyden_iteration_++;

        if (s_.size() >= options_.get_int("BROYDEN_MAXITER")) {
            s_.pop_front();
            y_.pop_front();
            p_.pop_front();
        }
    }

    // Set up status message
    std::stringstream ss;
    ss << "DFP " << (broyden_iteration_);
    broyden_status_ = ss.str();
}
void BroydenRHF::bfgs_seed()
{
    // Form orthonormal guess 
    form_C();
    
    // Sizing may have changed
    nocc_ = nalphapi_;
    nvir_ = nmopi_;
    nvir_ -= nocc_;

    // Grab new preconditioner
    eps_occ_ = epsilon_a_subset("SO","OCC");
    eps_vir_ = epsilon_a_subset("SO","VIR");

    // Reshape F/X registers
    Fia_ = SharedMatrix(new Matrix("Fia", nocc_, nvir_));
    Xia_ = SharedMatrix(new Matrix("Xia", nocc_, nvir_));

    // New orthonormal guess 
    C0_ = SharedMatrix(Ca_->clone());

    // Random stuff 
    broyden_iteration_ = 0;
    V_ = SharedMatrix(new Matrix("V", nmopi_, nmopi_));

    // L-BFGS state vectors
    y_.clear();
    s_.clear();
    p_.clear();

    Fia_old_.reset();
    Xia_old_.reset();
}
void BroydenRHF::bfgs_step()
{
    if (broyden_iteration_ == 0) {
        bfgs_seed();
        broyden_iteration_++;
    } else {
        // Gradient
        build_Fia();

        // Build new S/Y vectors
        if (Fia_old_) {
            y_.push_back(SharedMatrix(Fia_->clone()));
            y_[y_.size() - 1]->subtract(Fia_old_);
            s_.push_back(SharedMatrix(Xia_->clone()));
            s_[s_.size() - 1]->subtract(Xia_old_);
            p_.push_back(1.0 / y_[y_.size() - 1]->vector_dot(s_[s_.size() - 1]));
            Fia_old_->copy(Fia_);
            Xia_old_->copy(Xia_);
        } else {
            Fia_old_ = SharedMatrix(Fia_->clone());
            Xia_old_ = SharedMatrix(Xia_->clone());
        }        

        std::deque<double> alpha;
        std::deque<double> beta;

        // Produce the effective gradient
        for (int i = ((int)s_.size()) - 1; i >= 0; i--) {
            alpha.push_front(p_[i] * s_[i]->vector_dot(Fia_));
            Xia_->copy(y_[i]);
            Xia_->scale(alpha[0]);
            Fia_->subtract(Xia_);
        }

        // Apply the preconditioner (can be dynamic)
        precondition_Fia();

        for (int i = 0; i < s_.size(); i++) {
            beta.push_back(p_[i] * s_[i]->vector_dot(Fia_));
            Xia_->copy(s_[i]);
            Xia_->scale(alpha[i] - beta[i]);
            Fia_->add(Xia_);
        }

        // And here's the step
        Xia_->copy(Fia_);
        Xia_->scale(-1.0);
        Xia_->add(Xia_old_);

        rotate_orbitals();
        broyden_iteration_++;

        if (s_.size() >= options_.get_int("BROYDEN_MAXITER")) {
            s_.pop_front();
            y_.pop_front();
            p_.pop_front();
        }
    }

    // Set up status message
    std::stringstream ss;
    ss << "BFGS " << (broyden_iteration_);
    broyden_status_ = ss.str();
}
void BroydenRHF::broyden_seed()
{
    // Form orthonormal guess 
    form_C();
    
    // Sizing may have changed
    nocc_ = nalphapi_;
    nvir_ = nmopi_;
    nvir_ -= nocc_;

    // Grab new preconditioner
    eps_occ_ = epsilon_a_subset("SO","OCC");
    eps_vir_ = epsilon_a_subset("SO","VIR");

    // Reshape F/X registers
    Fia_ = SharedMatrix(new Matrix("Fia", nocc_, nvir_));
    Xia_ = SharedMatrix(new Matrix("Xia", nocc_, nvir_));

    // New orthonormal guess 
    C0_ = SharedMatrix(Ca_->clone());

    // Random stuff 
    broyden_iteration_ = 0;
    V_ = SharedMatrix(new Matrix("V", nmopi_, nmopi_));

    // Broyden state vectors
    s_.clear();
    s2_.clear();
}
void BroydenRHF::broyden_step()
{
    if (broyden_iteration_ == 0) {
        broyden_seed();
        broyden_iteration_++;
    } else {
        build_Fia();
        precondition_Fia();

        printf("||F_ia||_RMS = %24.16E\n", Fia_->rms());
        
        SharedMatrix Z(Fia_->clone());
        Z->scale(-1.0);

        for (int j = 0; j < ((int)(s_.size()) - 1); j++) {
            double val = Z->vector_dot(s_[j]) / s2_[j];
            Fia_->copy(s_[j+1]);
            Fia_->scale(val);
            Z->add(Fia_);
        }
        
        double norm = 1.0;
        if (s_.size()) {
            norm = 1 - Z->vector_dot(s_[s_.size() - 1]) / s2_[s_.size() - 1];
        }
        Z->scale(1.0 / norm);

        s_.push_back(SharedMatrix(Z->clone()));
        s2_.push_back(s_[s_.size() - 1]->vector_dot(s_[s_.size() - 1]));
        
        Xia_->add(s_[s_.size() - 1]);
        rotate_orbitals();
        broyden_iteration_++;

        if (s_.size() >= options_.get_int("BROYDEN_MAXITER")) {
            s_.pop_front();
            s2_.pop_front();
        }
    }

    // Set up status message
    std::stringstream ss;
    ss << "BROYDEN " << (broyden_iteration_);
    broyden_status_ = ss.str();
}
void BroydenRHF::precondition_Fia()
{
    // Diagonal 
    for (int h = 0; h < Fia_->nirrep(); h++) {
        int ni = nocc_[h];
        int na = nvir_[h];
        if (!ni || !na) continue;

        double** F2p = Fia_->pointer(h);
        double* eps_a = eps_vir_->pointer(h);
        double* eps_i = eps_occ_->pointer(h);

        for (int i = 0; i < ni; i++) {
            for (int a = 0; a < na; a++) {
                F2p[i][a] /= (eps_a[a] - eps_i[i]);
            }
        }

    }
}
void BroydenRHF::build_Fia()
{
    // Compute -Fia = C_im F_mn C_na (Orbital Gradient)
    double* temp = new double[nocc_.max() * nmopi_.max()];
    for (int h = 0; h < Fia_->nirrep(); h++) {
        int ni = nocc_[h];
        int na = nvir_[h];
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        if (!ni || !na) continue;

        double** Cp = Ca_->pointer(h);
        double** Fp = Fa_->pointer(h);
        double** F2p = Fia_->pointer(h);

        C_DGEMM('T','N',ni,nso,nso,1.0,&Cp[0][0],nmo,Fp[0],nso,0.0,temp,nso);
        C_DGEMM('N','N',ni,na,nso,-1.0,temp,nso,&Cp[0][ni],nmo,0.0,F2p[0],na); 

    }
    delete[] temp;
} 
void BroydenRHF::rotate_orbitals()
{
    // Rotate Ca by exp(X0) (Uses truncated Pade expansion)
    V_->zero();
    for (int h = 0; h < V_->nirrep(); h++) {
        int ni = nocc_[h];
        int na = nvir_[h];
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        if (!ni || !na) continue;
       
        double** Xp = Xia_->pointer(h);
        double** Yp = V_->pointer(h);
        
        for (int i = 0; i < ni; i++) {
            C_DAXPY(na, 1.0,Xp[i],1,&Yp[i][ni],1);
            C_DAXPY(na,-1.0,Xp[i],1,&Yp[ni][i],nmo);
        } 
    }

    //V_->print();
    V_->expm(options_.get_int("BROYDEN_PADE_N"));
    //V_->print();
    Ca_->zero();
    for (int h = 0; h < V_->nirrep(); h++) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        if (!nmo) continue;

        double** C0p = C0_->pointer(h);
        double** Vp = V_->pointer(h);
        double** C2p = Ca_->pointer(h); 
        
        C_DGEMM('N','N',nso,nmo,nmo,1.0,C0p[0],nmo,Vp[0],nmo,0.0,C2p[0],nmo);
    }

    //Ca_->print();
}

}} // Namespaces
