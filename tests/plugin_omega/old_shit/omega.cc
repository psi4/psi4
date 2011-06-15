#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <string>
#include <cstring>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include <libfunctional/superfunctional.h>
#include <libscf_solver/ks.h>
#include <libscf_solver/integralfunctors.h>
#include <libscf_solver/omegafunctors.h>

#include "omega.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;
using namespace psi::functional;
using namespace boost;

namespace psi{ namespace scf {

OmegaKS::OmegaKS(Options & options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) :
    UKS(options, psio, chkpt)
{
    common_init();
}
OmegaKS::OmegaKS(Options & options, boost::shared_ptr<PSIO> psio) :
    UKS(options, psio)
{
    common_init();
}
OmegaKS::~OmegaKS()
{
}
void OmegaKS::common_init()
{
    E_N_ = 0.0;
    Eold_N_ = 0.0;
    Drms_N_ = 0.0;

    ::memset(static_cast<void*>(doccpi_N_), '\0', 8 * sizeof(int));
    ::memset(static_cast<void*>(soccpi_N_), '\0', 8 * sizeof(int));
    ::memset(static_cast<void*>(nalphapi_N_), '\0', 8 * sizeof(int));
    ::memset(static_cast<void*>(nbetapi_N_), '\0', 8 * sizeof(int));

    // Doublet if singlet, m+1 otherwise
    if (multiplicity_ == 1) {
        nalpha_N_ = nalpha_; 
        nbeta_N_ = nbeta_ - 1;
    } else {
        nalpha_N_ = nalpha_ - 1; 
        nbeta_N_ = nbeta_;
    }

    epsilon_a_N_ = boost::shared_ptr<Vector>(new Vector("Epsilon_a_N", nirrep_, nmopi_));
    epsilon_b_N_ = boost::shared_ptr<Vector>(new Vector("Epsilon_b_N", nirrep_, nmopi_));

    Ca_N_ = boost::shared_ptr<Matrix>(new Matrix("Ca_N", nirrep_, nsopi_, nmopi_));
    Cb_N_ = boost::shared_ptr<Matrix>(new Matrix("Cb_N", nirrep_, nsopi_, nmopi_));
    Da_N_ = boost::shared_ptr<Matrix>(new Matrix("Da_N", nirrep_, nsopi_, nsopi_));
    Db_N_ = boost::shared_ptr<Matrix>(new Matrix("Db_N", nirrep_, nsopi_, nsopi_));
    Dt_N_ = boost::shared_ptr<Matrix>(new Matrix("Dt_N", nirrep_, nsopi_, nsopi_));
    Dtold_N_ = boost::shared_ptr<Matrix>(new Matrix("Dtold_N", nirrep_, nsopi_, nsopi_));
    J_N_ = boost::shared_ptr<Matrix>(new Matrix("J_N", nirrep_, nsopi_, nsopi_));
    Ka_N_ = boost::shared_ptr<Matrix>(new Matrix("Ka_N", nirrep_, nsopi_, nsopi_));
    Kb_N_ = boost::shared_ptr<Matrix>(new Matrix("Kb_N", nirrep_, nsopi_, nsopi_));
    wKa_N_ = boost::shared_ptr<Matrix>(new Matrix("wKa_N", nirrep_, nsopi_, nsopi_));
    wKb_N_ = boost::shared_ptr<Matrix>(new Matrix("wKb_N", nirrep_, nsopi_, nsopi_));
    Va_N_ = boost::shared_ptr<Matrix>(new Matrix("Va_N", nirrep_, nsopi_, nsopi_));
    Vb_N_ = boost::shared_ptr<Matrix>(new Matrix("Vb_N", nirrep_, nsopi_, nsopi_));
    Ga_N_ = boost::shared_ptr<Matrix>(new Matrix("Ga_N", nirrep_, nsopi_, nsopi_));
    Gb_N_ = boost::shared_ptr<Matrix>(new Matrix("Gb_N", nirrep_, nsopi_, nsopi_));
    Fa_N_ = boost::shared_ptr<Matrix>(new Matrix("Fa_N", nirrep_, nsopi_, nsopi_));
    Fb_N_ = boost::shared_ptr<Matrix>(new Matrix("Fb_N", nirrep_, nsopi_, nsopi_));
}
void OmegaKS::initialize_N()
{
    ::memcpy(static_cast<void*>(doccpi_N_), static_cast<void*>(doccpi_), 8 * sizeof(int));
    ::memcpy(static_cast<void*>(soccpi_N_), static_cast<void*>(soccpi_), 8 * sizeof(int));
    ::memcpy(static_cast<void*>(nalphapi_N_), static_cast<void*>(nalphapi_), 8 * sizeof(int));
    ::memcpy(static_cast<void*>(nbetapi_N_), static_cast<void*>(nbetapi_), 8 * sizeof(int));

    // Use the N-electron wavefunction as a guess
    Ca_N_->copy(Ca_);
    Cb_N_->copy(Cb_);

    find_occupation_N();

    form_D_N();

    E_N_ = E_;
}

double OmegaKS::compute_energy()
{

    if (print_ & (Communicator::world->me() == 0))
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
        fprintf(outfile, "  ==> Burn-In N-Electron Iterations <==\n\n");
        fprintf(outfile, "                        Total Energy        Delta E      Density RMS\n\n");
    }
    fflush(outfile);

    // Neither of these are idempotent
    if ((KS::options_.get_str("GUESS") == "SAD") || (KS::options_.get_str("GUESS") == "READ"))
        iteration_ = -1;
    else
        iteration_ = 0;

    bool converged = false;
    diis_performed_ = false;
   
    do {
        iteration_++;

        save_density_and_energy();

        timer_on("Form G");
        form_G();
        timer_off("Form G");

        // Reset fractional SAD occupation
        if (iteration_ == 0 && KS::options_.get_str("GUESS") == "SAD")
            reset_SAD_occupation();

        timer_on("Form F");
        form_F();
        timer_off("Form F");

        if (print_>3) {
            Fa_->print(outfile);
            Fb_->print(outfile);
        }

        E_ = compute_E();

        timer_on("DIIS");
        if (diis_enabled_ && iteration_ > 0 && iteration_ >= diis_start_ )
            save_fock();
        if (diis_enabled_ == true && iteration_ >= diis_start_ + min_diis_vectors_ - 1) {
            diis_performed_ = diis();
        } else {
            diis_performed_ = false;
        }
        timer_off("DIIS");

        if (print_>4 && diis_performed_ & (Communicator::world->me() == 0)) {
            fprintf(outfile,"  After DIIS:\n");
            Fa_->print(outfile);
            Fb_->print(outfile);
        }

        std::string status;
        if (!diis_performed_) status = "";
        else status = "DIIS";

        if (Communicator::world->me() == 0) {
            fprintf(outfile, "   UKS iter %3d: %20.14f   % 10.5e   % 10.5e %s\n",
                              iteration_, E_, E_ - Eold_, Drms_, status.c_str());
            fflush(outfile);
        }

        timer_on("Form C");
        form_C();
        timer_off("Form C");

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

    } while (!converged && iteration_ < maxiter_ );

    if (!converged) throw PSIEXCEPTION("Burn-in N-electron SCF did not converge");

    initialize_N();

    if (Communicator::world->me() == 0) {
        fprintf(outfile, "  ==> Burn-In N-1-Electron Iterations <==\n\n");
        fprintf(outfile, "                        Total Energy        Delta E      Density RMS\n\n");
    }
    fflush(outfile);

    iteration_ = 0;
    converged = false;
    diis_performed_ = false;
   
    do {
        iteration_++;

        save_density_and_energy_N();

        timer_on("Form G");
        form_G_N();
        timer_off("Form G");

        timer_on("Form F");
        form_F_N();
        timer_off("Form F");

        if (print_>3) {
            Fa_N_->print(outfile);
            Fb_N_->print(outfile);
        }

        E_N_ = compute_E_N();
        
        // TODO: DIIS

        //timer_on("DIIS");
        //if (diis_enabled_ && iteration_ > 0 && iteration_ >= diis_start_ )
        //    save_fock();
        //if (diis_enabled_ == true && iteration_ >= diis_start_ + min_diis_vectors_ - 1) {
        //    diis_performed_ = diis();
        //} else {
        //    diis_performed_ = false;
        //}
        //timer_off("DIIS");

        if (print_>4 && diis_performed_ & (Communicator::world->me() == 0)) {
            fprintf(outfile,"  After DIIS:\n");
            Fa_N_->print(outfile);
            Fb_N_->print(outfile);
        }

        std::string status;
        if (!diis_performed_) status = "";
        else status = "DIIS";

        if (Communicator::world->me() == 0) {
            fprintf(outfile, "   UKS iter %3d: %20.14f   % 10.5e   % 10.5e %s\n",
                              iteration_, E_N_, E_N_ - Eold_N_, Drms_N_, status.c_str());
            fflush(outfile);
        }

        timer_on("Form C");
        form_C_N();
        timer_off("Form C");

        timer_on("Form D");
        form_D_N();
        timer_off("Form D");

        if (print_ > 3){
            Ca_N_->print(outfile);
            Cb_N_->print(outfile);
            Da_N_->print(outfile);
            Db_N_->print(outfile);
        }

        converged = test_convergency_N();

    } while (!converged && iteration_ < maxiter_ );

    if (!converged) throw PSIEXCEPTION("Burn-in N-1-electron SCF did not converge");

    // TODO: Main crazy loop

    if (Communicator::world->me() == 0) 
        fprintf(outfile, "\n  ==> Post-Iterations <==\n\n");

    if (converged) {
        // Need to recompute the Fock matrices, as they are modified during the SCF interation
        // and might need to be dumped to checkpoint later
        form_F();

        // Print the orbitals
        if(print_)
            print_orbitals();

        if (Communicator::world->me() == 0) {
            fprintf(outfile, "\n  Energy converged.\n");
            fprintf(outfile, "\n  UKS Final Energy: %20.14f", E_);
            if (perturb_h_) {
                fprintf(outfile, " with %f perturbation", lambda_);
            }
            fprintf(outfile, "\n");
        }

        // TODO
        // Properties
        //if (print_) {
        //    boost::shared_ptr<OEProp> oe(new OEProp(boost::shared_ptr<Wavefunction>(this)));
        //    oe->add("DIPOLE");

        //    if (print_ >= 2) {
        //        oe->add("QUADRUPOLE");
        //        oe->add("MULLIKEN_CHARGES");
        //    }

        //    if (Communicator::world->me() == 0)
        //        fprintf(outfile, "\n  ==> Properties <==\n");
        //    oe->compute();
        //}

        save_information();
    } else {
        if (Communicator::world->me() == 0) {
            fprintf(outfile, "\n  Failed to converged.\n");
            fprintf(outfile, "    NOTE: MO Coefficients will not be saved to Checkpoint.\n");
        }
        E_ = 0.0;
        if(KS::psio_->open_check(PSIF_CHKPT))
            KS::psio_->close(PSIF_CHKPT, 1);
    }

    // Orbitals are always saved, in case a dual basis is required later
    save_orbitals();

    // Clean memory off, handle diis closeout, etc
    finalize();

    fflush(outfile);
    return E_;
}
void OmegaKS::finalize()
{
    UKS::finalize();
}
void OmegaKS::save_density_and_energy_N()
{
    Dtold_N_ = Dt_N_;
    Eold_N_ = E_N_; 
}
bool OmegaKS::test_convergency_N()
{
    // energy difference
    double ediff = E_N_ - Eold_N_;

    // RMS of the density
    Matrix D_rms;
    D_rms.copy(Dt_N_);
    D_rms.subtract(Dtold_N_);
    Drms_N_ = D_rms.rms();

    if (fabs(ediff) < energy_threshold_ && Drms_N_ < density_threshold_)
        return true;
    else
        return false;
}
double OmegaKS::compute_E_N()
{
    // E_DFT = 2.0 D*H + 2.0 D*J - \alpha D*K + E_xc
    double one_electron_E = Da_N_->vector_dot(H_);
    one_electron_E += Db_N_->vector_dot(H_);
    double coulomb_E = Da_N_->vector_dot(J_N_);
    coulomb_E += Db_N_->vector_dot(J_N_);

    double exchange_E = 0.0;
    if (functional_->isHybrid()) {
        exchange_E = -functional_->getExactExchange()*Da_N_->vector_dot(Ka_N_);
        exchange_E = -functional_->getExactExchange()*Db_N_->vector_dot(Kb_N_);
    }
    if (functional_->isRangeCorrected()) {
        exchange_E = -Da_N_->vector_dot(wKa_N_);
        exchange_E = -Db_N_->vector_dot(wKb_N_);
    }
    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += 0.5*coulomb_E;
    Etotal += exchange_E;
    Etotal += quad_values_["E_xc"];
    if (functional_->isDashD()) {
        double dashD_E = functional_->getDashD()->computeEnergy(HF::molecule_);
        Etotal += dashD_E;
    }

    if (print_ > 2) {
        fprintf(outfile, "N-1 Nuclear Repulsion energy: %24.16f\n", nuclearrep_);
        fprintf(outfile, "N-1 One-electron energy:      %24.16f\n", one_electron_E);
        fprintf(outfile, "N-1 Coulomb energy:           %24.16f\n", 0.5*coulomb_E);
        fprintf(outfile, "N-1 Exchange energy:          %24.16f\n", exchange_E);
        fprintf(outfile, "N-1 Functional energy:        %24.16f\n", quad_values_["E_xc_N"]);
    }
    return Etotal;
}
void OmegaKS::form_G_N()
{
    form_V_N();
    if (functional_->isRangeCorrected()) {
        Omega_Ka_Kb_Functor k_builder(functional_->getOmega(),wKa_N_,wKb_N_,Da_N_,Db_N_,Ca_N_,Cb_N_,nalphapi_N_,nbetapi_N_);
        process_omega_tei<Omega_Ka_Kb_Functor>(k_builder);
    }
    if (!functional_->isHybrid()) {
        // This will build J (stored in G)
        J_Functor j_builder(Ga_N_, Da_N_);
        process_tei<J_Functor>(j_builder);
        J_N_->copy(Ga_N_);

        Gb_N_->copy(Ga_N_);
        Ga_N_->add(Va_N_);
        Gb_N_->add(Vb_N_);
        Ga_N_->subtract(wKa_N_);
        Gb_N_->subtract(wKb_N_);

    } else {
        // This will build J (stored in G) and K
        J_Ka_Kb_Functor jk_builder(Ga_N_, Ka_N_, Kb_N_, Da_N_, Db_N_, Ca_N_, Cb_N_, nalphapi_N_, nbetapi_N_);
        process_tei<J_Ka_Kb_Functor>(jk_builder);
        J_N_->copy(Ga_N_);
        Gb_N_->copy(Ga_N_);

        double alpha = functional_->getExactExchange();
        Ka_N_->scale(alpha);
        Kb_N_->scale(alpha);
        Ga_N_->subtract(Ka_N_);
        Gb_N_->subtract(Kb_N_);
        Ka_N_->scale(1.0/alpha);
        Kb_N_->scale(1.0/alpha);
        Ga_N_->add(Va_N_);
        Gb_N_->add(Vb_N_);
        Ga_N_->subtract(wKa_N_);
        Gb_N_->subtract(wKb_N_);
    }
}
void OmegaKS::form_F_N()
{
    Fa_N_->copy(H_);
    Fa_N_->add(Ga_N_);

    Fb_N_->copy(H_);
    Fb_N_->add(Gb_N_);

    if (debug_) {
        Fa_N_->print(outfile);
        Fb_N_->print(outfile);
    }
}
void OmegaKS::form_C_N()
{
    diagonalize_F(Fa_N_, Ca_N_, epsilon_a_N_);
    diagonalize_F(Fb_N_, Cb_N_, epsilon_b_N_);
    find_occupation_N();
    if (debug_) {
        Ca_N_->print(outfile);
        Cb_N_->print(outfile);
    }
}
void OmegaKS::form_D_N()
{
    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = nalphapi_N_[h];
        int nb = nbetapi_N_[h];
    
        if (nso == 0 || nmo == 0) continue;

        double** Ca = Ca_N_->pointer(h);
        double** Cb = Cb_N_->pointer(h);
        double** Da = Da_N_->pointer(h);
        double** Db = Db_N_->pointer(h);

        if (na == 0) 
            memset(static_cast<void*>(Da[0]), '\0', sizeof(double)*nso*nso);
        if (nb == 0) 
            memset(static_cast<void*>(Db[0]), '\0', sizeof(double)*nso*nso);

        C_DGEMM('N','T',nso,nso,na,1.0,Ca[0],nmo,Ca[0],nmo,0.0,Da[0],nso);
        C_DGEMM('N','T',nso,nso,nb,1.0,Cb[0],nmo,Cb[0],nmo,0.0,Db[0],nso);

    }

    Dt_N_->copy(Da_N_);
    Dt_N_->add(Db_N_);

    if (debug_) {
        fprintf(outfile, "in UHF::form_D_N:\n");
        Da_N_->print();
        Db_N_->print();
    }
}
void OmegaKS::find_occupation_N()
{
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<epsilon_a_N_->nirrep(); ++h) {
        for (int i=0; i<epsilon_a_N_->dimpi()[h]; ++i)
            pairs.push_back(make_pair(epsilon_a_N_->get(h, i), h));
    }
    sort(pairs.begin(),pairs.end());

    // Store the current occupation
    int old_socc[8];
    int old_docc[8];
    for(int h = 0; h < nirrep_; ++h){
        old_socc[h] = soccpi_N_[h];
        old_docc[h] = doccpi_N_[h];
    }

    memset(doccpi_N_, 0, sizeof(int) * epsilon_a_N_->nirrep());
    for (int i=0; i<nbeta_N_; ++i)
        doccpi_N_[pairs[i].second]++;
    memset(soccpi_N_, 0, sizeof(int) * epsilon_a_N_->nirrep());
    for (int i=nbeta_N_; i<nalpha_N_; ++i)
        soccpi_[pairs[i].second]++;

    for (int i=0; i<epsilon_a_N_->nirrep(); ++i) {
        nalphapi_N_[i] = doccpi_N_[i] + soccpi_N_[i];
        nbetapi_N_[i]  = doccpi_N_[i];
    }

    bool occ_changed = false;
    for(int h = 0; h < nirrep_; ++h){
        if( old_socc[h] != soccpi_N_[h] || old_docc[h] != doccpi_N_[h]){
            occ_changed = true;
            break;
        }
    }
}
void OmegaKS::form_V_N()
{
}

}} // End Namespaces
