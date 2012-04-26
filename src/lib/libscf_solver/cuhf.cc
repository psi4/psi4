#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include <libfock/jk.h>
#include "integralfunctors.h"

#include "cuhf.h"

using namespace std;
using namespace psi;
using namespace boost;

namespace psi { namespace scf {

CUHF::CUHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) : HF(options, psio, chkpt)
{
    common_init();
}

CUHF::CUHF(Options& options, boost::shared_ptr<PSIO> psio) : HF(options, psio)
{
    common_init();
}

CUHF::~CUHF()
{
}

void CUHF::common_init()
{

    Drms_ = 0.0;

    Fa_         = SharedMatrix(factory_->create_matrix("F alpha"));
    Fb_         = SharedMatrix(factory_->create_matrix("F beta"));
    Fp_         = SharedMatrix(factory_->create_matrix("F charge"));
    Fm_         = SharedMatrix(factory_->create_matrix("F spin"));
    Da_         = SharedMatrix(factory_->create_matrix("D alpha"));
    Db_         = SharedMatrix(factory_->create_matrix("D beta"));
    Dp_         = SharedMatrix(factory_->create_matrix("D charge"));
    Dt_         = SharedMatrix(factory_->create_matrix("D total"));
    Dtold_      = SharedMatrix(factory_->create_matrix("D total old"));
    Lagrangian_ = SharedMatrix(factory_->create_matrix("Lagrangian"));
    Ca_         = SharedMatrix(factory_->create_matrix("C alpha"));
    Cb_         = SharedMatrix(factory_->create_matrix("C beta"));
    Cno_        = SharedMatrix(factory_->create_matrix("C NOs"));
    Cno_temp_   = SharedMatrix(factory_->create_matrix("C NO temp"));
    J_          = SharedMatrix(factory_->create_matrix("J total"));
    Ka_         = SharedMatrix(factory_->create_matrix("K alpha"));
    Kb_         = SharedMatrix(factory_->create_matrix("K beta"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = SharedVector(factory_->create_vector());
    No_ = SharedVector(factory_->create_vector());

}

void CUHF::finalize()
{
    // Form lagrangian
    for (int h=0; h<nirrep_; ++h) {
        for (int m=0; m<Lagrangian_->rowdim(h); ++m) {
            for (int n=0; n<Lagrangian_->coldim(h); ++n) {
                double sum = 0.0;
                for (int i=0; i<doccpi_[h]; ++i) {
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i)
                        +  epsilon_b_->get(h, i) * Cb_->get(h, m, i) * Cb_->get(h, n, i);
                }
                for (int i=doccpi_[h]; i<doccpi_[h]+soccpi_[h]; ++i)
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i);

                Lagrangian_->set(h, m, n, sum);
            }
        }
    }

    Dt_.reset();
    Dtold_.reset();
    Dp_.reset();
    Fp_.reset();
    Fm_.reset();
    Cno_.reset();
    Cno_temp_.reset();
    No_.reset();

    HF::finalize();
}

void CUHF::save_density_and_energy()
{
    Dtold_->copy(Dt_);
    Eold_ = E_;
}

void CUHF::form_G()
{
    if (scf_type_ == "DF" || scf_type_ == "PS") {

        // Push the C matrix on
        std::vector<SharedMatrix> & C = jk_->C_left();
        C.clear();
        C.push_back(Ca_subset("SO", "OCC"));
        C.push_back(Cb_subset("SO", "OCC"));
        
        // Run the JK object
        jk_->compute();

        // Pull the J and K matrices off
        const std::vector<SharedMatrix> & J = jk_->J();
        const std::vector<SharedMatrix> & K = jk_->K();
        J_->copy(J[0]);
        J_->add(J[1]);
        Ka_ = K[0];
        Kb_ = K[1];
        
    } else {
        J_Ka_Kb_Functor jk_builder(J_, Ka_, Kb_, Da_, Db_, Ca_, Cb_, nalphapi_, nbetapi_);
        process_tei<J_Ka_Kb_Functor>(jk_builder);
    }
}

void CUHF::save_information()
{
}

void CUHF::compute_spin_contamination()
{
    double dN = 0.0;

    for (int h =0; h < S_->nirrep(); h++) {
        int nbf = S_->colspi()[h];
        int nmo = Ca_->colspi()[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];
        if (na == 0 || nb == 0 || nbf == 0 || nmo == 0)
            continue;

        SharedMatrix Ht (new Matrix("H Temp", nbf, nb));
        SharedMatrix Ft (new Matrix("F Temp", na, nb));

        double** Sp = S_->pointer(h);
        double** Cap = Ca_->pointer(h);
        double** Cbp = Cb_->pointer(h);
        double** Htp = Ht->pointer(0);
        double** Ftp = Ft->pointer(0);

        C_DGEMM('N','N',nbf,nb,nbf,1.0,Sp[0],nbf,Cbp[0],nmo,0.0,Htp[0],nb);
        C_DGEMM('T','N',na,nb,nbf,1.0,Cap[0],nmo,Htp[0],nb,0.0,Ftp[0],nb);

        for (long int ab = 0; ab < (long int)na*nb; ab++)
            dN += Ftp[0][ab]*Ftp[0][ab];
    }

    double dS = (double)nbeta_ - (double)dN;

    double nm = (nalpha_ - nbeta_) / 2.0;
    double S2 = nm * (nm + 1.0);

    fprintf(outfile, "\n  @Spin Contamination Metric: %8.5F\n", dS);
      fprintf(outfile, "  @S^2 Expected:              %8.5F\n", S2);
      fprintf(outfile, "  @S^2 Observed:              %8.5F\n", S2 + dS);

}

bool CUHF::test_convergency()
{
    double ediff = E_ - Eold_;

    // RMS of the density
    Matrix Drms;
    Drms.copy(Dt_);
    Drms.subtract(Dtold_);
    Drms_ = Drms.rms();

    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}

void CUHF::form_initialF()
{
    Fa_->copy(H_);
    Fb_->copy(H_);

    if (debug_) {
        fprintf(outfile, "Initial Fock alpha matrix:\n");
        Fa_->print(outfile);
        fprintf(outfile, "Initial Fock beta matrix:\n");
        Fb_->print(outfile);
    }
}

void CUHF::form_F()
{
    // Form (rho_a + rho_b) / 2
    Dp_->copy(Da_);
    Dp_->add(Db_);
    Dp_->scale(-0.5); // This is a hack to get the eigenvectors in the
                      // order that I want
    if (debug_) {
      fprintf(outfile, "Charge Density Matrix (SO Basis):\n");
      Dp_->print();
    }

    // Transfrom to an orthonormal basis, C_a is convenient
    Dp_->transform(S_);
    Dp_->transform(Ca_);
    if (debug_) {
      fprintf(outfile, "Charge Density Matrix (Alpha Basis):\n");
      Dp_->print();
    }

    // Diagonalize the charge density and form the natural orbitals
    Dp_->diagonalize(Cno_temp_,No_);
    if (debug_) {
      fprintf(outfile, "CUHF Natural Orbital Occupations:\n");
      No_->print();
    }
    Cno_->gemm(false, false, 1.0, Ca_, Cno_temp_, 0.0);

    // Now we form the contributions to the Fock matrix from
    // the charge and spin densities
    Fp_->copy(J_);
    Fp_->scale(2.0);
    Fp_->subtract(Ka_);
    Fp_->subtract(Kb_);
    Fp_->scale(0.5);

    Fm_->copy(Ka_);
    Fm_->subtract(Kb_);
    Fm_->scale(-0.5);

    // Transform the spin density contributions to the NO basis
    Fm_->transform(Cno_);

    // Zero the core-virtual contributions
    //
    //            [ Fm_cc Fm_co   0   ]
    // Fm_tilde = [ Fm_oc Fm_oo Fm_ov ]
    //            [   0   Fm_vo Fm_vv ]
    //
    for (int h = 0; h < nirrep_; ++h) {
      for (int i = 0; i < doccpi_[h]; ++i) {
        for (int j = doccpi_[h] + soccpi_[h]; j < nmopi_[h]; ++j) {
          Fm_->set(h, i, j, 0.0);
          Fm_->set(h, j, i, 0.0);
        }
      }
    }

    // Return to the SO basis
    Fm_->back_transform(Cno_);
    Fm_->transform(S_);

    // Build the modified alpha and beta Fock matrices
    Fa_->copy(H_);
    Fa_->add(Fp_);
    Fa_->add(Fm_);

    Fb_->copy(H_);
    Fb_->add(Fp_);
    Fb_->subtract(Fm_);

    if (debug_) {
        Fa_->print(outfile);
        Fb_->print(outfile);
    }
}

void CUHF::form_C()
{
    diagonalize_F(Fa_, Ca_, epsilon_a_);
    diagonalize_F(Fb_, Cb_, epsilon_b_);
    find_occupation();
    if (debug_) {
        Ca_->print(outfile);
        Cb_->print(outfile);
    }
}

void CUHF::form_D()
{
    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];

        if (nso == 0 || nmo == 0) continue;

        double** Ca = Ca_->pointer(h);
        double** Cb = Cb_->pointer(h);
        double** Da = Da_->pointer(h);
        double** Db = Db_->pointer(h);

        if (na == 0)
            memset(static_cast<void*>(Da[0]), '\0', sizeof(double)*nso*nso);
        if (nb == 0)
            memset(static_cast<void*>(Db[0]), '\0', sizeof(double)*nso*nso);

        C_DGEMM('N','T',nso,nso,na,1.0,Ca[0],nmo,Ca[0],nmo,0.0,Da[0],nso);
        C_DGEMM('N','T',nso,nso,nb,1.0,Cb[0],nmo,Cb[0],nmo,0.0,Db[0],nso);

    }

    Dt_->copy(Da_);
    Dt_->add(Db_);

    if (debug_) {
        fprintf(outfile, "in CUHF::form_D:\n");
        Da_->print();
        Db_->print();
    }
}

// TODO: Once Dt_ is refactored to D_ the only difference between this and RHF::compute_initial_E is a factor of 0.5
double CUHF::compute_initial_E()
{
    return nuclearrep_ + 0.5 * (Dt_->vector_dot(H_));
}

double CUHF::compute_E()
{
    double one_electron_E = Dt_->vector_dot(H_);
    double two_electron_E = 0.5 * (Da_->vector_dot(Fa_) + Db_->vector_dot(Fb_) - one_electron_E);   
 
    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = two_electron_E; 
    energies_["XC"] = 0.0;
    energies_["-D"] = 0.0;

    double DH  = Dt_->vector_dot(H_);
    double DFa = Da_->vector_dot(Fa_);
    double DFb = Db_->vector_dot(Fb_);
    double Eelec = 0.5 * (DH + DFa + DFb);
    // fprintf(outfile, "electronic energy = %20.14f\n", Eelec);
    double Etotal = nuclearrep_ + Eelec;
    return Etotal;
}

void CUHF::save_fock()
{
    SharedMatrix FDSmSDFa = form_FDSmSDF(Fa_, Da_);
    SharedMatrix FDSmSDFb = form_FDSmSDF(Fb_, Db_);
//    FDSmSDFa->add(FDSmSDFb);

    if (initialized_diis_manager_ == false) {
        diis_manager_ = boost::shared_ptr<DIISManager>(new DIISManager(
            max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError,
            DIISManager::OnDisk));
//        diis_manager_->set_error_vector_size(1, DIISEntry::Matrix,
//            FDSmSDFa.get());
        diis_manager_->set_error_vector_size(2, DIISEntry::Matrix,
            FDSmSDFa.get(), DIISEntry::Matrix, FDSmSDFb.get());
        diis_manager_->set_vector_size(2, DIISEntry::Matrix,
            Fa_.get(), DIISEntry::Matrix, Fb_.get());
        initialized_diis_manager_ = true;
    }

//    diis_manager_->add_entry(3, FDSmSDFa.get(), Fa_.get(), Fb_.get());
    diis_manager_->add_entry(4, FDSmSDFa.get(), FDSmSDFb.get(), Fa_.get(),
        Fb_.get());
}

bool CUHF::diis()
{
    return diis_manager_->extrapolate(2, Fa_.get(), Fb_.get());
}

void CUHF::stability_analysis()
{
    throw PSIEXCEPTION("CUHF stability analysis has not been implemented yet.  Sorry :(");
}

}}
