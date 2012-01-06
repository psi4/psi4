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
#include "integralfunctors.h"

#include <libmints/mints.h>
#include "uhf.h"

using namespace std;
using namespace psi;
using namespace boost;

namespace psi { namespace scf {

UHF::UHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt) : HF(options, psio, chkpt)
{
    common_init();
}

UHF::UHF(Options& options, boost::shared_ptr<PSIO> psio) : HF(options, psio)
{
    common_init();
}

UHF::~UHF()
{
}

void UHF::common_init()
{

    Drms_ = 0.0;

    Fa_     = SharedMatrix(factory_->create_matrix("F alpha"));
    Fb_     = SharedMatrix(factory_->create_matrix("F beta"));
    Da_     = SharedMatrix(factory_->create_matrix("D alpha"));
    Db_     = SharedMatrix(factory_->create_matrix("D beta"));
    Dt_     = SharedMatrix(factory_->create_matrix("D total"));
    Dtold_  = SharedMatrix(factory_->create_matrix("D total old"));
    Lagrangian_ = SharedMatrix(factory_->create_matrix("Lagrangian"));
    Ca_     = SharedMatrix(factory_->create_matrix("C alpha"));
    Cb_     = SharedMatrix(factory_->create_matrix("C beta"));
    Ga_     = SharedMatrix(factory_->create_matrix("G alpha"));
    Gb_     = SharedMatrix(factory_->create_matrix("G beta"));
    J_      = SharedMatrix(factory_->create_matrix("J total"));
    Ka_     = SharedMatrix(factory_->create_matrix("K alpha"));
    Kb_     = SharedMatrix(factory_->create_matrix("K beta"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = SharedVector(factory_->create_vector());

}

void UHF::finalize()
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
    Ga_.reset();
    Gb_.reset();

    HF::finalize();
}

void UHF::save_density_and_energy()
{
    Dtold_->copy(Dt_);
    Eold_ = E_;
}

void UHF::form_G()
{
    // This will build J (stored in G) and K
    J_Ka_Kb_Functor jk_builder(Ga_, Ka_, Kb_, Da_, Db_, Ca_, Cb_, nalphapi_, nbetapi_);
    process_tei<J_Ka_Kb_Functor>(jk_builder);
    Gb_->copy(Ga_);
    Ga_->subtract(Ka_);
    Gb_->subtract(Kb_);
}

void UHF::save_information()
{
    compute_spin_contamination();
}

void UHF::compute_spin_contamination()
{
    SharedMatrix S = SharedMatrix(factory_->create_matrix("S (Overlap)"));
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(basisset_,basisset_, basisset_,basisset_));
    boost::shared_ptr<OneBodySOInt> so_overlap(fact->so_overlap());
    so_overlap->compute(S);

    double dN = 0.0;

    for (int h =0; h < S->nirrep(); h++) {
        int nbf = S->colspi()[h];
        int nmo = Ca_->colspi()[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];
        if (na == 0 || nb == 0 || nbf == 0 || nmo == 0)
            continue;

        SharedMatrix Ht (new Matrix("H Temp", nbf, nb));
        SharedMatrix Ft (new Matrix("F Temp", na, nb));

        double** Sp = S->pointer(h);
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

bool UHF::test_convergency()
{
    double ediff = E_ - Eold_;

    // RMS of the density
    Matrix Drms;
    Drms.copy(Dt_);
    Drms.subtract(Dtold_);
    Drms_ = 0.5 * Drms.rms();

    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}

void UHF::form_initialF()
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

void UHF::form_F()
{
    Fa_->copy(H_);
    Fa_->add(Ga_);

    Fb_->copy(H_);
    Fb_->add(Gb_);

    if (debug_) {
        Fa_->print(outfile);
        Fb_->print(outfile);
    }
}

void UHF::form_C()
{
    diagonalize_F(Fa_, Ca_, epsilon_a_);
    diagonalize_F(Fb_, Cb_, epsilon_b_);
    find_occupation();
    if (debug_) {
        Ca_->print(outfile);
        Cb_->print(outfile);
    }
}

void UHF::form_D()
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
            ::memset(static_cast<void*>(Da[0]), '\0', sizeof(double)*nso*nso);
        if (nb == 0)
            ::memset(static_cast<void*>(Db[0]), '\0', sizeof(double)*nso*nso);

        C_DGEMM('N','T',nso,nso,na,1.0,Ca[0],nmo,Ca[0],nmo,0.0,Da[0],nso);
        C_DGEMM('N','T',nso,nso,nb,1.0,Cb[0],nmo,Cb[0],nmo,0.0,Db[0],nso);

    }

    Dt_->copy(Da_);
    Dt_->add(Db_);

    if (debug_) {
        fprintf(outfile, "in UHF::form_D:\n");
        Da_->print();
        Db_->print();
    }
}

// TODO: Once Dt_ is refactored to D_ the only difference between this and RHF::compute_initial_E is a factor of 0.5
double UHF::compute_initial_E()
{
    Dt_->copy(Da_);
    Dt_->add(Db_);
    return nuclearrep_ + 0.5 * (Dt_->vector_dot(H_));
}

double UHF::compute_E()
{
    double DH  = Dt_->vector_dot(H_);
    double DFa = Da_->vector_dot(Fa_);
    double DFb = Db_->vector_dot(Fb_);
    double Eelec = 0.5 * (DH + DFa + DFb);
    // fprintf(outfile, "electronic energy = %20.14f\n", Eelec);
    double Etotal = nuclearrep_ + Eelec;
    return Etotal;
}

void UHF::save_fock()
{
    SharedMatrix FDSmSDFa = form_FDSmSDF(Fa_, Da_);
    SharedMatrix FDSmSDFb = form_FDSmSDF(Fb_, Db_);

    if (initialized_diis_manager_ == false) {
        diis_manager_ = boost::shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::OnDisk));
        diis_manager_->set_error_vector_size(2,
                                             DIISEntry::Matrix, FDSmSDFa.get(),
                                             DIISEntry::Matrix, FDSmSDFb.get());
        diis_manager_->set_vector_size(2,
                                       DIISEntry::Matrix, Fa_.get(),
                                       DIISEntry::Matrix, Fb_.get());
        initialized_diis_manager_ = true;
    }

    diis_manager_->add_entry(4, FDSmSDFa.get(), FDSmSDFb.get(), Fa_.get(), Fb_.get());
}

bool UHF::diis()
{
    return diis_manager_->extrapolate(2, Fa_.get(), Fb_.get());
}

}}
