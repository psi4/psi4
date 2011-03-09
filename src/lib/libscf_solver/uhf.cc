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

namespace psi { namespace scf {

UHF::UHF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt) : HF(options, psio, chkpt)
{
    common_init();
}

UHF::UHF(Options& options, shared_ptr<PSIO> psio) : HF(options, psio)
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
    Ca_     = SharedMatrix(factory_->create_matrix("C alpha"));
    Cb_     = SharedMatrix(factory_->create_matrix("C beta"));
    Ga_     = SharedMatrix(factory_->create_matrix("G alpha"));
    Gb_     = SharedMatrix(factory_->create_matrix("G beta"));
    J_      = SharedMatrix(factory_->create_matrix("J total"));
    Ka_     = SharedMatrix(factory_->create_matrix("K alpha"));
    Kb_     = SharedMatrix(factory_->create_matrix("K beta"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = SharedVector(factory_->create_vector());

    fprintf(outfile, "  SCF Algorithm Type is %s.\n", scf_type_.c_str());
    fprintf(outfile, "  DIIS %s.\n", diis_enabled_ ? "enabled" : "disabled");
}

void UHF::finalize()
{
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
    // Print the final docc vector
    char **temp2 = molecule_->irrep_labels();

    if(print_ > 1){
        fprintf(outfile, "\n  Final doubly occupied vector = (");
        for (int h=0; h<factory_->nirrep(); ++h) {
            fprintf(outfile, "%2d %3s ", doccpi_[h], temp2[h]);
        }
        fprintf(outfile, ")\n");
        fprintf(outfile, "  Final singly occupied vector = (");
        for (int h=0; h<factory_->nirrep(); ++h) {
            fprintf(outfile, "%2d %3s ", soccpi_[h], temp2[h]);
        }
        fprintf(outfile, ")\n");
    }

    compute_spin_contamination();

    bool print_mos = options_.get_bool("PRINT_MOS");
    if (print_mos) {
        fprintf(outfile, "\n  Alpha Molecular orbitals:\n");
        Ca_->eivprint(epsilon_a_);

        fprintf(outfile, "\n  Beta Molecular orbitals:\n");
        Cb_->eivprint(epsilon_b_);
    }

    // Print out orbital energies.
    std::vector<std::pair<double, int> > pairsa, pairsb;
    for (int h=0; h<epsilon_a_->nirrep(); ++h) {
        for (int i=0; i<epsilon_a_->dimpi()[h]; ++i) {
            pairsa.push_back(make_pair(epsilon_a_->get(h, i), h));
            pairsb.push_back(make_pair(epsilon_b_->get(h, i), h));
        }
    }
    sort(pairsa.begin(),pairsa.end());
    sort(pairsb.begin(),pairsb.end());
    if(print_){
        fprintf(outfile, "\n  Orbital energies (a.u.):\n    Alpha occupied\n      ");
        for (int i=1; i<=nalpha_; ++i) {
            fprintf(outfile, "%12.6f %3s  ", pairsa[i-1].first, temp2[pairsa[i-1].second]);
            if (i % 4 == 0)
                fprintf(outfile, "\n      ");
        }
        fprintf(outfile, "\n");
        fprintf(outfile, "\n    Alpha unoccupied\n      ");
        for (int i=nalpha_+1; i<=nso(); ++i) {
            fprintf(outfile, "%12.6f %3s  ", pairsa[i-1].first, temp2[pairsa[i-1].second]);
            if ((i-nalpha_) % 4 == 0)
                fprintf(outfile, "\n      ");
        }
        fprintf(outfile, "\n");

        fprintf(outfile, "\n    Beta occupied\n      ");
        for (int i=1; i<=nbeta_; ++i) {
            fprintf(outfile, "%12.6f %3s  ", pairsb[i-1].first, temp2[pairsb[i-1].second]);
            if (i % 4 == 0)
                fprintf(outfile, "\n      ");
        }
        fprintf(outfile, "\n");
        fprintf(outfile, "\n    Beta unoccupied\n      ");
        for (int i=nalpha_+1; i<=nso(); ++i) {
            fprintf(outfile, "%12.6f %3s  ", pairsb[i-1].first, temp2[pairsb[i-1].second]);
            if ((i-nbeta_) % 4 == 0)
                fprintf(outfile, "\n      ");
        }
        fprintf(outfile, "\n");
    }
    for (int i=0; i<epsilon_a_->nirrep(); ++i)
        free(temp2[i]);
    free(temp2);
}

void UHF::compute_spin_contamination()
{
    shared_ptr<Matrix> S = shared_ptr<Matrix>(factory_->create_matrix("S (Overlap)"));
    shared_ptr<IntegralFactory> fact(new IntegralFactory(basisset_,basisset_, basisset_,basisset_));
    shared_ptr<OneBodySOInt> so_overlap(fact->so_overlap());
    so_overlap->compute(S);

    double dN = 0.0;

    for (int h =0; h < S->nirrep(); h++) {
        int nbf = S->colspi()[h];
        int nmo = Ca_->colspi()[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];
        if (na == 0 || nb == 0 || nbf == 0 || nmo == 0)
            continue;

        shared_ptr<Matrix> Ht (new Matrix("H Temp", nbf, nb));
        shared_ptr<Matrix> Ft (new Matrix("F Temp", na, nb));

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

    fprintf(outfile, "\n  @Spin Contamination Metric: %7.4E\n", dS);
      fprintf(outfile, "  @S^2 Expected:              %7.4E\n", S2);
      fprintf(outfile, "  @S^2 Observed:              %7.4E\n", S2 + dS);

}

bool UHF::test_convergency()
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

void UHF::form_initialF()
{
    Fa_->copy(H_);
    Fb_->copy(H_);

    // Transform the Focks
    // Fa_->transform(Shalf_);
    // Fb_->transform(Shalf_);

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
    SharedMatrix eigvec = factory_->create_shared_matrix();

    Fa_->transform(Shalf_);
    Fa_->diagonalize(eigvec, epsilon_a_);

    Ca_->gemm(false, false, 1.0, Shalf_, eigvec, 0.0);

    Fb_->transform(Shalf_);
    Fb_->diagonalize(eigvec, epsilon_b_);

    Cb_->gemm(false, false, 1.0, Shalf_, eigvec, 0.0);

    if (debug_) {
        Ca_->print(outfile);
        Cb_->print(outfile);
    }
}

void UHF::form_D()
{
    int h, i, j, m;
    int *opi = Da_->rowspi();
    int nirreps = Da_->nirrep();
    double val;
    for (h=0; h<nirreps; ++h) {
        for (i=0; i<opi[h]; ++i) {
            for (j=0; j<opi[h]; ++j) {
                val = 0.0;
                for (m=0; m<nalphapi_[h]; ++m)
                    val += Ca_->get(h, i, m) * Ca_->get(h, j, m);
                Da_->set(h, i, j, val);

                val = 0.0;
                for (m=0; m<nbetapi_[h]; ++m)
                    val += Cb_->get(h, i, m) * Cb_->get(h, j, m);
                Db_->set(h, i, j, val);
            }
        }
    }

    // Form total density
    // TODO: Refactor Dt_ to D_ (found in HF)
    Dt_->copy(Da_);
    Dt_->add(Db_);

    if (debug_) {
        Da_->print(outfile);
        Db_->print(outfile);
    }
}

// TODO: Once Dt_ is refactored to D_ the only difference between this and RHF::compute_initial_E is a factor of 0.5
double UHF::compute_initial_E()
{
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
    initialized_diis_manager_ = false;
    if (initialized_diis_manager_ == false) {
        diis_manager_ = shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::OnDisk, psio_));
        diis_manager_->set_error_vector_size(2,
                                             DIISEntry::Matrix, Fa_.get(),
                                             DIISEntry::Matrix, Fb_.get());
        diis_manager_->set_vector_size(2,
                                       DIISEntry::Matrix, Fa_.get(),
                                       DIISEntry::Matrix, Fb_.get());
        initialized_diis_manager_ = true;
    }

    // Determine error matrix for this Fock
    SharedMatrix FaDaS(factory_->create_matrix()), DaS(factory_->create_matrix());
    SharedMatrix SDaFa(factory_->create_matrix()), DaFa(factory_->create_matrix());
    SharedMatrix FbDbS(factory_->create_matrix()), DbS(factory_->create_matrix());
    SharedMatrix SDbFb(factory_->create_matrix()), DbFb(factory_->create_matrix());

    // FDS = F_ * D_ * S_; Alpha
    DaS->gemm(false, false, 1.0, Da_, S_, 0.0);
    FaDaS->gemm(false, false, 1.0, Fa_, DaS, 0.0);
    // SDF = S_ * D_ * F_;
    DaFa->gemm(false, false, 1.0, Da_, Fa_, 0.0);
    SDaFa->gemm(false, false, 1.0, S_, DaFa, 0.0);

    // FDS = F_ * D_ * S_; Beta
    DbS->gemm(false, false, 1.0, Db_, S_, 0.0);
    FbDbS->gemm(false, false, 1.0, Fb_, DbS, 0.0);
    // SDF = S_ * D_ * F_;
    DbFb->gemm(false, false, 1.0, Db_, Fb_, 0.0);
    SDbFb->gemm(false, false, 1.0, S_, DbFb, 0.0);

    Matrix FaDaSmSDaFa;
    FaDaSmSDaFa.copy(FaDaS);
    FaDaSmSDaFa.subtract(SDaFa);
    FaDaSmSDaFa.transform(Shalf_);

    Matrix FbDbSmSDbFb;
    FbDbSmSDbFb.copy(FbDbS);
    FbDbSmSDbFb.subtract(SDbFb);
    FbDbSmSDbFb.transform(Shalf_);

    diis_manager_->add_entry(4, &FaDaSmSDaFa, &FbDbSmSDbFb, Fa_.get(), Fb_.get());
}

bool UHF::diis()
{
    return diis_manager_->extrapolate(2, Fa_.get(), Fb_.get());
}

}}
