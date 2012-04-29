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
#include "libtrans/integraltransform.h"
#include "libdpd/dpd.h"

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
        
        Ga_->copy(J_);

    } else {
        J_Ka_Kb_Functor jk_builder(Ga_, Ka_, Kb_, Da_, Db_, Ca_, Cb_, nalphapi_, nbetapi_);
        process_tei<J_Ka_Kb_Functor>(jk_builder);
    }

    Gb_->copy(Ga_);
    Ga_->subtract(Ka_);
    Gb_->subtract(Kb_);
}

void UHF::save_information()
{
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

void UHF::stability_analysis()
{
    if(scf_type_ == "DF"){
        throw PSIEXCEPTION("Stability analysis has not been implemented for density fitted wavefunctions yet.");
    }else{
        // Build the Fock Matrix
        SharedMatrix aMoF(new Matrix("Alpha MO basis fock matrix", nmopi_, nmopi_));
        SharedMatrix bMoF(new Matrix("Beta MO basis fock matrix", nmopi_, nmopi_));
        aMoF->transform(Fa_, Ca_);
        bMoF->transform(Fb_, Cb_);

        std::vector<boost::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::occ);
        spaces.push_back(MOSpace::vir);
        // Ref wfn is really "this"
        boost::shared_ptr<Wavefunction> wfn = Process::environment.reference_wavefunction();
#define ID(x) ints.DPD_ID(x)
        IntegralTransform ints(wfn, spaces, IntegralTransform::Unrestricted, IntegralTransform::DPDOnly,
                               IntegralTransform::QTOrder, IntegralTransform::None);
        ints.set_keep_dpd_so_ints(true);
        ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
        ints.transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
        dpd_set_default(ints.get_dpd_id());
        dpdbuf4 Aaa, Aab, Abb, I;
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                      ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
        // A_IA_jb = 2 (IA|jb)
        dpd_buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "UHF Hessian (IA|jb)", 2.0);
        dpd_buf4_close(&I);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
        // A_IA_JB = 2 (IA|JB)
        dpd_buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "UHF Hessian (IA|JB)", 2.0);
        // A_IA_JB -= (IB|JA)
        dpd_buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, psrq,
                           ID("[O,V]"), ID("[O,V]"), "UHF Hessian (IA|JB)", -1.0);
        dpd_buf4_close(&I);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                      ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
        // A_ia_jb = 2 (ia|jb)
        dpd_buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "UHF Hessian (ia|jb)", 2.0);
        // A_ia_jb -= (ib|ja)
        dpd_buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, psrq,
                           ID("[o,v]"), ID("[o,v]"), "UHF Hessian (ia|jb)", -1.0);
        dpd_buf4_close(&I);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
        // A_IA_JB -= (IJ|AB)
        dpd_buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, prqs,
                           ID("[O,V]"), ID("[O,V]"), "UHF Hessian (IA|JB)", -1.0);
        dpd_buf4_close(&I);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
        // A_ia_jb -= (ij|ab)
        dpd_buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, prqs,
                           ID("[o,v]"), ID("[o,v]"), "UHF Hessian (ia|jb)", -1.0);
        dpd_buf4_close(&I);

        dpd_buf4_init(&Aaa, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "UHF Hessian (IA|JB)");
        for(int h = 0; h < Aaa.params->nirreps; ++h){
            dpd_buf4_mat_irrep_init(&Aaa, h);
            dpd_buf4_mat_irrep_rd(&Aaa, h);
            for(int ia = 0; ia < Aaa.params->rowtot[h]; ++ia){
                int iabs = Aaa.params->roworb[h][ia][0];
                int aabs = Aaa.params->roworb[h][ia][1];
                int isym = Aaa.params->psym[iabs];
                int asym = Aaa.params->qsym[aabs];
                int irel = iabs - Aaa.params->poff[isym];
                int arel = aabs - Aaa.params->qoff[asym] + soccpi_[asym] + doccpi_[asym];
                for(int jb = 0; jb < Aaa.params->coltot[h]; ++jb){
                    int jabs = Aaa.params->colorb[h][jb][0];
                    int babs = Aaa.params->colorb[h][jb][1];
                    int jsym = Aaa.params->rsym[jabs];
                    int bsym = Aaa.params->ssym[babs];
                    int jrel = jabs - Aaa.params->roff[jsym];
                    int brel = babs - Aaa.params->soff[bsym] + soccpi_[asym] + doccpi_[asym];
                    // A_IA_JB += delta_IJ F_AB - delta_AB F_IJ
                    if((iabs == jabs) && (asym == bsym))
                        Aaa.matrix[h][ia][jb] += aMoF->get(asym, arel, brel);
                    if((aabs == babs) && (isym == jsym))
                        Aaa.matrix[h][ia][jb] -= aMoF->get(isym, irel, jrel);
                }
            }
            dpd_buf4_mat_irrep_wrt(&Aaa, h);
        }
        dpd_buf4_close(&Aaa);

        dpd_buf4_init(&Abb, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                      ID("[o,v]"), ID("[o,v]"), 0, "UHF Hessian (ia|jb)");

        for(int h = 0; h < Abb.params->nirreps; ++h){
            dpd_buf4_mat_irrep_init(&Abb, h);
            dpd_buf4_mat_irrep_rd(&Abb, h);
            for(int ia = 0; ia < Abb.params->rowtot[h]; ++ia){
                int iabs = Abb.params->roworb[h][ia][0];
                int aabs = Abb.params->roworb[h][ia][1];
                int isym = Abb.params->psym[iabs];
                int asym = Abb.params->qsym[aabs];
                int irel = iabs - Abb.params->poff[isym];
                int arel = aabs - Abb.params->qoff[asym] + doccpi_[asym];
                for(int jb = 0; jb < Abb.params->coltot[h]; ++jb){
                    int jabs = Abb.params->colorb[h][jb][0];
                    int babs = Abb.params->colorb[h][jb][1];
                    int jsym = Abb.params->rsym[jabs];
                    int bsym = Abb.params->ssym[babs];
                    int jrel = jabs - Abb.params->roff[jsym];
                    int brel = babs - Abb.params->soff[bsym] + doccpi_[asym];
                    // A_ia_jb += delta_ij F_ab - delta_ab F_ij
                    if((iabs == jabs) && (asym == bsym))
                        Abb.matrix[h][ia][jb] += bMoF->get(asym, arel, brel);
                    if((aabs == babs) && (isym == jsym))
                        Abb.matrix[h][ia][jb] -= bMoF->get(isym, irel, jrel);
                }
            }
            dpd_buf4_mat_irrep_wrt(&Abb, h);
        }
        dpd_buf4_close(&Abb);

        /*
         *  Perform the stability analysis
         */
        std::vector<std::pair<double, int> >eval_sym;

        dpd_buf4_init(&Aaa, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "UHF Hessian (IA|JB)");
        dpd_buf4_init(&Aab, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                      ID("[O,V]"), ID("[o,v]"), 0, "UHF Hessian (IA|jb)");
        dpd_buf4_init(&Abb, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                      ID("[o,v]"), ID("[o,v]"), 0, "UHF Hessian (ia|jb)");
        for(int h = 0; h < Aaa.params->nirreps; ++h) {
            int aDim = Aaa.params->rowtot[h];
            int bDim = Abb.params->rowtot[h];
            int dim =  aDim + bDim;
            if(dim == 0) continue;
            double *evals = init_array(dim);
            double **evecs = block_matrix(dim, dim);
            double **A = block_matrix(dim, dim);

            // Alpha-alpha contribution to the Hessian
            dpd_buf4_mat_irrep_init(&Aaa, h);
            dpd_buf4_mat_irrep_rd(&Aaa, h);
            for(int ia = 0; ia < aDim; ++ia)
                for(int jb = 0; jb < aDim; ++jb)
                    A[ia][jb] = Aaa.matrix[h][ia][jb];
            dpd_buf4_mat_irrep_close(&Aaa, h);

            // Alpha-beta and beta-alpha contribution to the Hessian
            dpd_buf4_mat_irrep_init(&Aab, h);
            dpd_buf4_mat_irrep_rd(&Aab, h);
            for(int ia = 0; ia < aDim; ++ia)
                for(int jb = 0; jb < bDim; ++jb)
                    A[ia][jb + aDim] = A[jb + aDim][ia] = Aab.matrix[h][ia][jb];
            dpd_buf4_mat_irrep_close(&Aab, h);
            // Beta-beta contribution to the Hessian
            dpd_buf4_mat_irrep_init(&Abb, h);
            dpd_buf4_mat_irrep_rd(&Abb, h);
            for(int ia = 0; ia < bDim; ++ia)
                for(int jb = 0; jb < bDim; ++jb)
                    A[ia + aDim][jb + aDim] = Abb.matrix[h][ia][jb];
            dpd_buf4_mat_irrep_close(&Abb, h);

            sq_rsp(dim, dim, A, evals, 1, evecs, 1e-12);

            int mindim = dim < 5 ? dim : 5;
            for(int i = 0; i < mindim; i++)
                eval_sym.push_back(std::make_pair(evals[i], h));

            free_block(A);
            free_block(evecs);
            delete [] evals;
        }

        fprintf(outfile, "\tLowest UHF->UHF stability eigenvalues:-\n");
        print_stability_analysis(eval_sym);

        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }
}


}}
