#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
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

#include <libmints/view.h>
#include "rohf.h"
#include <psi4-dec.h>

#define _DEBUG

using namespace std;
using namespace psi;
using namespace boost;

namespace psi { namespace scf {

ROHF::ROHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt)
    : HF(options, psio, chkpt)
{
    common_init();
}

ROHF::ROHF(Options& options, boost::shared_ptr<PSIO> psio)
    : HF(options, psio)
{
    common_init();
}

ROHF::~ROHF() {
}

void ROHF::common_init()
{
    Fa_      = SharedMatrix(factory_->create_matrix("Alpha Fock Matrix"));
    Fb_      = SharedMatrix(factory_->create_matrix("Beta Fock Matrix"));
    Feff_    = SharedMatrix(factory_->create_matrix("F effective (MO basis)"));
    soFeff_  = SharedMatrix(factory_->create_matrix("F effective (orthogonalized SO basis)"));
    Ct_      = SharedMatrix(factory_->create_matrix("Orthogonalized Molecular orbitals"));
    Ca_      = SharedMatrix(factory_->create_matrix("Molecular orbitals"));
    Cb_      = Ca_;
    Da_      = SharedMatrix(factory_->create_matrix("Alpha density matrix"));
    Db_      = SharedMatrix(factory_->create_matrix("Beta density matrix"));
    Lagrangian_ = SharedMatrix(factory_->create_matrix("Lagrangian matrix"));
    Ka_      = SharedMatrix(factory_->create_matrix("K alpha"));
    Kb_      = SharedMatrix(factory_->create_matrix("K beta"));
    Ga_      = SharedMatrix(factory_->create_matrix("G alpha"));
    Gb_      = SharedMatrix(factory_->create_matrix("G beta"));
    Dt_old_  = SharedMatrix(factory_->create_matrix("D alpha old"));
    Dt_      = SharedMatrix(factory_->create_matrix("D beta old"));
    moFa_    = SharedMatrix(factory_->create_matrix("MO Basis alpha Fock Matrix"));
    moFb_    = SharedMatrix(factory_->create_matrix("MO Basis beta Fock Matrix"));

    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = epsilon_a_;
    restricted_ = true;
}

void ROHF::semicanonicalize()
{
    // Quick sanity checks
    if(Fa_ == 0 || Fb_ == 0)
        throw PSIEXCEPTION("Wavefunction: Semicanonicalize called, but Fock matrices are not initialized.");
    if(Ca_ == 0 || Cb_ == 0)
        throw PSIEXCEPTION("Wavefunction: Semicanonicalize called, but orbitals are not initialized.");
    if(Ca_ != Cb_)
        throw PSIEXCEPTION("Wavefunction: Semicanonicalize called, but orbitals are not the same.");
    if(Fa_ == Fb_)
        throw PSIEXCEPTION("Wavefunction: Semicanonicalize called, but Fock matrices are the same.");
    if(epsilon_a_ != epsilon_b_)
        throw PSIEXCEPTION("Wavefunction: Semicanonicalize called, but orbital energies are not the same.");
    // Now, make space for the new orbitals
    Cb_ = SharedMatrix(Ca_->clone());
    epsilon_b_ = SharedVector(epsilon_b_->clone());
    Ca_->set_name("Alpha semicanonical orbitals");
    Cb_->set_name("Beta semicanonical orbitals");
    epsilon_a_->set_name("Alpha semicanonical orbital energies");
    epsilon_b_->set_name("Beta semicanonical orbital energies");
    restricted_ = false;

    SharedMatrix Crohf(Ca_->clone());
    SharedMatrix evecs;
    SharedVector evals;

    // Transform the Fock matrix to the MO basis
    SharedMatrix moFa(Fa_->clone());
    SharedMatrix moFb(Fb_->clone());
    moFa->transform(Ca_);
    moFb->transform(Ca_);

    // Pick out occ-occ, and vir-vir subsets of the Fock matrices
    Dimension aoccpi = doccpi_ + soccpi_;
    Dimension boccpi = doccpi_;
    Dimension avirpi = nmopi_ - aoccpi;
    Dimension bvirpi = nmopi_ - boccpi;
    View aOO(moFa, aoccpi, aoccpi);
    View aVV(moFa, avirpi, avirpi, aoccpi, aoccpi);
    View bOO(moFb, boccpi, boccpi);
    View bVV(moFb, bvirpi, bvirpi, boccpi, boccpi);
    SharedMatrix aFOO = aOO();
    SharedMatrix aFVV = aVV();
    SharedMatrix bFOO = bOO();
    SharedMatrix bFVV = bVV();

    // Canonicalize the Alpha occ-occ block
    evecs = SharedMatrix(new Matrix(nirrep_, aoccpi, aoccpi));
    evals = SharedVector(new Vector(nirrep_, aoccpi));
    aFOO->diagonalize(evecs, evals);
    for(int h = 0; h < nirrep_; ++h){
        double **pC  = Crohf->pointer(h);
        double **pCa = Ca_->pointer(h);
        double **pR  = evecs->pointer(h);
        if(aoccpi[h]){
            C_DGEMM('n', 'n', nsopi_[h], aoccpi[h], aoccpi[h], 1.0, pC[0],
                    nmopi_[h], pR[0], aoccpi[h], 0.0, pCa[0], nmopi_[h]);
            for(int p = 0; p < aoccpi[h]; ++p){
                double epsilon = evals->get(h, p);
                epsilon_a_->set(h, p, epsilon);
                for(int q = 0; q < aoccpi[h]; ++q){
                    moFa->set(h, p, q, p==q ? epsilon : 0.0);
                }
            }
        }
    }
    // Canonicalize the Alpha vir-vir block
    evecs = SharedMatrix(new Matrix(nirrep_, avirpi, avirpi));
    evals = SharedVector(new Vector(nirrep_, avirpi));
    aFVV->diagonalize(evecs, evals);
    for(int h = 0; h < nirrep_; ++h){
        double **pC  = Crohf->pointer(h);
        double **pCa = Ca_->pointer(h);
        double **pR  = evecs->pointer(h);
        if(avirpi[h]){
            C_DGEMM('n', 'n', nsopi_[h], avirpi[h], avirpi[h], 1.0, &(pC[0][aoccpi[h]]),
                    nmopi_[h], pR[0], avirpi[h], 0.0, &(pCa[0][aoccpi[h]]), nmopi_[h]);
            for(int p = 0; p < avirpi[h]; ++p){
                double epsilon = evals->get(h, p);
                epsilon_a_->set(h, p+aoccpi[h], epsilon);
                for(int q = 0; q < avirpi[h]; ++q){
                    moFa->set(h, p+aoccpi[h], q+aoccpi[h], p==q ? epsilon : 0.0);
                }
            }
        }
    }
    // Canonicalize the Beta occ-occ block
    evecs = SharedMatrix(new Matrix(nirrep_, boccpi, boccpi));
    evals = SharedVector(new Vector(nirrep_, boccpi));
    bFOO->diagonalize(evecs, evals);
    for(int h = 0; h < nirrep_; ++h){
        double **pC  = Crohf->pointer(h);
        double **pCb = Cb_->pointer(h);
        double **pR  = evecs->pointer(h);
        if(boccpi[h]){
            C_DGEMM('n', 'n', nsopi_[h], boccpi[h], boccpi[h], 1.0, pC[0],
                    nmopi_[h], pR[0], boccpi[h], 0.0, pCb[0], nmopi_[h]);
            for(int p = 0; p < boccpi[h]; ++p){
                double epsilon = evals->get(h, p);
                epsilon_b_->set(h, p, epsilon);
                for(int q = 0; q < boccpi[h]; ++q){
                    moFb->set(h, p, q, p==q ? epsilon : 0.0);
                }
            }
        }
    }
    // Canonicalize the Beta vir-vir block
    evecs = SharedMatrix(new Matrix(nirrep_, bvirpi, bvirpi));
    evals = SharedVector(new Vector(nirrep_, bvirpi));
    bFVV->diagonalize(evecs, evals);
    for(int h = 0; h < nirrep_; ++h){
        double **pC  = Crohf->pointer(h);
        double **pCb = Cb_->pointer(h);
        double **pR  = evecs->pointer(h);
        if(bvirpi[h]){
            C_DGEMM('n', 'n', nsopi_[h], bvirpi[h], bvirpi[h], 1.0, &(pC[0][boccpi[h]]),
                    nmopi_[h], pR[0], bvirpi[h], 0.0, &(pCb[0][boccpi[h]]), nmopi_[h]);
            for(int p = 0; p < bvirpi[h]; ++p){
                double epsilon = evals->get(h, p);
                epsilon_b_->set(h, p+boccpi[h], epsilon);
                for(int q = 0; q < bvirpi[h]; ++q){
                    moFb->set(h, p+boccpi[h], q+boccpi[h], p==q ? epsilon : 0.0);
                }
            }
        }
    }
    // Given the invariance w.r.t. occ-occ rotations, the existing Fa and Fb matrices
    // are still valid, because the densities used to construct them in the old basis
    // are equivalent to those in the new basis.  If the user forward transforms them
    // using the semicanonical orbitals, the correct semicanonical basis Fa and Fb
    // will be obtained.
}

void ROHF::finalize()
{
    // Form Lagrangian
    //
    // In HF, the Lagrangian is Xpi = Fpi. For RHF and UHF, this reduces
    // to Xii = ei. In the AO basis (where we want it), Xmn = Cmi ei Cni.
    // For ROHF, the effective Fock matrix is diagonal, not Fa and Fb (as
    // in UHF). So we need to form the Lagrangian as: Xmn = Cmp Fpi Cni.
    //
    // --EGH
    //
    // Let's build the MO Lagrangian in Feff_
    // ...I'm assuming these bitches are square
    Feff_->zero();
    moFa_->transform(Fa_, Ca_);
    moFb_->transform(Fb_, Ca_);
    for (int h=0; h<nirrep_; ++h) {
        for (int m=0; m<Feff_->rowdim(h); ++m) {
            double tval;
            for (int i=0; i<doccpi_[h]; ++i) {
                tval = moFa_->get(h, m, i);
                tval += moFb_->get(h, m, i);
                Feff_->set(h, m, i, tval);
            }
            for (int i=doccpi_[h]; i<doccpi_[h]+soccpi_[h]; ++i) {
                tval = moFa_->get(h, m, i);
                Feff_->set(h, m, i, tval);
            }
        }
    }
    Lagrangian_->back_transform(Feff_, Ca_);

    Feff_.reset();
    Ka_.reset();
    Kb_.reset();
    Ga_.reset();
    Gb_.reset();
    Dt_old_.reset();
    Dt_.reset();
    moFa_.reset();
    moFb_.reset();

    HF::finalize();
}

void ROHF::save_density_and_energy()
{
    Dt_old_->copy(Dt_);
    Eold_ = E_; // save previous energy
}

void ROHF::save_information()
{
}

void ROHF::save_fock()
{
    if (initialized_diis_manager_ == false) {
        diis_manager_ = boost::shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::OnDisk));
        diis_manager_->set_error_vector_size(1, DIISEntry::Matrix, soFeff_.get());
        diis_manager_->set_vector_size(1, DIISEntry::Matrix, soFeff_.get());
        initialized_diis_manager_ = true;
    }

    SharedMatrix errvec(Feff_);
    errvec->zero_diagonal();
    errvec->back_transform(Ct_);
    diis_manager_->add_entry(2, errvec.get(), soFeff_.get());
}

bool ROHF::diis()
{
    return diis_manager_->extrapolate(1, soFeff_.get());
}

bool ROHF::test_convergency()
{
    // energy difference
    double ediff = E_ - Eold_;

    // RMS of the density
    Matrix D_rms;
    D_rms.copy(Dt_);
    D_rms.subtract(Dt_old_);
    Drms_ = D_rms.rms();

    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_)
        return true;
    else
        return false;
}

void ROHF::form_initialF()
{
    // Form the initial Fock matrix, closed and open variants
    Fa_->copy(H_);
    Fa_->transform(X_);
    Fb_->copy(Fa_);

#ifdef _DEBUG
    if (debug_) {
        fprintf(outfile, "Initial alpha Fock matrix:\n");
        Fa_->print(outfile);
        fprintf(outfile, "Initial beta Fock matrix:\n");
        Fb_->print(outfile);
    }
#endif
}

void ROHF::form_F()
{
    // Start by constructing the standard Fa and Fb matrices encountered in UHF
    Fa_->copy(H_);
    Fb_->copy(H_);
    Fa_->add(Ga_);
    Fb_->add(Gb_);

    moFa_->transform(Fa_, Ca_);
    moFb_->transform(Fb_, Ca_);

    /*
     * Fo = open-shell fock matrix = 0.5 Fa
     * Fc = closed-shell fock matrix = 0.5 (Fa + Fb)
     *
     * Therefore
     *
     * 2(Fc-Fo) = Fb
     * 2Fo = Fa
     *
     * Form the effective Fock matrix, too
     * The effective Fock matrix has the following structure
     *          |  closed     open    virtual
     *  ----------------------------------------
     *  closed  |    Fc     2(Fc-Fo)    Fc
     *  open    | 2(Fc-Fo)     Fc      2Fo
     *  virtual |    Fc       2Fo       Fc
     */
    Feff_->copy(moFa_);
    Feff_->add(moFb_);
    Feff_->scale(0.5);
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = doccpi_[h]; i < doccpi_[h] + soccpi_[h]; ++i) {
            // Set the open/closed portion
            for (int j = 0; j < doccpi_[h]; ++j) {
                double val = moFb_->get(h, i, j);
                Feff_->set(h, i, j, val);
                Feff_->set(h, j, i, val);
            }
            // Set the open/virtual portion
            for (int j = doccpi_[h] + soccpi_[h]; j < nmopi_[h]; ++j) {
                double val = moFa_->get(h, i, j);
                Feff_->set(h, i, j, val);
                Feff_->set(h, j, i, val);
            }
        }
    }

    // Form the orthogonalized SO basis Feff matrix, for use in DIIS
    soFeff_->copy(Feff_);
    soFeff_->back_transform(Ct_);

    if (debug_) {
        Fa_->print();
        Fb_->print();
        moFa_->print();
        moFb_->print();
        Feff_->print();
        soFeff_->print();
    }
}

void ROHF::form_C()
{
    soFeff_->diagonalize(Ct_, epsilon_a_);
    //Form C = XC'
    Ca_->gemm(false, false, 1.0, X_, Ct_, 0.0);

    find_occupation();

    if (debug_) {
        Ca_->print(outfile);
        fprintf(outfile, "In ROHF::form_C:\n");
        Ct_->eivprint(epsilon_a_);
    }
}

void ROHF::form_initial_C()
{
    // In ROHF the creation of the C matrix depends on the previous iteration's C
    // matrix. Here we use Fa to generate the first C, where Fa was set by guess()
    // to either H or the GWH Hamiltonian.
    Fa_->transform(X_);
    Fa_->diagonalize(Ct_, epsilon_a_);
    find_occupation();
    Ca_->gemm(false, false, 1.0, X_, Ct_, 0.0);

    if (print_ > 3)
        Ca_->print(outfile, "initial C");
}

void ROHF::form_D()
{
    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];

        if (nso == 0 || nmo == 0) continue;

        double** Ca = Ca_->pointer(h);
        double** Da = Da_->pointer(h);
        double** Db = Db_->pointer(h);

        if (na == 0)
            memset(static_cast<void*>(Da[0]), '\0', sizeof(double)*nso*nso);
        if (nb == 0)
            memset(static_cast<void*>(Db[0]), '\0', sizeof(double)*nso*nso);


        C_DGEMM('N','T',nso,nso,na,1.0,Ca[0],nmo,Ca[0],nmo,0.0,Da[0],nso);
        C_DGEMM('N','T',nso,nso,nb,1.0,Ca[0],nmo,Ca[0],nmo,0.0,Db[0],nso);

    }

    Dt_->copy(Da_);
    Dt_->add(Db_);

    if (debug_) {
        fprintf(outfile, "in ROHF::form_D:\n");
        Da_->print();
        Db_->print();
    }
}

double ROHF::compute_initial_E()
{
    return 0.5 * (compute_E() + nuclearrep_);
}

double ROHF::compute_E() 
{
    double one_electron_E = Da_->vector_dot(H_) + Db_->vector_dot(H_);;
    double two_electron_E = 0.5 * (Da_->vector_dot(Fa_) + Db_->vector_dot(Fb_) - one_electron_E);   
 
    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] = two_electron_E; 
    energies_["XC"] = 0.0;
    energies_["-D"] = 0.0;

    double DH  = Da_->vector_dot(H_);
    DH += Db_->vector_dot(H_);
    double DFa = Da_->vector_dot(Fa_);
    double DFb = Db_->vector_dot(Fb_);
    double Eelec = 0.5 * (DH + DFa + DFb);
    double Etotal = nuclearrep_ + Eelec;
    return Etotal;
}

void ROHF::form_G()
{
    /*
     * This just builds the same Ga and Gb matrices used in UHF
     */
    if (scf_type_ == "DF" || scf_type_ == "PS") {

        std::vector<SharedMatrix> & C = jk_->C_left();
        C.clear();
        C.push_back(Ca_subset("SO", "OCC"));
        C.push_back(Cb_subset("SO", "OCC"));
        
        // Run the JK object
        jk_->compute();

        // Pull the J and K matrices off
        const std::vector<SharedMatrix> & J = jk_->J();
        const std::vector<SharedMatrix> & K = jk_->K();
        Ga_->copy(J[0]);
        Ga_->add(J[1]);
        Ka_ = K[0];
        Kb_ = K[1];
        
    } else {
        J_Ka_Kb_Functor jk_builder(Ga_, Ka_, Kb_, Da_, Db_, Ca_, Cb_, nalphapi_, nbetapi_);
        process_tei<J_Ka_Kb_Functor>(jk_builder);
    }

    Gb_->copy(Ga_);
    Ga_->subtract(Ka_);
    Gb_->subtract(Kb_);
}

void ROHF::stability_analysis()
{
    if(scf_type_ == "DF"){
        throw PSIEXCEPTION("Stability analysis has not been implemented for density fitted wavefunctions yet.");
    }else{
        // Build the Fock Matrix
        SharedMatrix aMoF(new Matrix("Alpha MO basis fock matrix", nmopi_, nmopi_));
        SharedMatrix bMoF(new Matrix("Beta MO basis fock matrix", nmopi_, nmopi_));
        aMoF->transform(Fa_, Ca_);
        bMoF->transform(Fb_, Ca_);

        std::vector<boost::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::occ);
        spaces.push_back(MOSpace::vir);
        // Ref wfn is really "this"
        boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
#define ID(x) ints.DPD_ID(x)
        IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly,
                               IntegralTransform::QTOrder, IntegralTransform::None);
        ints.set_keep_dpd_so_ints(true);
        ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
        ints.transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
        dpd_set_default(ints.get_dpd_id());
        dpdbuf4 Aaa, Aab, Aba, Abb, I, A;
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
        // A_IA_JB = (IA|JB)
        dpd_buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "ROHF Hessian (IA|JB)", 1.0);
        // A_IA_jb = (IA|jb)
        dpd_buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "ROHF Hessian (IA|jb)", 1.0);

        // A_IA_JB -= 0.5 (IB|JA)
        dpd_buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, psrq,
                           ID("[O,V]"), ID("[O,V]"), "ROHF Hessian (IA|JB)", -0.5);
        dpd_buf4_close(&I);

        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
        // A_IA_JB -= 0.5 (IJ|AB)
        dpd_buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, prqs,
                           ID("[O,V]"), ID("[O,V]"), "ROHF Hessian (IA|JB)", -0.5);
        dpd_buf4_close(&I);

        dpd_buf4_init(&Aaa, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "ROHF Hessian (IA|JB)");

        // A_ia_jb = A_IA_JB
        dpd_buf4_copy(&Aaa, PSIF_LIBTRANS_DPD, "ROHF Hessian (ia|jb)");
        for(int h = 0; h < Aaa.params->nirreps; ++h){
            dpd_buf4_mat_irrep_init(&Aaa, h);
            dpd_buf4_mat_irrep_rd(&Aaa, h);
            for(int ia = 0; ia < Aaa.params->rowtot[h]; ++ia){
                int iabs = Aaa.params->roworb[h][ia][0];
                int aabs = Aaa.params->roworb[h][ia][1];
                int isym = Aaa.params->psym[iabs];
                int asym = Aaa.params->qsym[aabs];
                int irel = iabs - Aaa.params->poff[isym];
                int arel = aabs - Aaa.params->qoff[asym];
                for(int jb = 0; jb < Aaa.params->coltot[h]; ++jb){
                    int jabs = Aaa.params->colorb[h][jb][0];
                    int babs = Aaa.params->colorb[h][jb][1];
                    int jsym = Aaa.params->rsym[jabs];
                    int bsym = Aaa.params->ssym[babs];
                    int jrel = jabs - Aaa.params->roff[jsym];
                    int brel = babs - Aaa.params->soff[bsym];
                    double val = Aaa.matrix[h][ia][jb];
                    // A_IA_JB += 0.5 delta_IJ F_AB - 0.5 delta_AB F_IJ
                    if((iabs == jabs) && (asym == bsym))
                        val += 0.5 * aMoF->get(asym, arel + doccpi_[asym], brel + doccpi_[bsym]);
                    if((aabs == babs) && (isym == jsym))
                        val -= 0.5 * aMoF->get(isym, irel, jrel);
                    // Zero out any socc-socc terms
                    if(arel < (soccpi_[asym]) || brel < (soccpi_[bsym]))
                        val = 0.0;
                    Aaa.matrix[h][ia][jb] = val;
                }
            }
            dpd_buf4_mat_irrep_wrt(&Aaa, h);
        }
        dpd_buf4_close(&Aaa);

        dpd_buf4_init(&Abb, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "ROHF Hessian (ia|jb)");
        for(int h = 0; h < Abb.params->nirreps; ++h){
            dpd_buf4_mat_irrep_init(&Abb, h);
            dpd_buf4_mat_irrep_rd(&Abb, h);
            for(int ia = 0; ia < Abb.params->rowtot[h]; ++ia){
                int iabs = Abb.params->roworb[h][ia][0];
                int aabs = Abb.params->roworb[h][ia][1];
                int isym = Abb.params->psym[iabs];
                int asym = Abb.params->qsym[aabs];
                int irel = iabs - Abb.params->poff[isym];
                int arel = aabs - Abb.params->qoff[asym];
                for(int jb = 0; jb < Abb.params->coltot[h]; ++jb){
                    int jabs = Abb.params->colorb[h][jb][0];
                    int babs = Abb.params->colorb[h][jb][1];
                    int jsym = Abb.params->rsym[jabs];
                    int bsym = Abb.params->ssym[babs];
                    int jrel = jabs - Abb.params->roff[jsym];
                    int brel = babs - Abb.params->soff[bsym];
                    double val = Abb.matrix[h][ia][jb];
                    // A_ia_jb += 0.5 delta_ij F_ab - 0.5 delta_ab F_ij
                    if((iabs == jabs) && (asym == bsym))
                        val += 0.5 * bMoF->get(asym, arel + doccpi_[asym], brel + doccpi_[bsym]);
                    if((aabs == babs) && (isym == jsym))
                        val -= 0.5 * bMoF->get(isym, irel, jrel);
                    // Zero out any socc-socc terms
                    if(irel >= (doccpi_[isym]) || jrel >= (doccpi_[jsym]))
                        val = 0.0;
                    Abb.matrix[h][ia][jb] = val;
                }
            }
            dpd_buf4_mat_irrep_wrt(&Abb, h);
        }
        dpd_buf4_close(&Abb);

        dpd_buf4_init(&Aab, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "ROHF Hessian (IA|jb)");
        for(int h = 0; h < Aab.params->nirreps; ++h){
            dpd_buf4_mat_irrep_init(&Aab, h);
            dpd_buf4_mat_irrep_rd(&Aab, h);
            for(int ia = 0; ia < Aab.params->rowtot[h]; ++ia){
                int iabs = Aab.params->roworb[h][ia][0];
                int aabs = Aab.params->roworb[h][ia][1];
                int isym = Aab.params->psym[iabs];
                int asym = Aab.params->qsym[aabs];
                int irel = iabs - Aab.params->poff[isym];
                int arel = aabs - Aab.params->qoff[asym];
                for(int jb = 0; jb < Aab.params->coltot[h]; ++jb){
                    int jabs = Aab.params->colorb[h][jb][0];
                    int babs = Aab.params->colorb[h][jb][1];
                    int jsym = Aab.params->rsym[jabs];
                    int bsym = Aab.params->ssym[babs];
                    int jrel = jabs - Aab.params->roff[jsym];
                    int brel = babs - Aab.params->soff[bsym];
                    // A_IA_jb += 0.5 delta_Ib F(beta)_jA
                    // Don't forget to account for the 0-based numbering of b
                    if((irel == (brel + doccpi_[bsym])) && (jsym == asym))
                        Aab.matrix[h][ia][jb] += 0.5 * bMoF->get(jsym, jrel, arel + doccpi_[asym]);
                    // Zero out any socc-socc terms
                    if(jrel >= (doccpi_[jsym]) || arel < (soccpi_[asym]))
                        Aab.matrix[h][ia][jb] = 0.0;
                }
            }
            dpd_buf4_mat_irrep_wrt(&Aab, h);
        }

        // A_ia_JB = A_IA_jb
        dpd_buf4_sort(&Aab, PSIF_LIBTRANS_DPD, rspq, ID("[O,V]"), ID("[O,V]"), "ROHF Hessian (ia|JB)");
        dpd_buf4_close(&Aab);

        // Alpha-Alpha
        dpd_buf4_init(&Aaa, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "ROHF Hessian (IA|JB)");
        dpd_buf4_copy(&Aaa, PSIF_LIBTRANS_DPD, "ROHF Hessian");
        dpd_buf4_close(&Aaa);
        dpd_buf4_init(&A, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "ROHF Hessian");
        // Alpha-Beta
        dpd_buf4_init(&Aab, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "ROHF Hessian (IA|jb)");
        dpd_buf4_axpy(&Aab, &A, 1.0);
        dpd_buf4_close(&Aab);
        // Beta-Alpha
        dpd_buf4_init(&Aba, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "ROHF Hessian (ia|JB)");
        dpd_buf4_axpy(&Aba, &A, 1.0);
        dpd_buf4_close(&Aba);
        // Beta-Beta
        dpd_buf4_init(&Abb, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "ROHF Hessian (ia|jb)");
        dpd_buf4_axpy(&Abb, &A, 1.0);
        dpd_buf4_close(&Abb);
        dpd_buf4_close(&A);

        /*
         *  Perform the stability analysis
         */
        std::vector<std::pair<double, int> >eval_sym;
        dpd_buf4_init(&A, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "ROHF Hessian");
        for(int h = 0; h < A.params->nirreps; ++h) {
            int npairs = A.params->rowtot[h];
            if(npairs == 0) continue;

            // Store the row indices, for convenience
            unsigned int rank = 0;
            double **U = block_matrix(npairs, npairs);
            for(int ia = 0; ia < npairs; ++ia){
                int iabs = Aab.params->roworb[h][ia][0];
                int aabs = Aab.params->roworb[h][ia][1];
                int isym = Aab.params->psym[iabs];
                int asym = Aab.params->qsym[aabs];
                int irel = iabs - Aab.params->poff[isym];
                int arel = aabs - Aab.params->qoff[asym];
                if(arel >= soccpi_[asym] || irel < doccpi_[isym]){
                    U[ia][ia] = 1.0;
                    rank++;
                }
            }
            int lastcol = npairs - 1;
            for(int ia = 0; ia < npairs; ++ia){
                if(U[ia][ia] == 0.0){
                    while(U[lastcol][lastcol] == 0.0) lastcol--;
                    if(lastcol > ia){
                        U[lastcol][ia] = U[ia][lastcol] = 1.0;
                        U[lastcol][lastcol] = 0.0;
                    }
                }
            }
            if(rank == 0) continue;

            dpd_buf4_mat_irrep_init(&A, h);
            dpd_buf4_mat_irrep_rd(&A, h);

            // Use the transformation matrix to rearrange the columns
            double **temp = block_matrix(npairs, npairs);
            C_DGEMM('n', 'n', npairs, npairs, npairs, 1.0, A.matrix[h][0], npairs,
                    U[0], npairs, 0.0, temp[0], npairs);
            C_DGEMM('n', 'n', npairs, npairs, npairs, 1.0, U[0], npairs,
                    temp[0], npairs, 0.0, A.matrix[h][0], npairs);

            double *evals = init_array(rank);
            double **evecs = block_matrix(rank, rank);

            sq_rsp(rank, rank, A.matrix[h], evals, 1, evecs, 1e-12);
            dpd_buf4_mat_irrep_close(&A, h);
            int mindim = rank < 15 ? rank : 15;
            for(int i = 0; i < mindim; i++)
                eval_sym.push_back(std::make_pair(evals[i], h));

            free_block(evecs);
            delete [] evals;
        }
        fprintf(outfile, "\tLowest ROHF->ROHF stability eigenvalues:-\n");
        print_stability_analysis(eval_sym);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }
}


}}
