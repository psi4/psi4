/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/physconst.h"

#include "psi4/libmints/basisset_parser.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/factory.h"
#include "psi4/libfock/v.h"
#include "psi4/libfock/jk.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libdpd/dpd.h"
#include "rhf.h"

using namespace std;

namespace psi { namespace scf {

RHF::RHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func)
    : HF(ref_wfn, func, Process::environment.options, PSIO::shared_object())
{
    common_init();
}

RHF::RHF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func,
         Options& options, std::shared_ptr<PSIO> psio)
    : HF(ref_wfn, func, options, psio)
{
    common_init();
}

RHF::~RHF()
{
}

void RHF::common_init()
{
    if (multiplicity_ != 1) throw PSIEXCEPTION("RHF: RHF reference is only for singlets.");
    Drms_ = 0.0;

    // Allocate matrix memory
    Fa_        = SharedMatrix(factory_->create_matrix("F"));
    Fb_        = Fa_;
    Ca_        = SharedMatrix(factory_->create_matrix("C"));
    Cb_        = Ca_;
    epsilon_a_ = SharedVector(factory_->create_vector());
    epsilon_b_ = epsilon_a_;
    Da_        = SharedMatrix(factory_->create_matrix("SCF density"));
    Db_        = Da_;
    Lagrangian_ = SharedMatrix(factory_->create_matrix("X"));
    D_         = Da_;
    Dold_      = SharedMatrix(factory_->create_matrix("D old"));
    Va_        = SharedMatrix(factory_->create_matrix("V"));
    Vb_        = Va_;
    G_         = SharedMatrix(factory_->create_matrix("G"));
    J_         = SharedMatrix(factory_->create_matrix("J"));
    K_         = SharedMatrix(factory_->create_matrix("K"));
    wK_        = SharedMatrix(factory_->create_matrix("wK"));

    same_a_b_dens_ = true;
    same_a_b_orbs_ = true;
}

void RHF::finalize()
{
    // Form lagrangian
    for (int h=0; h<nirrep_; ++h) {
        for (int m=0; m<Lagrangian_->rowdim(h); ++m) {
            for (int n=0; n<Lagrangian_->coldim(h); ++n) {
                double sum = 0.0;
                for (int i=0; i<doccpi_[h]; ++i) {
                    sum += epsilon_a_->get(h, i) * Ca_->get(h, m, i) * Ca_->get(h, n, i);
                }
                Lagrangian_->set(h, m, n, sum);
            }
        }
    }

    Dold_.reset();
    G_.reset();
    J_.reset();
    K_.reset();
    wK_.reset();

    HF::finalize();
}

SharedMatrix RHF::Da() const
{
    return D_;
}

void RHF::save_density_and_energy()
{
    Dold_->copy(D_);  // Save previous density
    Eold_ = E_;       // Save previous energy
}

void forPermutation(int depth, vector<int>& array,
      vector<int>& indices,int curDepth,vector<vector<int> >& finalindex) {
   int length=array.size();
   if(curDepth == 0) {
        finalindex.push_back(indices);
        return;
    }
    for (int i = 0; i < length; i++) {
       bool isgood=true;
       for(int j=length-1;j>=curDepth&&isgood;j--){
          if(indices[j]==array[i])isgood=false;
       }
       if(isgood){
          indices[curDepth-1]= array[i];
          forPermutation(depth, array,indices, curDepth - 1,finalindex);
       }
    }
}
void RHF::form_V()
{
    // Push the C matrix on
    std::vector<SharedMatrix> & C = potential_->C();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));

    // Run the potential object
    potential_->compute();

    // Pull the V matrices off
    const std::vector<SharedMatrix> & V = potential_->V();
    Va_ = V[0];
    Vb_ = Va_;
}
void RHF::form_G()
{
    if (functional_->needs_xc()) {
        timer_on("RKS: Form V");
        form_V();
        G_->copy(Va_);
        timer_off("RKS: Form V");
    } else {
        G_->zero();
    }

    /// Push the C matrix on
    std::vector<SharedMatrix> & C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));

    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedMatrix> & J = jk_->J();
    const std::vector<SharedMatrix> & K = jk_->K();
    const std::vector<SharedMatrix> & wK = jk_->wK();
    J_ = J[0];
    if (functional_->is_x_hybrid()) {
        K_ = K[0];
    }
    if (functional_->is_x_lrc()) {
        wK_ = wK[0];
    }

    G_->axpy(2.0, J_);

    double alpha = functional_->x_alpha();
    double beta = 1.0 - alpha;

    if (alpha != 0.0) {
        G_->axpy(-alpha, K_);
    } else {
        K_->zero();
    }

    if (functional_->is_x_lrc()) {
        G_->axpy(-beta, wK_);
    } else {
        wK_->zero();
    }
}

void RHF::save_information()
{
}

void RHF::compute_orbital_gradient(bool save_fock)
{
    // Conventional DIIS (X'[FDS - SDF]X, where X levels things out)
    SharedMatrix gradient = form_FDSmSDF(Fa_, Da_);
    Drms_ = gradient->rms();

    if(save_fock){
        if (initialized_diis_manager_ == false) {
            if (scf_type_ == "direct")
                diis_manager_ = std::shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::InCore));
            else
                diis_manager_ = std::shared_ptr<DIISManager>(new DIISManager(max_diis_vectors_, "HF DIIS vector", DIISManager::LargestError, DIISManager::OnDisk));
            diis_manager_->set_error_vector_size(1, DIISEntry::Matrix, gradient.get());
            diis_manager_->set_vector_size(1, DIISEntry::Matrix, Fa_.get());
            initialized_diis_manager_ = true;
        }
        diis_manager_->add_entry(2, gradient.get(), Fa_.get());
    }
}


bool RHF::diis()
{
    return diis_manager_->extrapolate(1, Fa_.get());
}

bool RHF::test_convergency()
{
    // energy difference
    double ediff = E_ - Eold_;

    // Drms was computed earlier
    if (fabs(ediff) < energy_threshold_ && Drms_ < density_threshold_){
       return true;
    }
    else
        return false;
}


void RHF::form_F()
{
    Fa_->copy(H_);
    Fa_->add(G_);

    if (debug_) {
        Fa_->print();
        J_->print();
        K_->print();
        if (functional_->needs_xc()){
            Va_->print();
        }
        G_->print();

    }
}

void RHF::form_C()
{
    diagonalize_F(Fa_, Ca_, epsilon_a_);
    find_occupation();
}

void RHF::form_D()
{
    for (int h = 0; h < nirrep_; ++h) {
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int na = doccpi_[h];

        if (nso == 0 || nmo == 0) continue;

        double** Ca = Ca_->pointer(h);
        double** D = D_->pointer(h);

        if (na == 0)
            memset(static_cast<void*>(D[0]), '\0', sizeof(double)*nso*nso);

        C_DGEMM('N','T',nso,nso,na,1.0,Ca[0],nmo,Ca[0],nmo,0.0,D[0],nso);

    }

    if (debug_) {
        outfile->Printf( "in RHF::form_D:\n");
        D_->print();
    }
}

void RHF::damp_update()
{
  D_->scale(1.0 - damping_percentage_);
  D_->axpy(damping_percentage_, Dold_);
}

double RHF::compute_initial_E()
{
    double Etotal = nuclearrep_ + D_->vector_dot(H_);
    return Etotal;
}

double RHF::compute_E()
{
    double one_electron_E = 2.0 * D_->vector_dot(H_);
    double coulomb_E = 2.0 * D_->vector_dot(J_);

    double XC_E = 0.0;
    if (functional_->needs_xc()) {
        XC_E = potential_->quadrature_values()["FUNCTIONAL"];
    }

    double exchange_E = 0.0;
    double alpha = functional_->x_alpha();
    double beta = 1.0 - alpha;
    if (functional_->is_x_hybrid()) {
        exchange_E -= alpha * Da_->vector_dot(K_);
    }
    if (functional_->is_x_lrc()) {
        exchange_E -=  beta * Da_->vector_dot(wK_);
    }


    double two_electron_E = D_->vector_dot(Fa_) - 0.5 * one_electron_E;

    energies_["Nuclear"] = nuclearrep_;
    energies_["One-Electron"] = one_electron_E;
    energies_["Two-Electron"] =  coulomb_E + exchange_E;
    energies_["XC"] = XC_E;
    energies_["-D"] = variables_["-D Energy"];
    double dashD_E = energies_["-D"];

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += coulomb_E;
    Etotal += exchange_E;
    Etotal += XC_E;
    Etotal += dashD_E;

    return Etotal;
}

void RHF::Hx(SharedMatrix x, SharedMatrix IFock, SharedMatrix Cocc, SharedMatrix Cvir, SharedMatrix ret)
{
    if (functional_->needs_xc()){
        throw PSIEXCEPTION("SCF: Cannot yet compute DFT Hessian-vector prodcuts.\n");
    }
    // => Effective one electron part <= //
    Dimension virpi = Cvir->colspi();
    for (size_t h=0; h<nirrep_; h++){
        if (!doccpi_[h] || !virpi[h]) continue;
        double** IFp = IFock->pointer(h);
        double** retp = ret->pointer(h);
        double** xp = x->pointer(h);

        // ret_ia = F_ij X_ja
        C_DGEMM('N','N',doccpi_[h],virpi[h],doccpi_[h],1.0,
                IFp[0],nmopi_[h],
                xp[0],virpi[h],0.0,retp[0],virpi[h]);

        // ret_ia -= X_ib F_ba
        C_DGEMM('N','N',doccpi_[h],virpi[h],virpi[h],-1.0,
                xp[0],virpi[h],
                (IFp[doccpi_[h]]+doccpi_[h]),nmopi_[h],1.0,retp[0],virpi[h]);
    }

    // => Two electron part <= //
    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    Cl.push_back(Cocc);

    SharedMatrix R = Matrix::doublet(Cvir, x, false, true);
    R->scale(-1.0);
    Cr.push_back(R);

    jk_->compute();

    // Just in case someone only clears out Cleft and gets very strange errors
    Cl.clear();
    Cr.clear();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    // D_nm = Cocc_ni x_ia Cvir_ma
    // Cocc_ni (4 * J[D]_nm - K[D]_nm - K[D]_mn) C_vir_ma
    J[0]->scale(4.0);
    J[0]->subtract(K[0]);
    J[0]->subtract(K[0]->transpose());
    R->gemm(false, false, 1.0, J[0], Cocc, 0.0);
    ret->gemm(true, false, 1.0, R, Cvir, 1.0);
    ret->scale(-4.0);

    // Cleaup
    R.reset();
}

int RHF::soscf_update()
{
    time_t start, stop;
    start = time(NULL);

    // => Build gradient and preconditioner <= //

    // Grab occ and vir orbitals
    SharedMatrix Cocc = Ca_subset("SO", "OCC");
    SharedMatrix Cvir = Ca_subset("SO", "VIR");
    Dimension virpi = Cvir->colspi();

    // MO Fock Matrix (Inactive Fock in Helgaker's language)
    SharedMatrix IFock = Matrix::triplet(Ca_, Fa_, Ca_, true, false, false);

    SharedMatrix Gradient = SharedMatrix(new Matrix("Gradient", nirrep_, doccpi_, virpi));
    SharedMatrix Precon = SharedMatrix(new Matrix("Precon", nirrep_, doccpi_, virpi));

    for (size_t h=0; h<nirrep_; h++){

        if (!doccpi_[h] || !virpi[h]) continue;
        double* gp = Gradient->pointer(h)[0];
        double* denomp = Precon->pointer(h)[0];
        double** fp = IFock->pointer(h);

        for (size_t i=0, target=0; i<doccpi_[h]; i++){
            for (size_t a=doccpi_[h]; a < nmopi_[h]; a++){
                gp[target] = -4.0 * fp[i][a];
                denomp[target++] = -4.0 * (fp[i][i] - fp[a][a]);
            }
        }
    }

    // Make sure the MO gradient is reasonably small
    if (Gradient->absmax() > 0.3){
        if (print_ > 1){
            outfile->Printf("    Gradient element too large for SOSCF, using DIIS.\n");
        }
        return 0;
    }

    if (soscf_print_){
        outfile->Printf("\n");
        outfile->Printf("    ==> SORHF Iterations <==\n");
        outfile->Printf("    Maxiter     = %11d\n", soscf_max_iter_);
        outfile->Printf("    Miniter     = %11d\n", soscf_min_iter_);
        outfile->Printf("    Convergence = %11.3E\n", soscf_conv_);
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("    %-4s   %11s     %10s\n", "Iter", "Residual RMS", "Time [s]");
        outfile->Printf("    ---------------------------------------\n");
    }

    // => Initial CG guess <= //
    SharedMatrix x = Gradient->clone();
    x->apply_denominator(Precon);

    // Calc hessian vector product, find residual and conditioned residual
    SharedMatrix r = Gradient->clone();
    SharedMatrix Ap = SharedMatrix(new Matrix("Ap", nirrep_, doccpi_, virpi));
    Hx(x, IFock, Cocc, Cvir, Ap);
    r->subtract(Ap);

    // Print iteration 0 timings and rms
    double rconv = r->sum_of_squares();
    double grad_rms = Gradient->sum_of_squares();
    if (grad_rms < 1.e-14){
        grad_rms = 1.e-14; // Prevent rel denom from being too small
    }
    double rms = sqrt(rconv / grad_rms);
    stop = time(NULL);
    if (soscf_print_){
        outfile->Printf("    %-5s %11.3E %10ld\n", "Guess", rms, stop-start);
    }

    // Build new p and z vectors
    SharedMatrix z = r->clone();
    z->apply_denominator(Precon);
    SharedMatrix p = z->clone();

    // => CG iterations <= //
    int fock_builds = 1;
    for (int cg_iter=1; cg_iter<soscf_max_iter_; cg_iter++) {

        // Calc hessian vector product
        Hx(p, IFock, Cocc, Cvir, Ap);
        fock_builds += 1;

        // Find factors and scale
        double rzpre = r->vector_dot(z);
        double alpha = rzpre / p->vector_dot(Ap);
        if (std::isnan(alpha)){
            outfile->Printf("RHF::SOSCF Warning CG alpha is zero/nan. Stopping CG.\n");
            alpha = 0.0;
        }

        x->axpy(alpha, p);
        r->axpy(-alpha, Ap);

        // Get residual
        double rconv = r->sum_of_squares();
        double rms = sqrt(rconv / grad_rms);
        stop = time(NULL);
        if (soscf_print_){
            outfile->Printf("    %-5d %11.3E %10ld\n", cg_iter, rms, stop-start);
        }

        // Check convergence
        if (((rms < soscf_conv_) && (cg_iter >= soscf_min_iter_)) || (alpha==0.0)) {
            break;
        }

        // Update p and z
        z->copy(r);
        z->apply_denominator(Precon);

        double beta = r->vector_dot(z) / rzpre;

        p->scale(beta);
        p->add(z);

    }
    if (soscf_print_){
        outfile->Printf("    ---------------------------------------\n");
        outfile->Printf("\n");
    }

    // => Rotate orbitals <= //
    rotate_orbitals(Ca_, x);

    // => Cleanup <= //
    Cocc.reset();
    Cvir.reset();
    IFock.reset();
    Precon.reset();
    Gradient.reset();
    Ap.reset();
    z.reset();
    r.reset();
    p.reset();

    return fock_builds;
}

bool RHF::stability_analysis()
{
    if(scf_type_ == "DF" || scf_type_ == "CD"){
        throw PSIEXCEPTION("Stability analysis has not been implemented for density fitted wavefunctions yet.");
    }else{
#define ID(x) ints.DPD_ID(x)
        // Build the Fock Matrix
        SharedMatrix moF(new Matrix("MO basis fock matrix", nmopi_, nmopi_));
        moF->transform(Fa_, Ca_);

        std::vector<std::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::occ);
        spaces.push_back(MOSpace::vir);
        IntegralTransform ints(shared_from_this(), spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly,
                               IntegralTransform::QTOrder, IntegralTransform::None);
        ints.set_keep_dpd_so_ints(true);
        ints.transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
        ints.transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
        dpd_set_default(ints.get_dpd_id());
        dpdbuf4 Asing, Atrip,I;
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
        // Singlet A_ia_jb = 4 (ia|jb)
        global_dpd_->buf4_scmcopy(&I, PSIF_LIBTRANS_DPD, "RHF Singlet Hessian (IA|JB)", 4.0);
        // Triplet A_ia_jb = -(ib|ja)
        global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, psrq,
                           ID("[O,V]"), ID("[O,V]"), "RHF Triplet Hessian (IA|JB)", -1.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
        // Triplet A_ia_jb -= (ij|ab)
        global_dpd_->buf4_sort_axpy(&I, PSIF_LIBTRANS_DPD, prqs,
                           ID("[O,V]"), ID("[O,V]"), "RHF Triplet Hessian (IA|JB)", -1.0);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_init(&Atrip, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "RHF Triplet Hessian (IA|JB)");
        for(int h = 0; h < Atrip.params->nirreps; ++h){
            global_dpd_->buf4_mat_irrep_init(&Atrip, h);
            global_dpd_->buf4_mat_irrep_rd(&Atrip, h);
            for(int ia = 0; ia < Atrip.params->rowtot[h]; ++ia){
                int iabs = Atrip.params->roworb[h][ia][0];
                int aabs = Atrip.params->roworb[h][ia][1];
                int isym = Atrip.params->psym[iabs];
                int asym = Atrip.params->qsym[aabs];
                int irel = iabs - Atrip.params->poff[isym];
                int arel = aabs - Atrip.params->qoff[asym] + doccpi_[asym];
                for(int jb = 0; jb < Atrip.params->coltot[h]; ++jb){
                    int jabs = Atrip.params->colorb[h][jb][0];
                    int babs = Atrip.params->colorb[h][jb][1];
                    int jsym = Atrip.params->rsym[jabs];
                    int bsym = Atrip.params->ssym[babs];
                    int jrel = jabs - Atrip.params->roff[jsym];
                    int brel = babs - Atrip.params->soff[bsym] + doccpi_[bsym];
                    // Triplet A_ia_jb += delta_ij F_ab - delta_ab F_ij
                    if((iabs == jabs) && (asym == bsym))
                        Atrip.matrix[h][ia][jb] += moF->get(asym, arel, brel);
                    if((aabs == babs) && (isym == jsym))
                        Atrip.matrix[h][ia][jb] -= moF->get(isym, irel, jrel);
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&Atrip, h);
        }
        // Singlet A += Triplet A
        global_dpd_->buf4_init(&Asing, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "RHF Singlet Hessian (IA|JB)");
        global_dpd_->buf4_axpy(&Atrip, &Asing, 1.0);
        global_dpd_->buf4_close(&Atrip);
        global_dpd_->buf4_close(&Asing);

        /*
         *  Perform the stability analysis
         */
        std::vector<std::pair<double, int> >singlet_eval_sym;
        std::vector<std::pair<double, int> >triplet_eval_sym;

        global_dpd_->buf4_init(&Asing, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "RHF Singlet Hessian (IA|JB)");
        global_dpd_->buf4_init(&Atrip, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "RHF Triplet Hessian (IA|JB)");
        for(int h = 0; h < Asing.params->nirreps; ++h) {
            int dim = Asing.params->rowtot[h];
            if(dim == 0) continue;
            double *evals = init_array(dim);
            double **evecs = block_matrix(dim, dim);

            global_dpd_->buf4_mat_irrep_init(&Asing, h);
            global_dpd_->buf4_mat_irrep_rd(&Asing, h);
            sq_rsp(dim, dim, Asing.matrix[h], evals, 1, evecs, 1e-12);
            global_dpd_->buf4_mat_irrep_close(&Asing, h);

            int mindim = dim < 5 ? dim : 5;
            for(int i = 0; i < mindim; i++)
                singlet_eval_sym.push_back(std::make_pair(evals[i], h));

            zero_arr(evals, dim);
            zero_mat(evecs, dim, dim);

            global_dpd_->buf4_mat_irrep_init(&Atrip, h);
            global_dpd_->buf4_mat_irrep_rd(&Atrip, h);
            sq_rsp(dim, dim, Atrip.matrix[h], evals, 1, evecs, 1e-12);
            global_dpd_->buf4_mat_irrep_close(&Atrip, h);

            for(int i = 0; i < mindim; i++)
                triplet_eval_sym.push_back(std::make_pair(evals[i], h));

            free_block(evecs);
            delete [] evals;
        }

        outfile->Printf( "    Lowest singlet (RHF->RHF) stability eigenvalues:-\n");
        print_stability_analysis(singlet_eval_sym);
        outfile->Printf( "    Lowest triplet (RHF->UHF) stability eigenvalues:-\n");
        print_stability_analysis(triplet_eval_sym);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }

    // FOLLOW is not implemented for RHF
    return false;
}

}}
