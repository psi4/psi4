/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include <psi4-dec.h>

#include "dcft.h"
#include "defines.h"
#include <vector>
#include <liboptions/liboptions.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include <libdpd/dpd.h>
#include <libdiis/diismanager.h>

#include <lib3index/3index.h>
#include <libmints/mints.h>
#include <libmints/sieve.h>
#include <libfock/jk.h>
#include <libfock/apps.h>
#include <physconst.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;

namespace psi { namespace dcft{

/**
  * Build the density-fitting tensor: b(Q|mn) in SO-basis
  * b(Q|mn) = Sum_P (mn|P) [J^-1/2]_PQ
  * where J is the matrix of (P|Q)
  */
void DCFTSolver::df_build_b_so()
{
    dcft_timer_on("DCFTSolver::df_build_b_so()");

    outfile->Printf( "\n\n\t                  ************************************************\n");
    outfile->Printf(     "\t                  *        Density Fitting Module in DCFT        *\n");
    outfile->Printf(     "\t                  *                by Xiao Wang                  *\n");
    outfile->Printf(     "\t                  ************************************************\n");
    outfile->Printf( "\n");

    primary_ = BasisSet::pyconstruct_orbital(molecule_,
                                             "BASIS", options_.get_str("BASIS"));
    auxiliary_ = BasisSet::pyconstruct_auxiliary(molecule_,
                                                 "DF_BASIS_DCFT", options_.get_str("DF_BASIS_DCFT"),
                                                 "RIFIT", options_.get_str("BASIS"));
    boost::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());

    nn_ = primary_->nbf();
    nQ_ = auxiliary_->nbf();

    // Print memory
    // TODO: print memory information for UHF
    df_memory();

    // Form J(P,Q)^-1/2
    dcft_timer_on("DCFTSolver::Form J^-1/2");
    formJm12(auxiliary_, zero);
    dcft_timer_off("DCFTSolver::Form J^-1/2");

    // Form B(Q, mu, nu)
    dcft_timer_on("DCFTSolver::Form B(Q,mn)");
    formbso(primary_, auxiliary_, zero);
    dcft_timer_off("DCFTSolver::Form B(Q,mn)");

    dcft_timer_off("DCFTSolver::df_build_b_so()");
}

/**
  * Form J(P,Q)^-1/2
  */
void DCFTSolver::formJm12(boost::shared_ptr<BasisSet> auxiliary, boost::shared_ptr<BasisSet> zero)
{
//    outfile->Printf("\tForming J(P,Q)^-1/2 ...\n\n");
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    double **J = block_matrix(nQ_, nQ_);
    Jm12_ = block_matrix(nQ_, nQ_);

    // => Integrals <= //
    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary, zero, auxiliary, zero));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > Jint;
    std::vector<const double*> buffer;
    for (int t = 0; t < nthreads; t++){
        Jint.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri()));
        buffer.push_back(Jint[t]->buffer());
    }

    std::vector<std::pair<int, int> > PQ_pairs;
    for (int P = 0; P < auxiliary->nshell(); P++){
        for (int Q = 0; Q <= P; Q++){
            PQ_pairs.push_back(std::pair<int, int>(P, Q));
        }
    }

    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++){
        int P = PQ_pairs[PQ].first;
        int Q = PQ_pairs[PQ].second;

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        Jint[thread]->compute_shell(P, 0, Q, 0);

        int nP = auxiliary->shell(P).nfunction();
        int oP = auxiliary->shell(P).function_index();

        int nQ = auxiliary->shell(Q).nfunction();
        int oQ = auxiliary->shell(Q).function_index();

        int index = 0;
        for (int p = 0; p < nP; p++){
            for (int q = 0; q < nQ; q++, ++index){
                J[p + oP][q + oQ] = buffer[thread][index];
            }
        }
    }

    // First, diagonalize J(P,Q)
    int lwork = nQ_ * 3;
    double* eigval = init_array(nQ_);
    double* work = init_array(lwork);
    int status = C_DSYEV('v', 'u', nQ_, J[0], nQ_, eigval, work, lwork);
    if(status){
        throw PsiException("Diagonalization of J failed", __FILE__, __LINE__);
    }
    free(work);

    // Now J contains eigenvectors of the original J
    double **J_copy = block_matrix(nQ_, nQ_);
    C_DCOPY(nQ_ * nQ_, J[0], 1, J_copy[0], 1);

    // Now form J^-1/2 = U(T) * j^-1/2 * U
    // where j^-1/2 is the diagonal matrix of inverse square
    // of the eigenvalues, and U is the matrix of eigenvectors of J
    for (int i = 0; i < nQ_; ++i){
        eigval[i] = (eigval[i] < 1.0E-10) ? 0.0 : 1.0 / sqrt(eigval[i]);
        // scale one set of eigenvectors by the diagonal elements j^-1/2
        C_DSCAL(nQ_, eigval[i], J[i], 1);
    }
    free(eigval);

    // J^-1/2 = J_copy(T) * J
    C_DGEMM('t', 'n', nQ_, nQ_, nQ_, 1.0, J_copy[0], nQ_, J[0], nQ_, 0.0, Jm12_[0], nQ_);
    free_block(J);
    free_block(J_copy);

    // Save J^-1/2
    Jmhalf_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS Jmhalf <P|Q>", nQ_, nQ_));
    Jmhalf_->set(Jm12_);
}

/**
  * Form b(Q|mn)
  */
void DCFTSolver::formbso(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<BasisSet> auxiliary, boost::shared_ptr<BasisSet> zero)
{
//    outfile->Printf("\tForming b(Q|mu,nu) ...\n\n");
    bQso_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|mn)", nQ_, nso_, nso_));
    double **Ap = block_matrix(nQ_, nso_ * nso_);
    double **Bp = block_matrix(nQ_, nso_ * nso_);

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    boost::shared_ptr<ERISieve> sieve = boost::shared_ptr<ERISieve>(new ERISieve(primary, 1.0E-20));
    const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //
    int max_rows;
    max_rows = auxiliary->nshell();

    // => Block Sizing <= //
    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < auxiliary->nshell(); P++) {
        int nP = auxiliary->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(auxiliary->nshell());

    // => Integrals <= //
    boost::shared_ptr<IntegralFactory> rifactory2(new IntegralFactory(auxiliary, zero, primary, primary));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int t = 0; t < nthreads; t++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory2->eri()));
        buffer.push_back(eri[t]->buffer());
    }

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {

        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop  = Pstarts[block+1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary->shell(Pstart).function_index();
        int pstop  = (Pstop == auxiliary->nshell() ? nQ_ : auxiliary->shell(Pstop ).function_index());
        int np = pstop - pstart;

        // > Integrals < //
        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int P =  PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell(P,0,M,N);

            int nP = auxiliary->shell(P).nfunction();
            int oP = auxiliary->shell(P).function_index();

            int nM = primary->shell(M).nfunction();
            int oM = primary->shell(M).function_index();

            int nN = primary->shell(N).nfunction();
            int oN = primary->shell(N).function_index();

            int index = 0;
            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++, index++) {
                         Bp[p + oP][(m + oM) * nso_ + (n + oN)] = buffer[thread][index];
                         Bp[p + oP][(n + oN) * nso_ + (m + oM)] = buffer[thread][index];
                    }
                }
            }
        }
    }

    C_DGEMM('N','N', nQ_, nso_ * nso_, nQ_, 1.0, Jm12_[0], nQ_, Bp[0], nso_ * nso_, 0.0, Ap[0], nso_ * nso_);
    bQso_->set(Ap);
}

/**
  * Calculate memory required for density-fitting
  */
void DCFTSolver::df_memory(){
    double memory = Process::environment.get_memory();
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    outfile->Printf("\t => Sizing <=\n\n");
    outfile->Printf("\t  Memory   = %11d MB\n", long(memory) / (1024L * 1024L));
    outfile->Printf("\t  Threads  = %11d\n", nthreads);
    outfile->Printf("\t  nn       = %11d\n", nn_);
    outfile->Printf("\t  nQ       = %11d\n\n", nQ_);
    outfile->Printf("\t => Primary Basis <=\n\n");
    primary_->print();
    outfile->Printf("\t => Auxiliary Basis <=\n\n");
    auxiliary_->print();

    // Memory requirements
    outfile->Printf("\t => Memory Requirement <=\n\n");
    double cost_df = 0.0;
    cost_df += nQ_ * nQ_; // store J(P|Q)-1/2
    cost_df += nQ_ * nso_ * nso_; // store b(Q|mn)
    cost_df += nQ_ * nso_ * navir_; // store b(Q|mV)
    cost_df += nQ_ * navir_ * navir_; // store b(Q|VV)
    cost_df *= sizeof(double);
    cost_df /= 1024.0 * 1024.0;
    outfile->Printf("\tMemory required for 3-index integrals   : %9.2lf MB \n", cost_df);

    double cost_abcd = 0.0;
    cost_abcd += 2.0 * pow(navir_, 4.0); // store g(ab|cd) and g<ac|bd>
    cost_abcd += 2.0 * nalpha_ * nalpha_ * navir_ * navir_; // store lambda<oo|vv> and G<ij|ab> intermediates
    cost_abcd *= sizeof(double);
    cost_abcd /= 1024.0 * 1024.0;
    outfile->Printf("\tMemory required for ABCD-type integrals : %9.2lf MB \n\n", cost_abcd);

    double memory_mb = (double)memory / (1024.0 * 1024.0);
    outfile->Printf("\tMinimum Memory required                 : %9.2lf MB \n", cost_df + cost_abcd);
    outfile->Printf("\tMemory available                        : %9.2lf MB \n\n", memory_mb);
//    if(cost_df >= memory_mb)
//            throw PSIEXCEPTION("There is NOT enough memory for ABCD-type contraction!");
}

/**
  * Transform b(Q|mn) -> b(Q|pq)
  */
void DCFTSolver::transform_b()
{
//    outfile->Printf("\tTransform b(Q|mn) -> b(Q|pq) ...\n\n");
    dcft_timer_on("DCFTSolver::Transform B(Q,mn) -> B(Q,pq)");

    CvirA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Alpha C(mu,a)", nso_, navir_));
    #pragma omp parallel for
    for (int mu = 0; mu < nso_; mu++){
        for (int a = 0; a < navir_; a++){
            CvirA_->set(mu, a, Ca_->get(mu, a + nalpha_));
        }
    }
    CoccA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Alpha C(mu,i)", nso_, nalpha_));
    #pragma omp parallel for
    for (int mu = 0; mu < nso_; mu++){
        for (int i = 0; i < nalpha_; i++){
            CoccA_->set(mu, i, Ca_->get(mu, i));
        }
    }

    if (options_.get_str("REFERENCE") == "UHF"){
        CvirB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Beta C(mu, a)", nso_, nbvir_));
        #pragma omp parallel for
        for (int mu = 0; mu < nso_; mu++){
            for (int a = 0; a < nbvir_; a++){
                CvirB_->set(mu, a, Cb_->get(mu, a + nbeta_));
            }
        }
        CoccB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Beta C(mu, i)", nso_, nbeta_));
        #pragma omp parallel for
        for (int mu = 0; mu < nso_; mu++){
            for (int i = 0; i < nbeta_; i++){
                CoccB_->set(mu, i, Cb_->get(mu, i));
            }
        }
    }

    formb_oo();

    formb_ov();

    formb_vv();

    formb_pq();

    dcft_timer_off("DCFTSolver::Transform B(Q,mn) -> B(Q,pq)");
}

/**
  * form b(Q,ij)
  */
void DCFTSolver::formb_oo()
{
    dcft_timer_on("DCFTSolver::b(Q|mn) -> b(Q|ij)");
    bQnoA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|mO)", nQ_, nso_ * nalpha_));
    bQnoA_->contract(false, false, nQ_ * nso_, nalpha_, nso_, bQso_, CoccA_, 1.0, 0.0);

    bQooA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|OO)", nQ_, nalpha_, nalpha_));
    bQooA_->contract233(true, false, nalpha_, nalpha_, CoccA_, bQnoA_, 1.0, 0.0);
    bQnoA_.reset();

    if (options_.get_str("REFERENCE") == "UHF"){
        bQnoB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|mo)", nQ_, nso_ * nbeta_));
        bQnoB_->contract(false, false, nQ_ * nso_, nbeta_, nso_, bQso_, CoccB_, 1.0, 0.0);

        bQooB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|oo)", nQ_, nbeta_, nbeta_));
        bQooB_->contract233(true, false, nbeta_, nbeta_, CoccB_, bQnoB_, 1.0, 0.0);
        bQnoB_.reset();
    }

    dcft_timer_off("DCFTSolver::b(Q|mn) -> b(Q|ij)");

}

/**
  * form b(Q,ia)
  */
void DCFTSolver::formb_ov()
{
    dcft_timer_on("DCFTSolver::b(Q|mn) -> b(Q|ia)");
    bQnvA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|mV)", nQ_, nso_ * navir_));
    bQnvA_->contract(false, false, nQ_ * nso_, navir_, nso_, bQso_, CvirA_, 1.0, 0.0);

    bQovA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|OV)", nQ_, nalpha_, navir_));
    bQovA_->contract233(true, false, nalpha_, navir_, CoccA_, bQnvA_, 1.0, 0.0);

    if (options_.get_str("REFERENCE") == "UHF"){
        bQnvB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|mv)", nQ_, nso_ * nbvir_));
        bQnvB_->contract(false, false, nQ_ * nso_, nbvir_, nso_, bQso_, CvirB_, 1.0, 0.0);

        bQovB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|ov)", nQ_, nbeta_, nbvir_));
        bQovB_->contract233(true, false, nbeta_, nbvir_, CoccB_, bQnvB_, 1.0, 0.0);
    }

    dcft_timer_off("DCFTSolver::b(Q|mn) -> b(Q|ia)");

}

/**
  * form b(Q,ab)
  */
void DCFTSolver::formb_vv()
{
    dcft_timer_on("DCFTSolver::b(Q|mn) -> b(Q|ab)");
    bQvvA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|VV)", nQ_, navir_, navir_));
    bQvvA_->contract233(true, false, navir_, navir_, CvirA_, bQnvA_, 1.0, 0.0);
    bQnvA_.reset();

    if (options_.get_str("REFERENCE") == "UHF"){
        bQvvB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|vv)", nQ_, nbvir_, nbvir_));
        bQvvB_->contract233(true, false, nbvir_, nbvir_, CvirB_, bQnvB_, 1.0, 0.0);
        bQnvB_.reset();
    }
    dcft_timer_off("DCFTSolver::b(Q|mn) -> b(Q|ab)");

    // Check:
//    bQvvA_->print();
}

/**
  * form b(Q,pq)
  */
void DCFTSolver::formb_pq()
{
    dcft_timer_on("DCFTSolver::b(Q|mn) -> b(Q|pq)");
    CallA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Alpha C(mu,p)", nso_, nso_));
    CallA_->set(Ca_);

    bQnqA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS Alpha B(Q|mq)", nQ_, nso_ * nso_));
    bQnqA_->contract(false, false, nQ_ * nso_, nso_, nso_, bQso_, CallA_, 1.0, 0.0);
    bQpqA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS Alpha B(Q|pq)", nQ_, nso_, nso_));
    bQpqA_->contract233(true, false, nso_, nso_, CallA_, bQnqA_, 1.0, 0.0);
    bQnqA_.reset();

    if (options_.get_str("REFERENCE") == "UHF"){
        CallB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Beta C(mu,p)", nso_, nso_));
        CallB_->set(Cb_);

        bQnqB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS Beta B(Q|mq)", nQ_, nso_ * nso_));
        bQnqB_->contract(false, false, nQ_ * nso_, nso_, nso_, bQso_, CallB_, 1.0, 0.0);
        bQpqB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS Beta B(Q|pq)", nQ_, nso_, nso_));
        bQpqB_->contract233(true, false, nso_, nso_, CallB_, bQnqB_, 1.0, 0.0);
        bQnqB_.reset();
    }
    dcft_timer_off("DCFTSolver::b(Q|mn) -> b(Q|pq)");

}

/**
 * Form density-fitted MO-basis TEI g(OV|OV)
 */
void DCFTSolver::form_df_g_ovov()
{
    dcft_timer_on("DCFTSolver::DF Transform_OVOV");
    dfoccwave::SharedTensor2d tei_mo;
    dpdbuf4 I;

    // g(ia|jb) = Sum_Q b(ia|Q) b(Q|jb)
    // Alpha-Alpha
    tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(IA|JB)",
                                                               nalpha_, navir_,
                                                               nalpha_, navir_));
    tei_mo->gemm(true, false, bQovA_, bQovA_, 1.0, 0.0);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&I, h);
        global_dpd_->buf4_mat_irrep_rd(&I, h);

        #pragma omp parallel for
        for (long int ia = 0; ia < I.params->rowtot[h]; ++ia){
            for (long int jb = 0; jb < I.params->coltot[h]; ++jb){
                I.matrix[h][ia][jb] = tei_mo->get(ia, jb);
            }

        }
        global_dpd_->buf4_mat_irrep_wrt(&I, h);
        global_dpd_->buf4_mat_irrep_close(&I, h);
    }
    global_dpd_->buf4_close(&I);
    tei_mo.reset();

    if (options_.get_str("REFERENCE") == "UHF"){
        // Alpha-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(IA|jb)",
                                                                   nalpha_, navir_,
                                                                   nbeta_, nbvir_));
        tei_mo->gemm(true, false, bQovA_, bQovB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                               ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ia = 0; ia < I.params->rowtot[h]; ++ia){
                for (long int jb = 0; jb < I.params->coltot[h]; ++jb){
                    I.matrix[h][ia][jb] = tei_mo->get(ia, jb);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();

        // Beta-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(ia|jb)",
                                                                   nbeta_, nbvir_,
                                                                   nbeta_, nbvir_));
        tei_mo->gemm(true, false, bQovB_, bQovB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                               ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ia = 0; ia < I.params->rowtot[h]; ++ia){
                for (long int jb = 0; jb < I.params->coltot[h]; ++jb){
                    I.matrix[h][ia][jb] = tei_mo->get(ia, jb);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();

    }

    dcft_timer_off("DCFTSolver::DF Transform_OVOV");

}

/**
 * Form density-fitted MO-basis TEI g(OO|OO)
 */
void DCFTSolver::form_df_g_oooo()
{
    dcft_timer_on("DCFTSolver::DF Transform_OOOO");
    dfoccwave::SharedTensor2d tei_mo;
    dpdbuf4 I;
    // g(ij|kl) = Sum_Q b(ij|Q) b(Q|jb)
    // Alpha-Alpha
    tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(IJ|KL)",
                                                               nalpha_, nalpha_,
                                                               nalpha_, nalpha_));
    tei_mo->gemm(true, false, bQooA_, bQooA_, 1.0, 0.0);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>=O]+"), ID("[O>=O]+"),
                           ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&I, h);
        global_dpd_->buf4_mat_irrep_rd(&I, h);
        # pragma omp parallel for
        for (long int ij = 0; ij < I.params->rowtot[h]; ++ij){
            for (long int kl = 0; kl < I.params->coltot[h]; ++kl){
                int i = I.params->roworb[h][ij][0];
                int j = I.params->roworb[h][ij][1];
                int k = I.params->colorb[h][kl][0];
                int l = I.params->colorb[h][kl][1];
                I.matrix[h][ij][kl] = tei_mo->get(nalpha_ * i + j, nalpha_ * k + l);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&I, h);
        global_dpd_->buf4_mat_irrep_close(&I, h);
    }
    global_dpd_->buf4_close(&I);
    tei_mo.reset();

    if (options_.get_str("REFERENCE") == "UHF"){
        // Alpha-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(IJ|kl)",
                                                                   nalpha_, nalpha_,
                                                                   nbeta_, nbeta_));
        tei_mo->gemm(true, false, bQooA_, bQooB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>=O]+"), ID("[o>=o]+"),
                               ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);
            # pragma omp parallel for
            for (long int ij = 0; ij < I.params->rowtot[h]; ++ij){
                for (long int kl = 0; kl < I.params->coltot[h]; ++kl){
                    int i = I.params->roworb[h][ij][0];
                    int j = I.params->roworb[h][ij][1];
                    int k = I.params->colorb[h][kl][0];
                    int l = I.params->colorb[h][kl][1];
                    I.matrix[h][ij][kl] = tei_mo->get(nalpha_ * i + j, nbeta_ * k + l);
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();

        // Beta-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(ij|kl)",
                                                                   nbeta_, nbeta_,
                                                                   nbeta_, nbeta_));
        tei_mo->gemm(true, false, bQooB_, bQooB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>=o]+"), ID("[o>=o]+"),
                               ID("[o>=o]+"), ID("[o>=o]+"), 0, "MO Ints (oo|oo)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);
            # pragma omp parallel for
            for (long int ij = 0; ij < I.params->rowtot[h]; ++ij){
                for (long int kl = 0; kl < I.params->coltot[h]; ++kl){
                    int i = I.params->roworb[h][ij][0];
                    int j = I.params->roworb[h][ij][1];
                    int k = I.params->colorb[h][kl][0];
                    int l = I.params->colorb[h][kl][1];
                    I.matrix[h][ij][kl] = tei_mo->get(nbeta_ * i + j, nbeta_ * k + l);
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();
    }
    dcft_timer_off("DCFTSolver::DF Transform_OOOO");

}

/**
 * Form density-fitted MO-basis TEI g(VV|OO)
 */
void DCFTSolver::form_df_g_vvoo()
{
    dcft_timer_on("DCFTSolver::DF Transform_OOVV");
    dfoccwave::SharedTensor2d tei_mo;
    dpdbuf4 I;

    if (options_.get_str("REFERENCE") == "RHF"){
        // g(ab|ij) = Sum_Q b(ab|Q) b(Q|ij)
        // Alpha-Alpha
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(VV|OO)",
                                                                   navir_, navir_,
                                                                   nalpha_, nalpha_));
        tei_mo->gemm(true, false, bQvvA_, bQooA_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V>=V]+"), ID("[O>=O]+"),
                               ID("[V>=V]+"), ID("[O>=O]+"), 0, "MO Ints (VV|OO)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ab = 0; ab < I.params->rowtot[h]; ++ab){
                for (long int ij = 0; ij < I.params->coltot[h]; ++ij){
                    int a = I.params->roworb[h][ab][0];
                    int b = I.params->roworb[h][ab][1];
                    int i = I.params->colorb[h][ij][0];
                    int j = I.params->colorb[h][ij][1];
                    I.matrix[h][ab][ij] = tei_mo->get(navir_*a + b, nalpha_*i + j);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();
    }
    else{
        // g(ab|ij) = Sum_Q b(ab|Q) b(Q|ij)
        // Alpha-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(VV|oo)",
                                                                   navir_, navir_,
                                                                   nbeta_, nbeta_));
        tei_mo->gemm(true, false, bQvvA_, bQooB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V>=V]+"), ID("[o>=o]+"),
                               ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ab = 0; ab < I.params->rowtot[h]; ++ab){
                for (long int ij = 0; ij < I.params->coltot[h]; ++ij){
                    int a = I.params->roworb[h][ab][0];
                    int b = I.params->roworb[h][ab][1];
                    int i = I.params->colorb[h][ij][0];
                    int j = I.params->colorb[h][ij][1];
                    I.matrix[h][ab][ij] = tei_mo->get(navir_*a + b, nbeta_*i + j);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();

        // g(ij|ab) = Sum_Q b(ij|Q) b(Q|ab)
        // Alpha-Alpha
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(OO|VV)",
                                                                   nalpha_, nalpha_,
                                                                   navir_, navir_));
        tei_mo->gemm(true, false, bQooA_, bQvvA_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>=O]+"), ID("[V>=V]+"),
                               ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ij = 0; ij < I.params->rowtot[h]; ++ij){
                for (long int ab = 0; ab < I.params->coltot[h]; ++ab){
                    int i = I.params->roworb[h][ij][0];
                    int j = I.params->roworb[h][ij][1];
                    int a = I.params->colorb[h][ab][0];
                    int b = I.params->colorb[h][ab][1];
                    I.matrix[h][ij][ab] = tei_mo->get(nalpha_*i + j, navir_*a + b);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();

        // Alpha-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(OO|vv)",
                                                                   nalpha_, nalpha_,
                                                                   nbvir_, nbvir_));
        tei_mo->gemm(true, false, bQooA_, bQvvB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>=O]+"), ID("[v>=v]+"),
                               ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ij = 0; ij < I.params->rowtot[h]; ++ij){
                for (long int ab = 0; ab < I.params->coltot[h]; ++ab){
                    int i = I.params->roworb[h][ij][0];
                    int j = I.params->roworb[h][ij][1];
                    int a = I.params->colorb[h][ab][0];
                    int b = I.params->colorb[h][ab][1];
                    I.matrix[h][ij][ab] = tei_mo->get(nalpha_*i + j, nbvir_*a + b);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();

        // Beta-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(oo|vv)",
                                                                   nbeta_, nbeta_,
                                                                   nbvir_, nbvir_));
        tei_mo->gemm(true, false, bQooB_, bQvvB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>=o]+"), ID("[v>=v]+"),
                               ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ij = 0; ij < I.params->rowtot[h]; ++ij){
                for (long int ab = 0; ab < I.params->coltot[h]; ++ab){
                    int i = I.params->roworb[h][ij][0];
                    int j = I.params->roworb[h][ij][1];
                    int a = I.params->colorb[h][ab][0];
                    int b = I.params->colorb[h][ab][1];
                    I.matrix[h][ij][ab] = tei_mo->get(nbeta_*i + j, nbvir_*a + b);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();
    }

    dcft_timer_off("DCFTSolver::DF Transform_OOVV");

}

/**
 * Form density-fitted MO-basis TEI g(VO|OO)
 */
void DCFTSolver::form_df_g_vooo()
{
    dcft_timer_on("DCFTSolver::DF Transform_VOOO");
    dfoccwave::SharedTensor2d tei_mo;
    dpdbuf4 I;

    // g(ai|jk) = Sum_Q b(ai|Q) (Q|jk)
    // Alpha-Alpha
    tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(AI|JK)",
                                                               navir_, nalpha_,
                                                               nalpha_, nalpha_));
    bQvoA_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|VO)", nQ_, navir_, nalpha_));
    bQvoA_->sort3a(132, nQ_, nalpha_, navir_, bQovA_, 1.0, 0.0);
    tei_mo->gemm(true, false, bQvoA_, bQooA_, 1.0, 0.0);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O>=O]+"),
                           ID("[V,O]"), ID("[O>=O]+"), 0, "MO Ints (VO|OO)");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&I, h);
        global_dpd_->buf4_mat_irrep_rd(&I, h);

        #pragma omp parallel for
        for (long int ai = 0; ai < I.params->rowtot[h]; ++ai){
            for (long int jk = 0; jk < I.params->coltot[h]; ++jk){
                int j = I.params->colorb[h][jk][0];
                int k = I.params->colorb[h][jk][1];
                I.matrix[h][ai][jk] = tei_mo->get(ai, nalpha_*j + k);
            }

        }
        global_dpd_->buf4_mat_irrep_wrt(&I, h);
        global_dpd_->buf4_mat_irrep_close(&I, h);
    }
    global_dpd_->buf4_close(&I);
    tei_mo.reset();

    if (options_.get_str("REFERENCE") == "UHF"){
        // g(ai|jk) = Sum_Q b(ai|Q) (Q|jk)
        // Alpha-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(AI|jk)",
                                                                   navir_, nalpha_,
                                                                   nbeta_, nbeta_));
        tei_mo->gemm(true, false, bQvoA_, bQooB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o>=o]+"),
                               ID("[V,O]"), ID("[o>=o]+"), 0, "MO Ints (VO|oo)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ai = 0; ai < I.params->rowtot[h]; ++ai){
                for (long int jk = 0; jk < I.params->coltot[h]; ++jk){
                    int j = I.params->colorb[h][jk][0];
                    int k = I.params->colorb[h][jk][1];
                    I.matrix[h][ai][jk] = tei_mo->get(ai, nbeta_*j + k);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();
        // Beta-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(ai|jk)",
                                                                   nbvir_, nbeta_,
                                                                   nbeta_, nbeta_));
        bQvoB_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF_BASIS B(Q|vo)", nQ_, nbvir_, nbeta_));
        bQvoB_->sort3a(132, nQ_, nbeta_, nbvir_, bQovB_, 1.0, 0.0);
        tei_mo->gemm(true, false, bQvoB_, bQooB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o>=o]+"),
                               ID("[v,o]"), ID("[o>=o]+"), 0, "MO Ints (vo|oo)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ai = 0; ai < I.params->rowtot[h]; ++ai){
                for (long int jk = 0; jk < I.params->coltot[h]; ++jk){
                    int j = I.params->colorb[h][jk][0];
                    int k = I.params->colorb[h][jk][1];
                    I.matrix[h][ai][jk] = tei_mo->get(ai, nbeta_*j + k);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();

        // g(jk|ai) = Sum_Q b(jk|Q) (Q|ai)
        // Alpha-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(JK|ai)",
                                                                   nalpha_, nalpha_,
                                                                   nbvir_, nbeta_));
        tei_mo->gemm(true, false, bQooA_, bQvoB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>=O]+"), ID("[v,o]"),
                               ID("[O>=O]+"), ID("[v,o]"), 0, "MO Ints (OO|vo)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int jk = 0; jk < I.params->rowtot[h]; ++jk){
                for (long int ai = 0; ai < I.params->coltot[h]; ++ai){
                    int j = I.params->roworb[h][jk][0];
                    int k = I.params->roworb[h][jk][1];
                    I.matrix[h][jk][ai] = tei_mo->get(nalpha_*j + k, ai);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();
    }

    dcft_timer_off("DCFTSolver::DF Transform_VOOO");

}

/**
 * Form density-fitted MO-basis TEI g(OV|VV)
 */
void DCFTSolver::form_df_g_ovvv()
{
    dcft_timer_on("DCFTSolver::DF Transform_OVVV");
    dfoccwave::SharedTensor2d tei_mo;
    dpdbuf4 I;

    // g(ia|bc) = Sum_Q b(ia|Q) (Q|bc)
    // Alpha-Alpha
    tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(IA|BC)",
                                                               nalpha_, navir_,
                                                               navir_, navir_));
    tei_mo->gemm(true, false, bQovA_, bQvvA_, 1.0, 0.0);
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V>=V]+"),
                           ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&I, h);
        global_dpd_->buf4_mat_irrep_rd(&I, h);

        #pragma omp parallel for
        for (long int ia = 0; ia < I.params->rowtot[h]; ++ia){
            for (long int bc = 0; bc < I.params->coltot[h]; ++bc){
                int b = I.params->colorb[h][bc][0];
                int c = I.params->colorb[h][bc][1];
                I.matrix[h][ia][bc] = tei_mo->get(ia, navir_*b + c);
            }

        }
        global_dpd_->buf4_mat_irrep_wrt(&I, h);
        global_dpd_->buf4_mat_irrep_close(&I, h);
    }
    global_dpd_->buf4_close(&I);
    tei_mo.reset();

    if (options_.get_str("REFERENCE") == "UHF"){
        // g(ia|bc) = Sum_Q b(ia|Q) (Q|bc)
        // Alpha-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(IA|bc)",
                                                                   nalpha_, navir_,
                                                                   nbvir_, nbvir_));
        tei_mo->gemm(true, false, bQovA_, bQvvB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v>=v]+"),
                               ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ia = 0; ia < I.params->rowtot[h]; ++ia){
                for (long int bc = 0; bc < I.params->coltot[h]; ++bc){
                    int b = I.params->colorb[h][bc][0];
                    int c = I.params->colorb[h][bc][1];
                    I.matrix[h][ia][bc] = tei_mo->get(ia, nbvir_*b + c);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();
        // Beta-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(ia|bc)",
                                                                   nbeta_, nbvir_,
                                                                   nbvir_, nbvir_));
        tei_mo->gemm(true, false, bQovB_, bQvvB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v>=v]+"),
                               ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int ia = 0; ia < I.params->rowtot[h]; ++ia){
                for (long int bc = 0; bc < I.params->coltot[h]; ++bc){
                    int b = I.params->colorb[h][bc][0];
                    int c = I.params->colorb[h][bc][1];
                    I.matrix[h][ia][bc] = tei_mo->get(ia, nbvir_*b + c);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();

        // g(bc|ia) = Sum_Q b(bc|Q) (Q|ia)
        // Alpha-Beta
        tei_mo = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(BC|ia)",
                                                                   navir_, navir_,
                                                                   nbeta_, nbvir_));
        tei_mo->gemm(true, false, bQvvA_, bQovB_, 1.0, 0.0);
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V>=V]+"), ID("[o,v]"),
                               ID("[V>=V]+"), ID("[o,v]"), 0, "MO Ints (VV|ov)");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            global_dpd_->buf4_mat_irrep_rd(&I, h);

            #pragma omp parallel for
            for (long int bc = 0; bc < I.params->rowtot[h]; ++bc){
                for (long int ia = 0; ia < I.params->coltot[h]; ++ia){
                    int b = I.params->roworb[h][bc][0];
                    int c = I.params->roworb[h][bc][1];
                    I.matrix[h][bc][ia] = tei_mo->get(navir_*b + c, ia);
                }

            }
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
        global_dpd_->buf4_close(&I);
        tei_mo.reset();
    }
    dcft_timer_off("DCFTSolver::DF Transform_OVVV");

}

/**
 * Compute the density-fitted <Vv||Vv> tensors in G intermediates
 * and contract with lambda_ijcd
 */

void DCFTSolver::build_DF_tensors_RHF()
{
    dcft_timer_on("DCFTSolver::build_df_tensors_RHF()");

    // Build Tau matrix in MO basis (All)
    mo_tauA_ = SharedMatrix(new Matrix("MO basis Tau", nirrep_, nmopi_, nmopi_));
    #pragma omp parallel for
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int j = 0; j < naoccpi_[h]; ++j){
                mo_tauA_->set(h, i, j, aocc_tau_->get(h, i, j));
            }
        }
    }

    #pragma omp parallel for
    for(int h = 0; h < nirrep_; ++h){
        for(int a = naoccpi_[h]; a < nmopi_[h]; ++a){
            for(int b = naoccpi_[h]; b < nmopi_[h]; ++b){
                mo_tauA_->set(h, a, b, avir_tau_->get(h, a - naoccpi_[h], b - naoccpi_[h]));
            }
        }
    }

    /* Build Intermediate G_ijab = Sum_cd lambda<ij|cd> gbar<ab|cd> */

    dfoccwave::SharedTensor2d gabcd, gacbd, lambda;
    dpdbuf4 L, G;

    // g(ab|cd) = Sum_Q b(ab|Q) b(Q|cd)
    gabcd = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(ab|cd)",
                                                              navir_, navir_,
                                                              navir_, navir_));
    gabcd->gemm(true, false, bQvvA_, bQvvA_, 1.0, 0.0);
    // g(ab|cd) -> g<ac|bd>
    gacbd = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G<ac|bd>",
                                                              navir_, navir_,
                                                              navir_, navir_));
    gacbd->sort(1324, gabcd, 1.0, 0.0);
    // Read lambda <ij|ab>
    lambda = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Lambda <OO|VV>",
                                                               nalpha_ * nalpha_,
                                                               navir_ * navir_));
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&L, h);
        global_dpd_->buf4_mat_irrep_rd(&L, h);
        lambda->set(L.matrix[h]);
        global_dpd_->buf4_mat_irrep_close(&L, h);
    }
    global_dpd_->buf4_close(&L);
    // Intermediate G_ijab = Sum_cd lambda <ij|cd> g<ab|cd>
    Gijab_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Tau(temp) <ij|ab>",
                                                              nalpha_ * nalpha_,
                                                              navir_, navir_));
    Gijab_->gemm(false, true, lambda, gacbd, 1.0, 0.0);

    // Write G_ijab to a dpdbuf4 object
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) SF <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);

        #pragma omp parallel for
        for (long int ij = 0; ij < G.params->rowtot[h]; ++ij){
            for (long int ab = 0; ab < G.params->coltot[h]; ++ab){
                double value = Gijab_->get(ij, ab);
                G.matrix[h][ij][ab] = value;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // Write MO Ints <AC|BD> to disk, if analytic gradients are requested
    if(options_.get_str("DERTYPE") == "FIRST") {

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                               ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&G, h);
            global_dpd_->buf4_mat_irrep_rd(&G, h);
            #pragma omp parallel for
            for (long int ab = 0; ab < G.params->rowtot[h]; ++ab){
                for (long int cd = 0; cd < G.params->coltot[h]; ++cd){
                    double value = gacbd->get(ab, cd);
                    G.matrix[h][ab][cd] = value;
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&G, h);
            global_dpd_->buf4_mat_irrep_close(&G, h);
        }
        global_dpd_->buf4_close(&G);
        psio_->close(PSIF_LIBTRANS_DPD, 1);

    }
    gabcd.reset();
    gacbd.reset();
    lambda.reset();
    Gijab_.reset();

    /* Build [Gbar*Gamma]<Q|P> */
    build_gbarGamma_RHF();

    dcft_timer_off("DCFTSolver::build_df_tensors_RHF()");
}

/**
 * Form MO-based contraction [Gbar*Gamma]<Q|P>
 * [Gbar*Gamma]<q|p> = Sum_rs Gbar<qs|pr> Gamma<r|s>
 */
void DCFTSolver::build_gbarGamma_RHF()
{
    dcft_timer_on("DCFTSolver::Gbar<QS|PR> Gamma<R|S>");

    // Form MO-based gamma<R|S> = kappa<R|S> + tau<R|S>
    mo_gammaA_ = SharedMatrix(new Matrix("MO-based Gamma Alpha", nso_, nso_));
    mo_gammaA_->copy(kappa_mo_a_);
    mo_gammaA_->add(mo_tauA_);
    mo_gbarGamma_A_ = SharedMatrix(new Matrix("MO-based Gbar_Gamma_A", nso_, nso_));

    // Form MO-based gbar*gamma
    /*
     *  Gbar_Gamma <q|p> = (2 g<qs|pr> - g<qs|rp>) Gamma<r|s>
     *                   = (2 g(qp|sr) - g(qr|sp)) Gamma<r|s>
     */

    dfoccwave::SharedTensor2d gqpsr, gtilde_qpsr, gamma_3ind, gbarGamma_3ind;
    // Build g<qs|pr> = g(qp|sr) = Sum_Q b(qp|Q) b(Q|sr)
    gqpsr = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(qp|sr)",
                                                              nso_, nso_,
                                                              nso_, nso_));
    dcft_timer_on("DCFTSolver::DF Transform <PQ|RS>");
    gqpsr->gemm(true, false, bQpqA_, bQpqA_, 1.0, 0.0);
    dcft_timer_off("DCFTSolver::DF Transform <PQ|RS>");

    // Build g~(qp|sr) = 2 g(qp|sr) - g(qr|sp)
    dcft_timer_on("DCFTSolver::Build g~(qp|sr)");
    gtilde_qpsr = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G~(qp|sr)",
                                                                   nso_, nso_,
                                                                   nso_, nso_));
    gtilde_qpsr->copy(gqpsr);
    gtilde_qpsr->sort(1432, gqpsr, -1.0, 2.0);
    dcft_timer_off("DCFTSolver::Build g~(qp|sr)");

    // Gamma<r|s> -> Gamma(1|rs)
    gamma_3ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Gamma (1|rs)", 1, nso_, nso_));
    gamma_3ind->copy(convert_2ind_to_3ind(mo_gammaA_));
    // [gbar_gamma](qp|1) = g~(qp|sr) Gamma(sr|1)
    gbarGamma_3ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Gbar*Gamma (qp|1)", nso_ * nso_, 1));
    gbarGamma_3ind->gemm(false, true, gtilde_qpsr, gamma_3ind, 1.0, 0.0);

    // [gbar_gamma](qp|1) -> [gbar_gamma]<q|p>
    mo_gbarGamma_A_->copy(backconvert_3ind_to_2ind(gbarGamma_3ind));
    gqpsr.reset();
    gtilde_qpsr.reset();
    gamma_3ind.reset();
    gbarGamma_3ind.reset();

    dcft_timer_off("DCFTSolver::Gbar<QS|PR> Gamma<R|S>");
}

/**
 * Compute the density-fitted <VV|VV>, <Vv||Vv>, and <vv|vv> tensors in G intermediates
 * and contract with lambda_ijcd
 */

void DCFTSolver::build_DF_tensors_UHF()
{
    dcft_timer_on("DCFTSolver::build_df_tensors_UHF");

    // Build Tau matrix in MO basis (All)
    // Alpha-Alpha
    mo_tauA_ = SharedMatrix(new Matrix("MO basis Tau Alpha", nirrep_, nmopi_, nmopi_));
    #pragma omp parallel for
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int j = 0; j < naoccpi_[h]; ++j){
                mo_tauA_->set(h, i, j, aocc_tau_->get(h, i, j));
            }
        }
    }
    #pragma omp parallel for
    for(int h = 0; h < nirrep_; ++h){
        for(int a = naoccpi_[h]; a < nmopi_[h]; ++a){
            for(int b = naoccpi_[h]; b < nmopi_[h]; ++b){
                mo_tauA_->set(h, a, b, avir_tau_->get(h, a - naoccpi_[h], b - naoccpi_[h]));
            }
        }
    }

    // Beta-Beta
    mo_tauB_ = SharedMatrix(new Matrix("MO basis Tau Beta", nirrep_, nmopi_, nmopi_));
    #pragma omp parallel for
    for (int h = 0; h < nirrep_; ++h){
        for (int i = 0; i < nboccpi_[h]; ++i){
            for (int j = 0; j < nboccpi_[h]; ++j){
                mo_tauB_->set(h, i, j, bocc_tau_->get(h, i, j));
            }
        }
    }
    #pragma omp parallel for
    for (int h = 0; h < nirrep_; ++h){
        for (int a = nboccpi_[h]; a < nmopi_[h]; ++a){
            for (int b = nboccpi_[h]; b < nmopi_[h]; ++b){
                mo_tauB_->set(h, a, b, bvir_tau_->get(h, a - nboccpi_[h], b - nboccpi_[h]));
            }
        }
    }

    /* Build Intermediate G_ijab = Sum_cd lambda<ij|cd> gbar<ab|cd> */

    dfoccwave::SharedTensor2d gabcd, gacbd, gbar_acbd, lambda;
    dpdbuf4 L, G;

    /********** Alpha-Alpha **********/
    // g(AB|CD) = Sum_Q b(AB|Q) b(Q|CD)
    gabcd = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(AB|CD)",
                                                              navir_, navir_,
                                                              navir_, navir_));
    gabcd->gemm(true, false, bQvvA_, bQvvA_, 1.0, 0.0);
    // g(AB|CD) -> g<AC|BD>
    gacbd = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G<AC|BD>",
                                                              navir_, navir_,
                                                              navir_, navir_));
    gacbd->sort(1324, gabcd, 1.0, 0.0);
    // gbar <AC||BD> = g<AC|BD> - g<AC|DB>
    gbar_acbd = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF Gbar <AC||BD>",
                                                                  navir_, navir_,
                                                                  navir_, navir_));
    gbar_acbd->copy(gacbd);
    gbar_acbd->sort(1243, gacbd, -1.0, 1.0);
    // Read Lambda <IJ|AB>
    lambda = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Lambda <OO|VV>",
                                                               nalpha_, nalpha_,
                                                               navir_, navir_));
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&L, h);
        global_dpd_->buf4_mat_irrep_rd(&L, h);
        lambda->set(L.matrix[h]);
        global_dpd_->buf4_mat_irrep_close(&L, h);
    }
    global_dpd_->buf4_close(&L);
    // Intermediate G_IJAB = 1/2 Sum_CD lambda<IJ|CD> gbar<AB|CD>
    Gijab_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Tau(temp) <IJ|AB>",
                                                               nalpha_, nalpha_,
                                                               navir_, navir_));
    Gijab_->gemm(false, true, lambda, gbar_acbd, 0.5, 0.0);
    // Write G_IJAB to a dpdbuf4 object
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) <OO|VV>");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for (long int ij = 0; ij < G.params->rowtot[h]; ++ij){
            for (long int ab = 0; ab < G.params->coltot[h]; ++ab){
                double value = Gijab_->get(ij, ab);
                G.matrix[h][ij][ab] = value;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);
    // Write MO Ints <VV|VV> to a dpdbuf4 object, if analytic gradients are requested
    if (options_.get_str("DERTYPE") == "FIRST"){
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                               ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&G, h);
            global_dpd_->buf4_mat_irrep_rd(&G, h);
            #pragma omp parallel for
            for (long int ab = 0; ab < G.params->rowtot[h]; ++ab){
                for (long int cd = 0; cd < G.params->coltot[h]; ++cd){
                    double value = gacbd->get(ab, cd);
                    G.matrix[h][ab][cd] = value;
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&G, h);
            global_dpd_->buf4_mat_irrep_close(&G, h);
        }
        global_dpd_->buf4_close(&G);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }
    gabcd.reset();
    gacbd.reset();
    gbar_acbd.reset();
    lambda.reset();
    Gijab_.reset();
    /********** End of Alpha-Alpha **********/

    /********** Alpha-Beta **********/
    // g(AB|cd) = Sum_Q b(AB|Q) b(Q|cd)
    gabcd = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(AB|cd)",
                                                              navir_, navir_,
                                                              nbvir_, nbvir_));
    gabcd->gemm(true, false, bQvvA_, bQvvB_, 1.0, 0.0);
    // g(AB|cd) -> g<Ac|Bd>
    gacbd = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G<Ac|Bd>",
                                                              navir_, nbvir_,
                                                              navir_, nbvir_));
    gacbd->sort(1324, gabcd, 1.0, 0.0);
    // gbar <Ac||Bd> = g<Ac|Bd>, so gbar_acbd is not needed
    // Read Lambda <Ij|Ab>
    lambda = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Lambda <Oo|Vv>",
                                                               nalpha_, nbeta_,
                                                               navir_, nbvir_));
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&L, h);
        global_dpd_->buf4_mat_irrep_rd(&L, h);
        lambda->set(L.matrix[h]);
        global_dpd_->buf4_mat_irrep_close(&L, h);
    }
    global_dpd_->buf4_close(&L);
    // Intermediate G_IjAb = Sum_Cd lambda<Ij|Cd> gbar<Ab|Cd>
    Gijab_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Tau(temp) <Ij|Ab>",
                                                               nalpha_, nbeta_,
                                                               navir_, nbvir_));
    Gijab_->gemm(false, true, lambda, gacbd, 1.0, 0.0);
    // Write G_IjAb to a dpdbuf4 object
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "tau(temp) <Oo|Vv>");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for (long int ij = 0; ij < G.params->rowtot[h]; ++ij){
            for (long int ab = 0; ab < G.params->coltot[h]; ++ab){
                double value = Gijab_->get(ij, ab);
                G.matrix[h][ij][ab] = value;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);
    // Write MO Ints <Vv|Vv> to a dpdbuf4 object, if analytic gradients are requested
    if (options_.get_str("DERTYPE") == "FIRST"){
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                               ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&G, h);
            global_dpd_->buf4_mat_irrep_rd(&G, h);
            #pragma omp parallel for
            for (long int ab = 0; ab < G.params->rowtot[h]; ++ab){
                for (long int cd = 0; cd < G.params->coltot[h]; ++cd){
                    double value = gacbd->get(ab, cd);
                    G.matrix[h][ab][cd] = value;
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&G, h);
            global_dpd_->buf4_mat_irrep_close(&G, h);
        }
        global_dpd_->buf4_close(&G);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }
    gabcd.reset();
    gacbd.reset();
    gbar_acbd.reset();
    lambda.reset();
    Gijab_.reset();
    /********** End of Alpha-Beta **********/

    /********** Beta-Beta **********/
    // g(ab|cd) = Sum_Q b(ab|Q) b(Q|cd)
    gabcd = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(ab|cd)",
                                                              nbvir_, nbvir_,
                                                              nbvir_, nbvir_));
    gabcd->gemm(true, false, bQvvB_, bQvvB_, 1.0, 0.0);
    // g(ab|cd) -> g<ac|bd>
    gacbd = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G<ac|bd>",
                                                              nbvir_, nbvir_,
                                                              nbvir_, nbvir_));
    gacbd->sort(1324, gabcd, 1.0, 0.0);
    // gbar <ac||bd> = g<ac|bd> - g<ac|db>
    gbar_acbd = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF Gbar <ac||bd>",
                                                                  nbvir_, nbvir_,
                                                                  nbvir_, nbvir_));
    gbar_acbd->copy(gacbd);
    gbar_acbd->sort(1243, gacbd, -1.0, 1.0);
    // Read Lambda <ij|ab>
    lambda = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Lambda <oo|vv>",
                                                               nbeta_, nbeta_,
                                                               nbvir_, nbvir_));
    global_dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&L, h);
        global_dpd_->buf4_mat_irrep_rd(&L, h);
        lambda->set(L.matrix[h]);
        global_dpd_->buf4_mat_irrep_close(&L, h);
    }
    global_dpd_->buf4_close(&L);
    // Intermediate G_ijab = 1/2 Sum_cd lambda<ij|cd> gbar<ab|cd>
    Gijab_ = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Tau(temp) <ij|ab>",
                                                               nbeta_, nbeta_,
                                                               nbvir_, nbvir_));
    Gijab_->gemm(false, true, lambda, gbar_acbd, 0.5, 0.0);
    // Write G_ijab to a dpdbuf4 object
    global_dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "tau(temp) <oo|vv>");
    for (int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for (long int ij = 0; ij < G.params->rowtot[h]; ++ij){
            for (long int ab = 0; ab < G.params->coltot[h]; ++ab){
                double value = Gijab_->get(ij, ab);
                G.matrix[h][ij][ab] = value;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);
    // Write MO Ints <vv|vv> to a dpdbuf4 object, if analytic gradients are requested
    if (options_.get_str("DERTYPE") == "FIRST"){
        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                               ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
        for (int h = 0; h < nirrep_; ++h){
            global_dpd_->buf4_mat_irrep_init(&G, h);
            global_dpd_->buf4_mat_irrep_rd(&G, h);
            #pragma omp parallel for
            for (long int ab = 0; ab < G.params->rowtot[h]; ++ab){
                for (long int cd = 0; cd < G.params->coltot[h]; ++cd){
                    double value = gacbd->get(ab, cd);
                    G.matrix[h][ab][cd] = value;
                }
            }
            global_dpd_->buf4_mat_irrep_wrt(&G, h);
            global_dpd_->buf4_mat_irrep_close(&G, h);
        }
        global_dpd_->buf4_close(&G);
        psio_->close(PSIF_LIBTRANS_DPD, 1);
    }
    gabcd.reset();
    gacbd.reset();
    gbar_acbd.reset();
    lambda.reset();
    Gijab_.reset();
    /********** End of Beta-Beta **********/

    /* Build [gbar*gamma]<q|p> */
    build_gbarGamma_UHF();

    dcft_timer_off("DCFTSolver::build_df_tensors_UHF");

}

/**
 * Form MO-based contraction [Gbar*Gamma]<q|p>
 * [Gbar*Gamma]<q|p> = Sum_rs Gbar<qs|pr> Gamma<r|s>
 */
void DCFTSolver::build_gbarGamma_UHF()
{
    dcft_timer_on("DCFTSolver::Gbar<QS|PR> Gamma<R|S>");

    // Form gamma<R|S> = kappa<R|S> + tau<R|S>
    mo_gammaA_ = SharedMatrix(new Matrix("MO-based Gamma Alpha", nso_, nso_));
    mo_gammaA_->copy(kappa_mo_a_);
    mo_gammaA_->add(mo_tauA_);
    mo_gammaB_ = SharedMatrix(new Matrix("MO-based Gamma Beta", nso_, nso_));
    mo_gammaB_->copy(kappa_mo_b_);
    mo_gammaB_->add(mo_tauB_);
    mo_gbarGamma_A_ = SharedMatrix(new Matrix("MO-based Gbar_Gamma_A", nso_, nso_));
    mo_gbarGamma_B_ = SharedMatrix(new Matrix("MO-based Gbar_Gamma_B", nso_, nso_));

    // Form MO-based gbar*gamma

    /*
     *  Gbar_Gamma <Q|P> = Gbar<QS|PR> Gamma<R|S> + Gbar<Qs|Pr> Gamma<r|s>
     *                   = (g<QS|PR> - g<QS|RP>) Gamma<R|S> + g<Qs|Pr> Gamma<r|s>
     *                   = (g(QP|SR) - g(QR|SP)) Gamma<R|S> + g(QP|sr) Gamma<r|s>
     */

    dfoccwave::SharedTensor2d gqpsr, gtilde_qpsr, gamma_3ind, gbarGamma_3ind;

    /********** Gbar<QS|PR> Gamma<R|S> **********/
    // Build g(QP|SR) = g<QS|PR> = Sum_Aux b(QP|Aux) b(Aux|SR)
    gqpsr = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(QP|SR)",
                                                              nso_, nso_,
                                                              nso_, nso_));
    gqpsr->gemm(true, false, bQpqA_, bQpqA_, 1.0, 0.0);
    // Build g~(QP|SR) = g(QP|SR) - g(QR|SP) = gbar<QS|PR>
    gtilde_qpsr = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G~(QP|SR)",
                                                                    nso_, nso_,
                                                                    nso_, nso_));
    gtilde_qpsr->copy(gqpsr);
    gtilde_qpsr->sort(1432, gqpsr, -1.0, 1.0);
    // Gamma<R|S> -> Gamma(1|RS)
    gamma_3ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Gamma (1|RS)", 1, nso_, nso_));
    gamma_3ind->copy(convert_2ind_to_3ind(mo_gammaA_));
    // [gbar_gamma](QP|1) = g~(QP|SR) Gamma(SR|1)
    gbarGamma_3ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Gbar*Gamma (QP|1)", nso_ * nso_, 1));
    gbarGamma_3ind->gemm(false, true, gtilde_qpsr, gamma_3ind, 1.0, 0.0);
    // [gbar_gamma](QP|1) -> [gbar_gamma]<Q|P>
    mo_gbarGamma_A_->copy(backconvert_3ind_to_2ind(gbarGamma_3ind));
    gqpsr.reset();
    gtilde_qpsr.reset();
    gamma_3ind.reset();
    gbarGamma_3ind.reset();

    /********** Gbar<Qs|Pr> Gamma<r|s> **********/
    // Build g(QP|sr) = g<Qs|Pr> = Sum_Aux b(QP|Aux) b(Aux|sr)
    gqpsr = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(QP|sr)",
                                                              nso_, nso_,
                                                              nso_, nso_));
    gqpsr->gemm(true, false, bQpqA_, bQpqB_, 1.0, 0.0);
    // Build g~(QP|sr) = g(QP|sr) = gbar<Qs|Pr>. So no need to build it.
    // Gamma<r|s> -> Gamma(1|rs)
    gamma_3ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Gamma (1|rs)", 1, nso_, nso_));
    gamma_3ind->copy(convert_2ind_to_3ind(mo_gammaB_));
    // [gbar_gamma](QP|1) = g~(QP|sr) Gamma(sr|1)
    gbarGamma_3ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Gbar Gamma (QP|1)", nso_ * nso_, 1));
    gbarGamma_3ind->gemm(false, true, gqpsr, gamma_3ind, 1.0, 0.0);
    // [gbar_gamma](QP|1) -> [gbar_gamma]<Q|P>
    mo_gbarGamma_A_->add(backconvert_3ind_to_2ind(gbarGamma_3ind));
    gqpsr.reset();
    gamma_3ind.reset();
    gbarGamma_3ind.reset();

    /*
     *  Gbar_Gamma <q|p> = Gbar<qs|pr> Gamma<r|s> + Gbar<qS|pR> Gamma<R|S>
     *                   = (g<qs|pr> - g<qs|rp>) Gamma<r|s> + g<qS|pR> Gamma<R|S>
     *                   = (g(qp|sr) - g(qr|sp)) Gamma<r|s> + g(qp|SR) Gamma<R|S>
     */

    /********** Gbar<qs|pr> Gamma<r|s> **********/
    // Build g(qp|sr) = g<qs|pr> = Sum_Aux b(qp|Aux) b(Aux|sr)
    gqpsr = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(qp|sr)",
                                                              nso_, nso_,
                                                              nso_, nso_));
    gqpsr->gemm(true, false, bQpqB_, bQpqB_, 1.0, 0.0);
    // Build g~(qp|sr) = g(qp|sr) - g(qr|sp) = gbar<qs|pr>
    gtilde_qpsr = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G~(qp|sr)",
                                                                    nso_, nso_,
                                                                    nso_, nso_));
    gtilde_qpsr->copy(gqpsr);
    gtilde_qpsr->sort(1432, gqpsr, -1.0, 1.0);
    // Gamma<r|s> -> Gamma(1|rs)
    gamma_3ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Gamma (1|rs)", 1, nso_, nso_));
    gamma_3ind->copy(convert_2ind_to_3ind(mo_gammaB_));
    // [gbar_gamma](qp|1) = g~(qp|sr) Gamma(sr|1)
    gbarGamma_3ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Gbar Gamma (qp|1)", nso_ * nso_, 1));
    gbarGamma_3ind->gemm(false, true, gtilde_qpsr, gamma_3ind, 1.0, 0.0);
    // [gbar_gamma](qp|1) -> [gbar_gamma]<q|p>
    mo_gbarGamma_B_->copy(backconvert_3ind_to_2ind(gbarGamma_3ind));
    gqpsr.reset();
    gtilde_qpsr.reset();
    gamma_3ind.reset();
    gbarGamma_3ind.reset();

    /********** Gbar<qS|pR> Gamma<R|S> **********/
    // Build g(qp|SR) = g<qS|pR> = Sum_Aux b(qp|Aux) b(Aux|SR)
    gqpsr = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("DF G(qp|SR)",
                                                              nso_, nso_,
                                                              nso_, nso_));
    gqpsr->gemm(true, false, bQpqB_, bQpqA_, 1.0, 0.0);
    // Build g~(qp|SR) = g(qp|SR) = gbar<qS|pR>. So no need to build it.
    // Gamma<R|S> -> Gamma(1|RS)
    gamma_3ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Gamma (1|rs)", 1, nso_, nso_));
    gamma_3ind->copy(convert_2ind_to_3ind(mo_gammaA_));
    // [gbar_gamma](qp|1) = g~(qp|SR) Gamma(SR|1)
    gbarGamma_3ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Gbar Gamma (qp|1)", nso_ * nso_, 1));
    gbarGamma_3ind->gemm(false, true, gqpsr, gamma_3ind, 1.0, 0.0);
    // [gbar_gamma](qp|1) -> [gbar_gamma]<q|p>
    mo_gbarGamma_B_->add(backconvert_3ind_to_2ind(gbarGamma_3ind));
    gqpsr.reset();
    gamma_3ind.reset();
    gbarGamma_3ind.reset();

    dcft_timer_off("DCFTSolver::Gbar<QS|PR> Gamma<R|S>");
}

/**
 * Convert 2-index <P|Q> to 3-index (1|PQ)
 *
 * @param X: X<P|Q>
 */
dfoccwave::SharedTensor2d DCFTSolver::convert_2ind_to_3ind(SharedMatrix X)
{
    // Debug use
    if (X->rowdim() == 1)
        throw PSIEXCEPTION("Dimension of X in convert_2ind_to_3ind(X) is wrong!");

    dfoccwave::SharedTensor2d three_ind = dfoccwave::SharedTensor2d(new dfoccwave::Tensor2d("Three-index tensor", 1, nso_, nso_));
    #pragma omp parallel for
    for (int p = 0; p < nso_; ++p){
        for (int q = 0; q < nso_; ++q){
            long int pq = p * nso_ + q;
            three_ind->set(0, pq, X->get(p, q));
        }
    }

    return three_ind;
}

/**
 * Back-convert 3-index (PQ|1) to 2-index <P|Q>
 *
 * @param X: X(PQ|1)
 */
SharedMatrix DCFTSolver::backconvert_3ind_to_2ind(dfoccwave::SharedTensor2d X)
{
    // Debug use
    if (X->dim2() != 1)
        throw PSIEXCEPTION("Dimension of X in backconvert_3ind_to_2ind(X) is wrong!");

    SharedMatrix two_ind = SharedMatrix(new Matrix("Two-index tensor", nso_, nso_));
    #pragma omp parallel for
    for (int p = 0; p < nso_; ++p){
        for (int q =  0; q < nso_; ++q){
            long int pq = p * nso_ + q;
            two_ind->set(p, q, X->get(pq, 0));
        }
    }

    return two_ind;
}

}} // End namespace
