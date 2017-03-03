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

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include "psi4/psifiles.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/view.h"
#include "psi4/libmints/integral.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/basisset.h"
#include "dcft.h"
#include "defines.h"
#include <vector>
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libdiis/diismanager.h"

#include "psi4/lib3index/3index.h"

#include "psi4/libmints/sieve.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/apps.h"
#include "psi4/physconst.h"

#ifdef _OPENMP
#include <omp.h>
#endif



namespace psi { namespace dcft{

/**
  * Build the density-fitting tensor: b(Q|mn) in AO-basis
  * b(Q|mn) = Sum_P (mn|P) [J^-1/2]_PQ
  * where J is the matrix of (P|Q)
  */
void DCFTSolver::df_build_b_ao()
{
    dcft_timer_on("DCFTSolver::df_build_b_ao()");

    outfile->Printf( "\n\n\t                  ************************************************\n");
    outfile->Printf(     "\t                  *        Density Fitting Module in DCFT        *\n");
    outfile->Printf(     "\t                  *                by Xiao Wang                  *\n");
    outfile->Printf(     "\t                  ************************************************\n");
    outfile->Printf( "\n");

    primary_ = get_basisset("ORBITAL");
    auxiliary_ = get_basisset("DF_BASIS_DCFT");
    auxiliary_scf_ = get_basisset("DF_BASIS_SCF");

    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());

    nn_ = primary_->nbf();
    nQ_ = auxiliary_->nbf();
    nQ_scf_ = auxiliary_scf_->nbf();

    // Print memory
    // TODO: print memory information for UHF
    df_memory();

    // Form J(P,Q)^-1/2
    dcft_timer_on("DCFTSolver::Form J^-1/2");
    formJm12(auxiliary_, zero);
    dcft_timer_off("DCFTSolver::Form J^-1/2");

    // Form B(Q, mu, nu)
    dcft_timer_on("DCFTSolver::Form B(Q,mn)");
    formb_ao(primary_, auxiliary_, zero);
    dcft_timer_off("DCFTSolver::Form B(Q,mn)");

    // Form J(P,Q)^-1/2 for SCF terms
    dcft_timer_on("DCFTSolver::Form J^-1/2 (SCF terms)");
    formJm12_scf(auxiliary_scf_, zero);
    dcft_timer_off("DCFTSolver::Form J^-1/2 (SCF terms)");

    // Form B(Q, mu, nu) for SCF terms
    dcft_timer_on("DCFTSolver::Form B(Q,mn) (SCF terms)");
    formb_ao_scf(primary_, auxiliary_scf_, zero);
    dcft_timer_off("DCFTSolver::Form B(Q,mn) (SCF terms)");

    dcft_timer_off("DCFTSolver::df_build_b_ao()");
}

/**
  * Form J(P,Q)^-1/2
  */
void DCFTSolver::formJm12(std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<BasisSet> zero)
{
//    outfile->Printf("\tForming J(P,Q)^-1/2 ...\n\n");
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    double **J = block_matrix(nQ_, nQ_);
    Jm12_ = block_matrix(nQ_, nQ_);

    // => Integrals <= //
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary, zero, auxiliary, zero));
    std::vector<std::shared_ptr<TwoBodyAOInt> > Jint;
    std::vector<const double*> buffer;
    for (int t = 0; t < nthreads; t++){
        Jint.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri()));
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

}

/**
  * Form b(Q|mn)
  */
void DCFTSolver::formb_ao(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<BasisSet> zero)
{
    bQmn_ao_ = SharedMatrix(new Matrix(nQ_, nso_ * nso_));
    double **Ap = bQmn_ao_->pointer();
    double **Bp = block_matrix(nQ_, nso_ * nso_);

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    std::shared_ptr<ERISieve> sieve = std::shared_ptr<ERISieve>(new ERISieve(primary, 1.0E-20));
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
    std::shared_ptr<IntegralFactory> rifactory2(new IntegralFactory(auxiliary, zero, primary, primary));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int t = 0; t < nthreads; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory2->eri()));
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
}

/**
  * Calculate memory required for density-fitting
  */
void DCFTSolver::df_memory(){
    double memory = Process::environment.get_memory();
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
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

    if (options_.get_str("REFERENCE") == "RHF"){
        cost_df += nQ_ * nQ_; // J(P|Q)-1/2
        cost_df += 2 * nQ_ * nso_ * nso_; // b(Q|mn)
        cost_df += nQ_ * nalpha_ * nalpha_; // b(Q|oo)
        cost_df += 2 * nQ_ * nalpha_ * navir_; // b(Q|ov) and b(Q|vo)
        cost_df += nQ_ * navir_ * navir_; // b(Q|vv)
        cost_df += nQ_ * nso_ * nso_; // b(Q|pq)
        cost_df += 2 * navirpi_.max() * navirpi_.max() * navirpi_.max(); // (V'V|VV)
    }
    else {
        cost_df += nQ_ * nQ_; // J(P|Q)-1/2
        cost_df += 2 * nQ_ * nso_ * nso_; // b(Q|mn)
        cost_df += 2 * nQ_ * nalpha_ * nalpha_; // b(Q|oo)
        cost_df += 4 * nQ_ * nalpha_ * navir_; // b(Q|ov) and b(Q|vo)
        cost_df += 2 * nQ_ * navir_ * navir_; // b(Q|vv)
        cost_df += 2 * nQ_ * nso_ * nso_; // b(Q|pq)
        cost_df += 2 * navirpi_.max() * navirpi_.max() * navirpi_.max(); // (V'V|VV)
    }

    cost_df *= sizeof(double);
    cost_df /= 1024.0 * 1024.0;

    double memory_mb = (double)memory / (1024.0 * 1024.0);
    outfile->Printf("\tMinimum Memory required                 : %9.2lf MB \n", cost_df);
    outfile->Printf("\tMemory available                        : %9.2lf MB \n\n", memory_mb);
//    if(cost_df >= memory_mb)
//            throw PSIEXCEPTION("There is NOT enough memory for ABCD-type contraction!");
}

/**
  * Transform b(Q|mn) -> b(Q|pq)
  */
void DCFTSolver::transform_b()
{
    dcft_timer_on("DCFTSolver::Transform B(Q,mn) -> B(Q,pq)");

    formb_oo();
    formb_ov();
    formb_vv();
    formb_pq();

    dcft_timer_off("DCFTSolver::Transform B(Q,mn) -> B(Q,pq)");
}

/**
  * Transform b(Q|mu,nu) from AO basis to SO basis
  */
void DCFTSolver::transform_b_ao2so()
{
    dcft_timer_on("DCFTSolver::Transform b(Q|mn) AO-basis -> SO-basis");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    double** bQmn_ao_p = bQmn_ao_->pointer();

    // Set up dimensions for SO-basis b(Q|mn)
    Dimension Q(nirrep_), mn(nirrep_);
    for (int hn = 0; hn < nirrep_; ++hn){
        Q[hn] = nQ_;
        for (int hm = 0; hm < nirrep_; ++hm){
            mn[hm ^ hn] += nsopi_[hm] * nsopi_[hn];
        }
    }
    bQmn_so_ = SharedMatrix(new Matrix("Fully-transformed b", Q, mn));

    std::vector<int> offset(nirrep_);
    for (int h = 0; h < nirrep_; ++h){
        offset.push_back(0);
    }

    // AO-basis b(Q|mn) -> SO-basis b(Q|mn)
    for (int h = 0; h < nirrep_; ++h){
        double** bQmn_so_p = bQmn_so_->pointer(h);
        for (int hm = 0; hm < nirrep_; ++hm){
            int hn = h ^ hm;
            if (nsopi_[hm] > 0 && nsopi_[hn] > 0){

                SharedMatrix tmp (new Matrix("Half-transformed b", nQ_, nso_ * nsopi_[hn]));
                double** tmpp = tmp->pointer();
                double** ao2so_n_p = reference_wavefunction()->aotoso()->pointer(hn);
                double** ao2so_m_p = reference_wavefunction()->aotoso()->pointer(hm);
                // First-half transformation
                C_DGEMM('N', 'N', nQ_ * nso_, nsopi_[hn], nso_, 1.0, bQmn_ao_p[0], nso_, ao2so_n_p[0], nsopi_[hn], 0.0, tmpp[0], nsopi_[hn]);
                // Second-half transformation
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int Q = 0; Q < nQ_; ++Q){
                    C_DGEMM('T', 'N', nsopi_[hm], nsopi_[hn], nso_, 1.0, ao2so_m_p[0], nsopi_[hm], tmpp[Q], nsopi_[hn], 0.0, bQmn_so_p[Q]+offset[h], nsopi_[hn]);
                }
            }
        offset[h] += nsopi_[hm] * nsopi_[hn];
        }
    }

    bQmn_ao_.reset();

    dcft_timer_off("DCFTSolver::Transform b(Q|mn) AO-basis -> SO-basis");
}

/**
  * form b(Q,ij)
  */
void DCFTSolver::formb_oo()
{
    dcft_timer_on("DCFTSolver::b(Q|mn) -> b(Q|ij)");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Set up dimensions for b(Q|IJ)
    Dimension OO(nirrep_), Q(nirrep_);
    for (int hI = 0; hI < nirrep_; ++hI){
        Q[hI] = nQ_;
        for (int hJ = 0; hJ < nirrep_; ++hJ){
            OO[hI ^ hJ] += naoccpi_[hI] * naoccpi_[hJ];
        }
    }
    bQijA_mo_ = SharedMatrix(new Matrix("b(Q|IJ)", Q, OO));

    std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
    for (int h = 0; h < nirrep_; ++h){
        offset_so.push_back(0);
        offset_mo.push_back(0);
    }

    for (int h = 0; h < nirrep_; ++h){
        double** bQmn_so_p = bQmn_so_->pointer(h);
        double** bQijA_mo_p = bQijA_mo_->pointer(h);
        for (int hI = 0; hI < nirrep_; ++hI){
            int hJ = h ^ hI;
            if (naoccpi_[hI] > 0 && naoccpi_[hJ] > 0){
                double** CaJp = Ca_->pointer(hJ);
                double** CaIp = Ca_->pointer(hI);
                SharedMatrix tmp (new Matrix("Half-transformed b_IJ", nQ_, nsopi_[hI] * naoccpi_[hJ]));
                double** tmpp = tmp->pointer();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int Q = 0; Q < nQ_; ++Q){
                    // First-half transformation
                    C_DGEMM('N', 'N', nsopi_[hI], naoccpi_[hJ], nsopi_[hJ], 1.0, bQmn_so_p[Q] + offset_so[h], nsopi_[hJ], CaJp[0], nsopi_[hJ], 0.0, tmpp[Q], naoccpi_[hJ]);
                    // Second-half transformation
                    C_DGEMM('T', 'N', naoccpi_[hI], naoccpi_[hJ], nsopi_[hI], 1.0, CaIp[0], nsopi_[hI], tmpp[Q], naoccpi_[hJ], 0.0, bQijA_mo_p[Q] + offset_mo[h], naoccpi_[hJ]);
                }
            }
            offset_so[h] += nsopi_[h ^ hI] * nsopi_[hI];
            offset_mo[h] += naoccpi_[h ^ hI] * naoccpi_[hI];
        }
    }

    if (options_.get_str("REFERENCE") != "RHF"){

        // Set up dimensions for b(Q|ij)
        Dimension oo(nirrep_), Q(nirrep_);
        for (int hi = 0; hi < nirrep_; ++hi){
            Q[hi] = nQ_;
            for (int hj = 0; hj < nirrep_; ++hj){
                oo[hi ^ hj] += nboccpi_[hi] * nboccpi_[hj];
            }
        }

        bQijB_mo_ = SharedMatrix(new Matrix("b(Q|ij)", Q, oo));

        std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
        for (int h = 0; h < nirrep_; ++h){
            offset_so.push_back(0);
            offset_mo.push_back(0);
        }

        for (int h = 0; h < nirrep_; ++h){
            double** bQmn_so_p = bQmn_so_->pointer(h);
            double** bQijB_mo_p = bQijB_mo_->pointer(h);
            for (int hi = 0; hi < nirrep_; ++hi){
                int hj = h ^ hi;
                if (nboccpi_[hi] > 0 && nboccpi_[hj] > 0){
                    double** Cbjp = Cb_->pointer(hj);
                    double** Cbip = Cb_->pointer(hi);
                    SharedMatrix tmp (new Matrix("Half-transformed b_ij", nQ_, nsopi_[hi] * nboccpi_[hj]));
                    double** tmpp = tmp->pointer();
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int Q = 0; Q < nQ_; ++Q){
                        // First-half transformation
                        C_DGEMM('N', 'N', nsopi_[hi], nboccpi_[hj], nsopi_[hj], 1.0, bQmn_so_p[Q] + offset_so[h], nsopi_[hj], Cbjp[0], nsopi_[hj], 0.0, tmpp[Q], nboccpi_[hj]);
                        // Second-half transformation
                        C_DGEMM('T', 'N', nboccpi_[hi], nboccpi_[hj], nsopi_[hi], 1.0, Cbip[0], nsopi_[hi], tmpp[Q], nboccpi_[hj], 0.0, bQijB_mo_p[Q] + offset_mo[h], nboccpi_[hj]);
                    }
                }
                offset_so[h] += nsopi_[h ^ hi] * nsopi_[hi];
                offset_mo[h] += nboccpi_[h ^ hi] * nboccpi_[hi];
            }
        }


    }

    dcft_timer_off("DCFTSolver::b(Q|mn) -> b(Q|ij)");
}

/**
  * form b(Q,ia)
  */
void DCFTSolver::formb_ov()
{
    dcft_timer_on("DCFTSolver::b(Q|mn) -> b(Q|ia)");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Set up dimensions for b(Q|IA)
    Dimension OV(nirrep_), Q(nirrep_);
    for (int hI = 0; hI < nirrep_; ++hI){
        Q[hI] = nQ_;
        for (int hA = 0; hA < nirrep_; ++hA){
            OV[hI ^ hA] += naoccpi_[hI] * navirpi_[hA];
        }
    }
    bQiaA_mo_ = SharedMatrix(new Matrix("b(Q|IA)", Q, OV));

    std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
    for (int h = 0; h < nirrep_; ++h){
        offset_so.push_back(0);
        offset_mo.push_back(0);
    }

    for (int h = 0; h < nirrep_; ++h){
        double** bQmn_so_p = bQmn_so_->pointer(h);
        double** bQiaA_mo_p = bQiaA_mo_->pointer(h);
        for (int hI = 0; hI < nirrep_; ++hI){
            int hA = h ^ hI;
            if (naoccpi_[hI] > 0 && navirpi_[hA] > 0){
                double** CaVp = Ca_->pointer(hA);
                double** CaOp = Ca_->pointer(hI);
                SharedMatrix tmp (new Matrix("Half-transformed b_OV", nQ_, nsopi_[hI] * navirpi_[hA]));
                double** tmpp = tmp->pointer();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int Q = 0; Q < nQ_; ++Q){
                    // First-half transformation
                    C_DGEMM('N', 'N', nsopi_[hI], navirpi_[hA], nsopi_[hA], 1.0, bQmn_so_p[Q] + offset_so[h], nsopi_[hA], CaVp[0] + naoccpi_[hA], nsopi_[hA], 0.0, tmpp[Q], navirpi_[hA]);
                    // Second-half transformation
                    C_DGEMM('T', 'N', naoccpi_[hI], navirpi_[hA], nsopi_[hI], 1.0, CaOp[0], nsopi_[hI], tmpp[Q], navirpi_[hA], 0.0, bQiaA_mo_p[Q] + offset_mo[h], navirpi_[hA]);
                }
            }
            offset_so[h] += nsopi_[h ^ hI] * nsopi_[hI];
            offset_mo[h] += navirpi_[h ^ hI] * naoccpi_[hI];
        }
    }

    if (options_.get_str("REFERENCE") != "RHF"){
        // Set up dimensions for b(Q|ia)
        Dimension ov(nirrep_), Q(nirrep_);
        for (int hi = 0; hi < nirrep_; ++hi){
            Q[hi] = nQ_;
            for (int ha = 0; ha < nirrep_; ++ha){
                ov[hi ^ ha] += nboccpi_[hi] * nbvirpi_[ha];
            }
        }
        bQiaB_mo_ = SharedMatrix (new Matrix("b(Q|ia)", Q, ov));

        std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
        for (int h = 0; h < nirrep_; ++h){
            offset_mo.push_back(0);
            offset_so.push_back(0);
        }

        for (int h = 0; h < nirrep_; ++h){
            double** bQmn_so_p = bQmn_so_->pointer(h);
            double** bQiaB_mo_p = bQiaB_mo_->pointer(h);
            for (int hi = 0; hi < nirrep_; ++hi){
                int ha = h ^ hi;
                if (nboccpi_[hi] > 0 && nbvirpi_[ha] > 0){
                    double** Cbvp = Cb_->pointer(ha);
                    double** Cbop = Cb_->pointer(hi);
                    SharedMatrix tmp (new Matrix("Half-transformed b_ov", nQ_, nsopi_[hi] * nbvirpi_[ha]));
                    double** tmpp = tmp->pointer();
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int Q = 0; Q < nQ_; ++Q){
                        // First-half transformation
                        C_DGEMM('N', 'N', nsopi_[hi], nbvirpi_[ha], nsopi_[ha], 1.0, bQmn_so_p[Q] + offset_so[h], nsopi_[ha], Cbvp[0] + nboccpi_[ha], nsopi_[ha], 0.0, tmpp[Q], nbvirpi_[ha]);
                        // Second-half transformation
                        C_DGEMM('T', 'N', nboccpi_[hi], nbvirpi_[ha], nsopi_[hi], 1.0, Cbop[0], nsopi_[hi], tmpp[Q], nbvirpi_[ha], 0.0, bQiaB_mo_p[Q] + offset_mo[h], nbvirpi_[ha]);
                    }
                }
                offset_so[h] += nsopi_[h ^ hi] * nsopi_[hi];
                offset_mo[h] += nbvirpi_[h ^ hi] * nboccpi_[hi];
            }
        }
    }
    dcft_timer_off("DCFTSolver::b(Q|mn) -> b(Q|ia)");
}

/**
  * form b(Q,ab)
  */
void DCFTSolver::formb_vv()
{
    dcft_timer_on("DCFTSolver::b(Q|mn) -> b(Q|ab)");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Set up dimensions for b(Q|AB)
    Dimension VV(nirrep_), Q(nirrep_);
    for (int hA = 0; hA < nirrep_; ++hA){
        Q[hA] = nQ_;
        for (int hB = 0; hB < nirrep_; ++hB){
            VV[hA ^ hB] += navirpi_[hA] * navirpi_[hB];
        }
    }
    bQabA_mo_ = SharedMatrix(new Matrix("b(Q|AB)", Q, VV));

    std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
    for (int h = 0; h < nirrep_; ++h){
        offset_so.push_back(0);
        offset_mo.push_back(0);
    }

    for (int h = 0; h < nirrep_; ++h){
        double** bQmn_so_p = bQmn_so_->pointer(h);
        double** bQabA_mo_p = bQabA_mo_->pointer(h);
        for (int hA = 0; hA < nirrep_; ++hA){
            int hB = h ^ hA;
            if (navirpi_[hA] > 0 && navirpi_[hB] > 0){
                double** CaBp = Ca_->pointer(hB);
                double** CaAp = Ca_->pointer(hA);
                SharedMatrix tmp (new Matrix("Half-transformed b_VV", nQ_, nsopi_[hA] * navirpi_[hB]));
                double** tmpp = tmp->pointer();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int Q = 0; Q < nQ_; ++Q){
                    // First-half transformation
                    C_DGEMM('N', 'N', nsopi_[hA], navirpi_[hB], nsopi_[hB], 1.0, bQmn_so_p[Q] + offset_so[h], nsopi_[hB], CaBp[0] + naoccpi_[hB], nsopi_[hB], 0.0, tmpp[Q], navirpi_[hB]);
                    // Second-half transformation
                    C_DGEMM('T', 'N', navirpi_[hA], navirpi_[hB], nsopi_[hA], 1.0, CaAp[0] + naoccpi_[hA], nsopi_[hA], tmpp[Q], navirpi_[hB], 0.0, bQabA_mo_p[Q] + offset_mo[h], navirpi_[hB]);
                }

            }
            offset_so[h] += nsopi_[h ^ hA] * nsopi_[hA];
            offset_mo[h] += navirpi_[h ^ hA] * navirpi_[hA];
        }
    }

    if (options_.get_str("REFERENCE") != "RHF"){

        // Set up dimensions for b(Q|ab)
        Dimension vv(nirrep_), Q(nirrep_);
        for (int ha = 0; ha < nirrep_; ++ha){
            Q[ha] = nQ_;
            for (int hb = 0; hb < nirrep_; ++hb){
                vv[ha ^ hb] += nbvirpi_[ha] * nbvirpi_[hb];
            }
        }
        bQabB_mo_ = SharedMatrix(new Matrix("b(Q|ab)", Q, vv));

        std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
        for (int h = 0; h < nirrep_; ++h){
            offset_so.push_back(0);
            offset_mo.push_back(0);
        }

        for (int h = 0; h < nirrep_; ++h){
            double** bQmn_so_p = bQmn_so_->pointer(h);
            double** bQabB_mo_p = bQabB_mo_->pointer(h);
            for (int ha = 0; ha < nirrep_; ++ha){
                int hb = h ^ ha;
                if (nbvirpi_[ha] > 0 && nbvirpi_[hb] > 0){
                    double** Cbbp = Cb_->pointer(hb);
                    double** Cbap = Cb_->pointer(ha);
                    SharedMatrix tmp (new Matrix("Half-transformed b_vv", nQ_, nsopi_[ha] * nbvirpi_[hb]));
                    double** tmpp = tmp->pointer();
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int Q = 0; Q < nQ_; ++Q){
                        // First-half transformation
                        C_DGEMM('N', 'N', nsopi_[ha], nbvirpi_[hb], nsopi_[hb], 1.0, bQmn_so_p[Q] + offset_so[h], nsopi_[hb], Cbbp[0] + nboccpi_[hb], nsopi_[hb], 0.0, tmpp[Q], nbvirpi_[hb]);
                        // Second-half transformation
                        C_DGEMM('T', 'N', nbvirpi_[ha], nbvirpi_[hb], nsopi_[ha], 1.0, Cbap[0] + nboccpi_[ha], nsopi_[ha], tmpp[Q], nbvirpi_[hb], 0.0, bQabB_mo_p[Q] + offset_mo[h], nbvirpi_[hb]);
                    }

                }
                offset_so[h] += nsopi_[h ^ ha] * nsopi_[ha];
                offset_mo[h] += nbvirpi_[h ^ ha] * nbvirpi_[ha];
            }
        }
    }

    dcft_timer_off("DCFTSolver::b(Q|mn) -> b(Q|ab)");

}

/**
  * form b(Q,pq)
  */
void DCFTSolver::formb_pq()
{
    dcft_timer_on("DCFTSolver::b(Q|mn) -> b(Q|pq)");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Set up dimensions for b(Aux|PQ)
    Dimension PQ(nirrep_), Aux(nirrep_);
    for (int hP = 0; hP < nirrep_; ++hP){
        Aux[hP] = nQ_;
        for (int hQ = 0; hQ < nirrep_; ++hQ){
            PQ[hP ^ hQ] += nsopi_[hP] * nsopi_[hQ];
        }
    }
    bQpqA_mo_ = SharedMatrix(new Matrix("b(Aux|PQ)", Aux, PQ));

    std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
    for (int h = 0; h < nirrep_; ++h){
        offset_so.push_back(0);
        offset_mo.push_back(0);
    }

    for (int h = 0; h < nirrep_; ++h){
        double** bQmn_so_p = bQmn_so_->pointer(h);
        double** bQpqA_mo_p = bQpqA_mo_->pointer(h);
        for (int hP = 0; hP < nirrep_; ++hP){
            int hQ = h ^ hP;
            if (nsopi_[hP] > 0 && nsopi_[hQ] > 0){
                double** Caqp = Ca_->pointer(hQ);
                double** Capp = Ca_->pointer(hP);
                SharedMatrix tmp (new Matrix("Half-transformed b_PQ", nQ_, nsopi_[hP] * nsopi_[hQ]));
                double** tmpp = tmp->pointer();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int Aux = 0; Aux < nQ_; ++Aux){
                    // First-half transformation
                    C_DGEMM('N', 'N', nsopi_[hP], nsopi_[hQ], nsopi_[hQ], 1.0, bQmn_so_p[Aux] + offset_so[h], nsopi_[hQ], Caqp[0], nsopi_[hQ], 0.0, tmpp[Aux], nsopi_[hQ]);
                    // Second-half transformation
                    C_DGEMM('T', 'N', nsopi_[hP], nsopi_[hQ], nsopi_[hP], 1.0, Capp[0], nsopi_[hP], tmpp[Aux], nsopi_[hQ], 0.0, bQpqA_mo_p[Aux] + offset_mo[h], nsopi_[hQ]);
                }

            }
            offset_so[h] += nsopi_[h ^ hP] * nsopi_[hP];
            offset_mo[h] += nsopi_[h ^ hP] * nsopi_[hP];
        }
    }

    if (options_.get_str("REFERENCE") != "RHF"){
        // Set up dimensions for b(Aux|pq)
        bQpqB_mo_ = SharedMatrix(new Matrix("b(Aux|pq)", Aux, PQ));

        std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
        for (int h = 0; h < nirrep_; ++h){
            offset_so.push_back(0);
            offset_mo.push_back(0);
        }

        for (int h = 0; h < nirrep_; ++h){
            double** bQmn_so_p = bQmn_so_->pointer(h);
            double** bQpqB_mo_p = bQpqB_mo_->pointer(h);
            for (int hp = 0; hp < nirrep_; ++hp){
                int hq = h ^ hp;
                if (nsopi_[hp] > 0 && nsopi_[hq] > 0){
                    double** Cbqp = Cb_->pointer(hq);
                    double** Cbpp = Cb_->pointer(hp);
                    SharedMatrix tmp (new Matrix("Half-transformed b_pq", nQ_, nsopi_[hp] * nsopi_[hq]));
                    double** tmpp = tmp->pointer();
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int Aux = 0; Aux < nQ_; ++Aux){
                        // First-half transformation
                        C_DGEMM('N', 'N', nsopi_[hp], nsopi_[hq], nsopi_[hq], 1.0, bQmn_so_p[Aux] + offset_so[h], nsopi_[hq], Cbqp[0], nsopi_[hq], 0.0, tmpp[Aux], nsopi_[hq]);
                        // Second-half transformation
                        C_DGEMM('T', 'N', nsopi_[hp], nsopi_[hq], nsopi_[hp], 1.0, Cbpp[0], nsopi_[hp], tmpp[Aux], nsopi_[hq], 0.0, bQpqB_mo_p[Aux] + offset_mo[h], nsopi_[hq]);
                    }
                }
                offset_so[h] += nsopi_[h ^ hp] * nsopi_[hp];
                offset_mo[h] += nsopi_[h ^ hp] * nsopi_[hp];
            }
        }
    }

    dcft_timer_off("DCFTSolver::b(Q|mn) -> b(Q|pq)");

}

/**
 * Form density-fitted MO-basis TEI g(OV|OV)
 */
void DCFTSolver::form_df_g_ovov()
{
    dcft_timer_on("DCFTSolver::DF Transform_OVOV");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Alpha-Alpha
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                           ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    for (int h = 0; h < nirrep_; ++h){
        if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
            global_dpd_->buf4_mat_irrep_init(&I, h);
            double** bQiaA_mo_p = bQiaA_mo_->pointer(h);
            C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQiaA_mo_p[0], bQiaA_mo_->coldim(h), bQiaA_mo_p[0], bQiaA_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
    }
    global_dpd_->buf4_close(&I);

    if (options_.get_str("REFERENCE") != "RHF"){
        // Alpha-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                               ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                global_dpd_->buf4_mat_irrep_init(&I, h);
                double** bQiaA_mo_p = bQiaA_mo_->pointer(h);
                double** bQiaB_mo_p = bQiaB_mo_->pointer(h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQiaA_mo_p[0], bQiaA_mo_->coldim(h), bQiaB_mo_p[0], bQiaB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

        // Beta-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                               ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                global_dpd_->buf4_mat_irrep_init(&I, h);
                double** bQiaB_mo_p = bQiaB_mo_->pointer(h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQiaB_mo_p[0], bQiaB_mo_->coldim(h), bQiaB_mo_p[0], bQiaB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);
    }

    dcft_timer_off("DCFTSolver::DF Transform_OVOV");


}

/**
 * Form density-fitted MO-basis TEI g(OO|OO)
 */
void DCFTSolver::form_df_g_oooo()
{
    dcft_timer_on("DCFTSolver::DF Transform_OOOO");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Alpha-Alpha
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                           ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
    for (int h = 0; h < nirrep_; ++h){
        if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
            double** bQijA_mo_p = bQijA_mo_->pointer(h);
            global_dpd_->buf4_mat_irrep_init(&I, h);
            C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQijA_mo_p[0], bQijA_mo_->coldim(h), bQijA_mo_p[0], bQijA_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
    }
    global_dpd_->buf4_close(&I);

    if (options_.get_str("REFERENCE") != "RHF"){
        // Alpha-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                               ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQijA_mo_p = bQijA_mo_->pointer(h);
                double** bQijB_mo_p = bQijB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQijA_mo_p[0], bQijA_mo_->coldim(h), bQijB_mo_p[0], bQijB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I ,h);
            }
        }
        global_dpd_->buf4_close(&I);

        // Beta-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                      ID("[o>=o]+"), ID("[o>=o]+"), 0, "MO Ints (oo|oo)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQijB_mo_p = bQijB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQijB_mo_p[0], bQijB_mo_->coldim(h), bQijB_mo_p[0], bQijB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

    }

    dcft_timer_off("DCFTSolver::DF Transform_OOOO");

}

/**
 * Form density-fitted MO-basis TEI g(VV|OO)
 */
void DCFTSolver::form_df_g_vvoo()
{
    dcft_timer_on("DCFTSolver::DF Transform_OOVV");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    dpdbuf4 I;

    if (options_.get_str("REFERENCE") == "RHF"){
        // g(AB|IJ) = Sum_Q b(AB|Q) b(Q|IJ)
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                               ID("[V>=V]+"), ID("[O>=O]+"), 0, "MO Ints (VV|OO)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQabA_mo_p = bQabA_mo_->pointer(h);
                double** bQijA_mo_p = bQijA_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQabA_mo_p[0], bQabA_mo_->coldim(h), bQijA_mo_p[0], bQijA_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

    }
    else{
        // g(ab|ij) = Sum_Q b(ab|Q) b(Q|ij)

        // Alpha-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                               ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQabA_mo_p = bQabA_mo_->pointer(h);
                double** bQijB_mo_p = bQijB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQabA_mo_p[0], bQabA_mo_->coldim(h), bQijB_mo_p[0], bQijB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

        // g(ij|ab) = Sum_Q b(ij|Q) b(Q|ab)

        // Alpha-Alpha
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                               ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQijA_mo_p = bQijA_mo_->pointer(h);
                double** bQabA_mo_p = bQabA_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQijA_mo_p[0], bQijA_mo_->coldim(h), bQabA_mo_p[0], bQabA_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

        // Alpha-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                               ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQijA_mo_p = bQijA_mo_->pointer(h);
                double** bQabB_mo_p = bQabB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQijA_mo_p[0], bQijA_mo_->coldim(h), bQabB_mo_p[0], bQabB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

        // Beta-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                               ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQijB_mo_p = bQijB_mo_->pointer(h);
                double** bQabB_mo_p = bQabB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQijB_mo_p[0], bQijB_mo_->coldim(h), bQabB_mo_p[0], bQabB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

    }

    dcft_timer_off("DCFTSolver::DF Transform_OOVV");

}

/**
 * Form density-fitted MO-basis TEI g(VO|OO)
 */
void DCFTSolver::form_df_g_vooo()
{
    dcft_timer_on("DCFTSolver::DF Transform_VOOO");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    dpdbuf4 I;

    /*** Form b(Q|AI) ***/

    // Put detailed information of b(Q|ai) block into 'block_Qai'
    // Put detailed information of b(Q|ia) block into 'block_Qia'
    std::vector<std::vector<std::pair<long int, long int>>> block_Qai, block_Qia;
    Dimension VO(nirrep_), Q(nirrep_);
    for (int h = 0; h < nirrep_; ++h){
        long int entrance_Qai = 0;
        long int entrance_Qia = 0;
        std::vector<std::pair<long int, long int>> subblock_Qai;
        std::vector<std::pair<long int, long int>> subblock_Qia;
        // b(Q|ai) subblocks
        for (int ha = 0; ha < nirrep_; ++ha){
            int hi = h ^ ha;
            std::pair<long int, long int> subsubblock(entrance_Qai, navirpi_[ha] * naoccpi_[hi]);
            subblock_Qai.push_back(subsubblock);
            entrance_Qai += subsubblock.second;
        }
        block_Qai.push_back(subblock_Qai);
        // b(Q|ia) subblocks
        for (int hi = 0; hi < nirrep_; ++hi){
            int ha = h ^ hi;
            std::pair<long int, long int> subsubblock(entrance_Qia, naoccpi_[hi] * navirpi_[ha]);
            subblock_Qia.push_back(subsubblock);
            entrance_Qia += subsubblock.second;
        }
        block_Qia.push_back(subblock_Qia);
        // Dimension of b(Q|ai)
        Q[h] = nQ_;
        VO[h] = entrance_Qai;
    }

    // Sort b(Q|IA) -> b(Q|AI)
    bQaiA_mo_ = SharedMatrix(new Matrix("b(Q|AI)", Q, VO));
    for (int h = 0; h < nirrep_; ++h){
        for (int hA = 0; hA < nirrep_; ++hA){
            int hI = h ^ hA;
            if (navirpi_[hA] > 0 && naoccpi_[hI] > 0){
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int A = 0; A < navirpi_[hA]; ++A){
                    for (int I = 0; I < naoccpi_[hI]; ++I){
                        long int IA = block_Qia[h][hI].first + I * navirpi_[hA] + A;
                        long int AI = block_Qai[h][hA].first + A * naoccpi_[hI] + I;
                        bQaiA_mo_->set_column(h, AI, bQiaA_mo_->get_column(h, IA));
                    }
                }
            }
        }
    }

    // g(ai|jk) = Sum_Q b(ai|Q) (Q|jk)

    // Alpha-Alpha
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                           ID("[V,O]"), ID("[O>=O]+"), 0, "MO Ints (VO|OO)");
    for (int h = 0; h < nirrep_; ++h){
        if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
            double** bQaiA_mo_p = bQaiA_mo_->pointer(h);
            double** bQijA_mo_p = bQijA_mo_->pointer(h);
            global_dpd_->buf4_mat_irrep_init(&I, h);
            C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQaiA_mo_p[0], bQaiA_mo_->coldim(h), bQijA_mo_p[0], bQijA_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
    }
    global_dpd_->buf4_close(&I);

    if (options_.get_str("REFERENCE") != "RHF"){

        /*** Form b(Q|ai) ***/

        // Put detailed information of b(Q|ai) block into 'block_Qai'
        // Put detailed information of b(Q|ia) block into 'block_Qia'
        std::vector<std::vector<std::pair<long int, long int>>> block_Qai, block_Qia;
        Dimension vo(nirrep_), Q(nirrep_);
        for (int h = 0; h < nirrep_; ++h){
            long int entrance_Qai = 0;
            long int entrance_Qia = 0;
            std::vector<std::pair<long int, long int>> subblock_Qai;
            std::vector<std::pair<long int, long int>> subblock_Qia;
            // b(Q|ai) subblocks
            for (int ha = 0; ha < nirrep_; ++ha){
                int hi = h ^ ha;
                std::pair<long int, long int> subsubblock(entrance_Qai, nbvirpi_[ha] * nboccpi_[hi]);
                subblock_Qai.push_back(subsubblock);
                entrance_Qai += subsubblock.second;
            }
            block_Qai.push_back(subblock_Qai);
            // b(Q|ia) subblocks
            for (int hi = 0; hi < nirrep_; ++hi){
                int ha = h ^ hi;
                std::pair<long int, long int> subsubblock(entrance_Qia, nboccpi_[hi] * nbvirpi_[ha]);
                subblock_Qia.push_back(subsubblock);
                entrance_Qia += subsubblock.second;
            }
            block_Qia.push_back(subblock_Qia);
            // Dimension of b(Q|ai)
            Q[h] = nQ_;
            vo[h] = entrance_Qai;
        }

        // Sort b(Q|ia) -> b(Q|ai)
        bQaiB_mo_ = SharedMatrix(new Matrix("b(Q|ai)", Q, vo));
        for (int h = 0; h < nirrep_; ++h){
            for (int ha = 0; ha < nirrep_; ++ha){
                int hi = h ^ ha;
                if (nbvirpi_[ha] > 0 && nboccpi_[hi] > 0){
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int a = 0; a < nbvirpi_[ha]; ++a){
                        for (int i = 0; i < nboccpi_[hi]; ++i){
                            long int ia = block_Qia[h][hi].first + i * nbvirpi_[ha] + a;
                            long int ai = block_Qai[h][ha].first + a * nboccpi_[hi] + i;
                            bQaiB_mo_->set_column(h, ai, bQiaB_mo_->get_column(h, ia));
                        }
                    }
                }
            }
        }


        // g(ai|jk) = Sum_Q b(ai|Q) (Q|jk)

        // Alpha-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,o]"),
                               ID("[V,O]"), ID("[o>=o]+"), 0, "MO Ints (VO|oo)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQaiA_mo_p = bQaiA_mo_->pointer(h);
                double** bQijB_mo_p = bQijB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQaiA_mo_p[0], bQaiA_mo_->coldim(h), bQijB_mo_p[0], bQijB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

        // Beta-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                               ID("[v,o]"), ID("[o>=o]+"), 0, "MO Ints (vo|oo)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQaiB_mo_p = bQaiB_mo_->pointer(h);
                double** bQijB_mo_p = bQijB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQaiB_mo_p[0], bQaiB_mo_->coldim(h), bQijB_mo_p[0], bQijB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

        // g(jk|ai) = Sum_Q b(jk|Q) (Q|ai)

        // Alpha-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,o]"),
                               ID("[O>=O]+"), ID("[v,o]"), 0, "MO Ints (OO|vo)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQijA_mo_p = bQijA_mo_->pointer(h);
                double** bQaiB_mo_p = bQaiB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQijA_mo_p[0], bQijA_mo_->coldim(h), bQaiB_mo_p[0], bQaiB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);

            }
        }
    }

    dcft_timer_off("DCFTSolver::DF Transform_VOOO");

}

/**
 * Form density-fitted MO-basis TEI g(OV|VV)
 */
void DCFTSolver::form_df_g_ovvv()
{
    dcft_timer_on("DCFTSolver::DF Transform_OVVV");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    dpdbuf4 I;

    // g(ia|bc) = Sum_Q b(ia|Q) (Q|bc)

    // Alpha-Alpha
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                           ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    for (int h = 0; h < nirrep_; ++h){
        if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
            double** bQiaA_mo_p = bQiaA_mo_->pointer(h);
            double** bQabA_mo_p = bQabA_mo_->pointer(h);
            global_dpd_->buf4_mat_irrep_init(&I, h);
            C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQiaA_mo_p[0], bQiaA_mo_->coldim(h), bQabA_mo_p[0], bQabA_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
    }
    global_dpd_->buf4_close(&I);

    if (options_.get_str("REFERENCE") != "RHF"){
        // g(ia|bc) = Sum_Q b(ia|Q) (Q|bc)

        // Alpha-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                               ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQiaA_mo_p = bQiaA_mo_->pointer(h);
                double** bQabB_mo_p = bQabB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQiaA_mo_p[0], bQiaA_mo_->coldim(h), bQabB_mo_p[0], bQabB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

        // Beta-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                               ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQiaB_mo_p = bQiaB_mo_->pointer(h);
                double** bQabB_mo_p = bQabB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQiaB_mo_p[0], bQiaB_mo_->coldim(h), bQabB_mo_p[0], bQabB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

        // g(bc|ia) = Sum_Q b(bc|Q) (Q|ia)

        // Alpha-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,v]"),
                               ID("[V>=V]+"), ID("[o,v]"), 0, "MO Ints (VV|ov)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQabA_mo_p = bQabA_mo_->pointer(h);
                double** bQiaB_mo_p = bQiaB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQabA_mo_p[0], bQabA_mo_->coldim(h), bQiaB_mo_p[0], bQiaB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);
    }

    dcft_timer_off("DCFTSolver::DF Transform_OVVV");

}

/**
 * Form density-fitted MO-basis TEI g(VV|VV)
 */
void DCFTSolver::form_df_g_vvvv()
{
    dcft_timer_on("DCFTSolver::DF Transform_VVVV");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    dpdbuf4 I;

    // g(ab|cd) = Sum_Q b(ab|Q) b(Q|cd)
    // Alpha-Alpha
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                           ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
    for (int h = 0; h < nirrep_; ++h){
        if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
            double** bQabA_mo_p = bQabA_mo_->pointer(h);
            global_dpd_->buf4_mat_irrep_init(&I, h);
            C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQabA_mo_p[0], bQabA_mo_->coldim(h), bQabA_mo_p[0], bQabA_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
            global_dpd_->buf4_mat_irrep_wrt(&I, h);
            global_dpd_->buf4_mat_irrep_close(&I, h);
        }
    }
    global_dpd_->buf4_close(&I);

    if (options_.get_str("REFERENCE") != "RHF"){
        // Alpha-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                               ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQabA_mo_p = bQabA_mo_->pointer(h);
                double** bQabB_mo_p = bQabB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQabA_mo_p[0], bQabA_mo_->coldim(h), bQabB_mo_p[0], bQabB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);

        // Beta-Beta
        global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                               ID("[v>=v]+"), ID("[v>=v]+"), 0, "MO Ints (vv|vv)");
        for (int h = 0; h < nirrep_; ++h){
            if (I.params->rowtot[h] > 0 && I.params->coltot[h] > 0){
                double** bQabB_mo_p = bQabB_mo_->pointer(h);
                global_dpd_->buf4_mat_irrep_init(&I, h);
                C_DGEMM('T', 'N', I.params->rowtot[h], I.params->coltot[h], nQ_, 1.0, bQabB_mo_p[0], bQabB_mo_->coldim(h), bQabB_mo_p[0], bQabB_mo_->coldim(h), 0.0, I.matrix[h][0], I.params->coltot[h]);
                global_dpd_->buf4_mat_irrep_wrt(&I, h);
                global_dpd_->buf4_mat_irrep_close(&I, h);
            }
        }
        global_dpd_->buf4_close(&I);
    }

    dcft_timer_off("DCFTSolver::DF Transform_VVVV");

}

/**
 * Compute the density-fitted ERI <vv||vv> tensors in G intermediates
 * and contract with lambda_ijcd.
 * Compute the density-fitted ERI <qs|pr> tensors
 * and contract with gamma<r|s>
 */
void DCFTSolver::build_DF_tensors_RHF()
{
    dcft_timer_on("DCFTSolver::build_df_tensors_RHF()");

    // Form gbar<AB|CD> lambda<CD|IJ>
    build_gbarlambda_RHF_v3mem();

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

    /* Build [Gbar*Gamma]<Q|P> */
    build_gbarGamma_RHF();

    dcft_timer_off("DCFTSolver::build_df_tensors_RHF()");
}

/**
  * Compute the contraction, gbar<ab|cd> lambda<ij|cd>, using density fitting.
  * Memory required: O(V^3)
  */
void DCFTSolver::build_gbarlambda_RHF_v3mem()
{
    dcft_timer_on("DCFTSolver::DF lambda<ij|cd> gbar<ab|cd> (v3 in memory)");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Put detailed information of b(Q|ab) block into 'block'
    std::vector<std::vector<std::pair<long int, long int>>> block;
    for (int hab = 0; hab < nirrep_; ++hab){
        long int entrance = 0;
        std::vector<std::pair<long int, long int>> subblock;
        for (int ha = 0; ha < nirrep_; ++ha){
            int hb = hab ^ ha;
            std::pair<long int, long int> subsubblock(entrance, navirpi_[ha] * navirpi_[hb]);
            subblock.push_back(subsubblock);
            entrance += subsubblock.second;
        }
        block.push_back(subblock);
    }

    /*
     * Intermediate G_SF_<IJ|AB> = lambda_SF_<IJ|CD> g<AB|CD>
     */

    dpdbuf4 Laa, Gaa;

    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "Lambda SF <OO|VV>");
    global_dpd_->buf4_init(&Gaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) SF <OO|VV>");
    global_dpd_->buf4_scm(&Gaa, 0.0);

    for (int hac = 0; hac < nirrep_; ++hac){
        for (int ha = 0; ha < nirrep_; ++ha){
            int hc = hac ^ ha;
            int hbd = hac;
            for (int hb = 0; hb < nirrep_; ++hb){
                int hd = hbd ^ hb;
                int hab = ha ^ hb;
                int hcd = hc ^ hd;
                int hij = hcd;

                if (Laa.params->rowtot[hij] > 0 && Laa.params->coltot[hcd] > 0 && Gaa.params->rowtot[hij] > 0 && Gaa.params->coltot[hab] > 0
                        && navirpi_[ha] > 0 && navirpi_[hc] > 0 && navirpi_[hb] > 0 && navirpi_[hd] > 0){
                    double** bQvvAp = bQabA_mo_->pointer(hac);

                    global_dpd_->buf4_mat_irrep_init(&Laa, hij);
                    global_dpd_->buf4_mat_irrep_rd(&Laa, hij);
                    global_dpd_->buf4_mat_irrep_init(&Gaa, hij);
                    global_dpd_->buf4_mat_irrep_rd(&Gaa, hij);

                    if (hb == hd){
                        std::vector<SharedMatrix> CBD;
                        for (int i = 0; i < nthreads; ++i){
                            CBD.push_back(SharedMatrix(new Matrix("g(A'C|BD)", navirpi_[hc], navirpi_[hb]*navirpi_[hd])));
                        }
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int A = 0; A < navirpi_[ha]; ++A){
                            int thread = 0;
                            #ifdef _OPENMP
                                 thread = omp_get_thread_num();
                            #endif
                            double** CBDp = CBD[thread]->pointer();

                            // g(A'C|BD) = b(A'C|Q) b(Q|BD)
                            C_DGEMM('T', 'N', navirpi_[hc], navirpi_[hb]*navirpi_[hd], nQ_, 1.0, bQvvAp[0]+block[hac][ha].first+A*navirpi_[hc], bQabA_mo_->coldim(hac), bQvvAp[0]+block[hbd][hb].first, bQabA_mo_->coldim(hbd), 0.0, CBDp[0], navirpi_[hb]*navirpi_[hd]);
                            // G<IJ|A'B> = lambda<IJ|CD> g(A'C|DB)
                            C_DGEMM('N', 'N', Gaa.params->rowtot[hij], navirpi_[hb], navirpi_[hc]*navirpi_[hd], 1.0, Laa.matrix[hij][0]+block[hcd][hc].first, Laa.params->coltot[hij], CBDp[0], navirpi_[hb], 1.0, Gaa.matrix[hij][0]+block[hab][ha].first+A*navirpi_[hb], Gaa.params->coltot[hij]);
                        }
                    }
                    else{
                        std::vector<SharedMatrix> CBD, CDB;
                        for (int i = 0; i < nthreads; ++i){
                            CBD.push_back(SharedMatrix(new Matrix("g(A'C|BD)", navirpi_[hc], navirpi_[hb]*navirpi_[hd])));
                            CDB.push_back(SharedMatrix(new Matrix("g(A'C|DB)", navirpi_[hc], navirpi_[hd]*navirpi_[hb])));
                        }
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int A = 0; A < navirpi_[ha]; ++A){
                            int thread = 0;
                            #ifdef _OPENMP
                                 thread = omp_get_thread_num();
                            #endif
                            double** CBDp = CBD[thread]->pointer();

                            // g(A'C|BD) = b(A'C|Q) b(Q|BD)
                            C_DGEMM('T', 'N', navirpi_[hc], navirpi_[hb]*navirpi_[hd], nQ_, 1.0, bQvvAp[0]+block[hac][ha].first+A*navirpi_[hc], bQabA_mo_->coldim(hac), bQvvAp[0]+block[hbd][hb].first, bQabA_mo_->coldim(hbd), 0.0, CBDp[0], navirpi_[hb]*navirpi_[hd]);

                            // g(A'C|BD) -> g(A'C|DB)
                            for (int B = 0; B < navirpi_[hb]; ++B){
                                for (int D = 0; D < navirpi_[hd]; ++D)
                                    CDB[thread]->set_column(0, D*navirpi_[hb]+B, CBD[thread]->get_column(0, B*navirpi_[hd]+D));
                            }

                            double** CDBp = CDB[thread]->pointer();

                            // G<IJ|A'B> = lambda<IJ|CD> g(A'C|DB)
                            C_DGEMM('N', 'N', Gaa.params->rowtot[hij], navirpi_[hb], navirpi_[hc]*navirpi_[hd], 1.0, Laa.matrix[hij][0]+block[hcd][hc].first, Laa.params->coltot[hij], CDBp[0], navirpi_[hb], 1.0, Gaa.matrix[hij][0]+block[hab][ha].first+A*navirpi_[hb], Gaa.params->coltot[hij]);
                        }

                    }
                    global_dpd_->buf4_mat_irrep_wrt(&Gaa, hij);
                    global_dpd_->buf4_mat_irrep_close(&Gaa, hij);

                    global_dpd_->buf4_mat_irrep_close(&Laa, hij);
                }
            }
        }
    }

    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Gaa);

    dcft_timer_off("DCFTSolver::DF lambda<ij|cd> gbar<ab|cd> (v3 in memory)");
}

/**
 * Form MO-based contraction [Gbar*Gamma]<q|p>
 * [Gbar*Gamma]<q|p> = Sum_rs Gbar<qs|pr> Gamma<r|s>
 */
void DCFTSolver::build_gbarGamma_RHF()
{
    dcft_timer_on("DCFTSolver::Gbar<QS|PR> Gamma<R|S> (FastBuilder)");

    build_gbarKappa_RHF();

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Form gamma<R|S> = kappa<R|S> + tau<R|S>
    mo_gammaA_ = SharedMatrix (new Matrix("MO-basis Gamma", nirrep_, nmopi_, nmopi_));
//    mo_gammaA_->copy(kappa_mo_a_);
//    mo_gammaA_->add(mo_tauA_);
    mo_gbarGamma_A_ = SharedMatrix (new Matrix("MO-basis Gbar*Gamma", nirrep_, nmopi_, nmopi_));
    mo_gammaA_->copy(mo_tauA_);

    // Put detailed information of b(Q|pq) block into 'block'
    std::vector<std::vector<std::pair<long int, long int>>> block;
    for (int hpq = 0; hpq < nirrep_; ++hpq){
        long int entrance = 0;
        std::vector<std::pair<long int, long int>> subblock;
        for (int hp = 0; hp < nirrep_; ++hp){
            int hq = hpq ^ hp;
            std::pair<long int, long int> subsubblock(entrance, nsopi_[hp] * nsopi_[hq]);
            subblock.push_back(subsubblock);
            entrance += subsubblock.second;
        }
        block.push_back(subblock);
    }

    /*
     *  f_tilde <Q|P> = gbar<QS|PR> gamma<R|S> + gbar<Qs|Pr> gamma<r|s>
     *                  = 2 g(QP|SR) gamma<R|S> - g(QR|SP) gamma<R|S>
     *                  = 2 b(QP|Aux) b(Aux|SR) gamma<R|S> - b(QR|Aux) b(Aux|SP) gamma<R|S>
     */

    // f_tilde <Q|P> = 2 b(QP|Aux) b(Aux|SR) gamma<R|S>
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int hq = 0; hq < nirrep_; ++hq){
        int hp = hq;
        if (nsopi_[hq] > 0){
            double** tFAp = mo_gbarGamma_A_->pointer(hq);
            double** bQpqAp = bQpqA_mo_->pointer(0);
            SharedMatrix Q (new Matrix("b(Q|SR)gamma<R|S>", 1, nQ_));
            double** Qp = Q->pointer();
            // (Q) = b(Q|SR) gamma<R|S>
            for (int hr = 0; hr < nirrep_; ++hr){
                int hs = hr;
                if (nsopi_[hr] > 0){
                    double** gamma_rs_p = mo_gammaA_->pointer(hr);
                    C_DGEMV('N', nQ_, nsopi_[hr]*nsopi_[hs], 1.0, bQpqAp[0]+block[0][hr].first, bQpqA_mo_->coldim(0), gamma_rs_p[0], 1, 1.0, Qp[0], 1);
                }
            }
            // tilde_f <Q|P> = 2 b(QP|Aux)*(Aux) where (Aux) = (Q)
            C_DGEMV('T', nQ_, nsopi_[hp]*nsopi_[hq], 2.0, bQpqAp[0]+block[0][hp].first, bQpqA_mo_->coldim(0), Qp[0], 1, 0.0, tFAp[0], 1);
        }
    }

    // f_tilde <Q|P> -= b(QR|Aux) b(Aux|SP) gamma<R|S>
    for (int hq = 0; hq < nirrep_; ++hq){
        int hp = hq;
        if (nsopi_[hq] > 0){
            for (int hr = 0; hr < nirrep_; ++hr){
                int hs = hr;
                if (nsopi_[hr] > 0){
                    double** bQpqAp = bQpqA_mo_->pointer(hq^hr);
                    double** gamma_rs_p = mo_gammaA_->pointer(hr);

                    std::vector<SharedMatrix> rs;
                    for (int i = 0; i < nthreads; ++i){
                        rs.push_back(SharedMatrix(new Matrix("<Q'P'|RS>", nsopi_[hr], nsopi_[hs])));
                    }

                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int q = 0; q < nsopi_[hq]; ++q){
                        for (int p = q; p < nsopi_[hp]; ++p){

                            int thread = 0;
                            #ifdef _OPENMP
                            thread = omp_get_thread_num();
                            #endif
                            double** rsp = rs[thread]->pointer();

                            // <Q'P'|RS> = b(Q'R|Aux) b(Aux|P'S)
                            C_DGEMM('T', 'N', nsopi_[hr], nsopi_[hs], nQ_, 1.0, bQpqAp[0]+block[hq^hr][hq].first+q*nsopi_[hr], bQpqA_mo_->coldim(hq^hr), bQpqAp[0]+block[hp^hs][hp].first+p*nsopi_[hs], bQpqA_mo_->coldim(hp^hs), 0.0, rsp[0], nsopi_[hs]);
                            // - <Q'P'|RS> * gamma<R|S>
                            double value = - C_DDOT(nsopi_[hr]*nsopi_[hs], rsp[0], 1, gamma_rs_p[0], 1);
                            mo_gbarGamma_A_->add(hp, q, p, value);
                            if (q != p){
                                mo_gbarGamma_A_->add(hp, p, q, value);
                            }
                        }
                    }
                }
            }
        }
    }

    mo_gbarGamma_A_->add(mo_gbarKappa_A_);

    dcft_timer_off("DCFTSolver::Gbar<QS|PR> Gamma<R|S> (FastBuilder)");

}

void DCFTSolver::build_gbarKappa_RHF()
{
    dcft_timer_on("DCFTSolver::Gbar<QS|PR> Kappa<R|S>");

    formb_pq_scf();

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    mo_gbarKappa_A_ = SharedMatrix (new Matrix("MO-basis Gbar*Kappa", nirrep_, nmopi_, nmopi_));

    // Put detailed information of b(Q|pq) block into 'block'
    std::vector<std::vector<std::pair<long int, long int>>> block;
    for (int hpq = 0; hpq < nirrep_; ++hpq){
        long int entrance = 0;
        std::vector<std::pair<long int, long int>> subblock;
        for (int hp = 0; hp < nirrep_; ++hp){
            int hq = hpq ^ hp;
            std::pair<long int, long int> subsubblock(entrance, nsopi_[hp] * nsopi_[hq]);
            subblock.push_back(subsubblock);
            entrance += subsubblock.second;
        }
        block.push_back(subblock);
    }

    /*
     *  f_tilde <Q|P> = gbar<QS|PR> gamma<R|S> + gbar<Qs|Pr> gamma<r|s>
     *                  = 2 g(QP|SR) gamma<R|S> - g(QR|SP) gamma<R|S>
     *                  = 2 b(QP|Aux) b(Aux|SR) gamma<R|S> - b(QR|Aux) b(Aux|SP) gamma<R|S>
     */

    // f_tilde <Q|P> = 2 b(QP|Aux) b(Aux|SR) gamma<R|S>
    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int hq = 0; hq < nirrep_; ++hq){
        int hp = hq;
        if (nsopi_[hq] > 0){
            double** tFAp = mo_gbarKappa_A_->pointer(hq);
            double** bQpqAp = bQpqA_mo_scf_->pointer(0);
            SharedMatrix Q (new Matrix("b(Q|SR)gamma<R|S>", 1, nQ_scf_));
            double** Qp = Q->pointer();
            // (Q) = b(Q|SR) gamma<R|S>
            for (int hr = 0; hr < nirrep_; ++hr){
                int hs = hr;
                if (nsopi_[hr] > 0){
                    double** gamma_rs_p = kappa_mo_a_->pointer(hr);
                    C_DGEMV('N', nQ_scf_, nsopi_[hr]*nsopi_[hs], 1.0, bQpqAp[0]+block[0][hr].first, bQpqA_mo_scf_->coldim(0), gamma_rs_p[0], 1, 1.0, Qp[0], 1);
                }
            }
            // tilde_f <Q|P> = 2 b(QP|Aux)*(Aux) where (Aux) = (Q)
            C_DGEMV('T', nQ_scf_, nsopi_[hp]*nsopi_[hq], 2.0, bQpqAp[0]+block[0][hp].first, bQpqA_mo_scf_->coldim(0), Qp[0], 1, 0.0, tFAp[0], 1);
        }
    }

    // f_tilde <Q|P> -= b(QR|Aux) b(Aux|SP) gamma<R|S>
    for (int hq = 0; hq < nirrep_; ++hq){
        int hp = hq;
        if (nsopi_[hq] > 0){
            for (int hr = 0; hr < nirrep_; ++hr){
                int hs = hr;
                if (nsopi_[hr] > 0){
                    double** bQpqAp = bQpqA_mo_scf_->pointer(hq^hr);
                    double** gamma_rs_p = kappa_mo_a_->pointer(hr);

                    std::vector<SharedMatrix> rs;
                    for (int i = 0; i < nthreads; ++i){
                        rs.push_back(SharedMatrix(new Matrix("<Q'P'|RS>", nsopi_[hr], nsopi_[hs])));
                    }

                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int q = 0; q < nsopi_[hq]; ++q){
                        for (int p = q; p < nsopi_[hp]; ++p){

                            int thread = 0;
                            #ifdef _OPENMP
                            thread = omp_get_thread_num();
                            #endif
                            double** rsp = rs[thread]->pointer();

                            // <Q'P'|RS> = b(Q'R|Aux) b(Aux|P'S)
                            C_DGEMM('T', 'N', nsopi_[hr], nsopi_[hs], nQ_scf_, 1.0, bQpqAp[0]+block[hq^hr][hq].first+q*nsopi_[hr], bQpqA_mo_scf_->coldim(hq^hr), bQpqAp[0]+block[hp^hs][hp].first+p*nsopi_[hs], bQpqA_mo_scf_->coldim(hp^hs), 0.0, rsp[0], nsopi_[hs]);
                            // - <Q'P'|RS> * gamma<R|S>
                            double value = - C_DDOT(nsopi_[hr]*nsopi_[hs], rsp[0], 1, gamma_rs_p[0], 1);
                            mo_gbarKappa_A_->add(hp, q, p, value);
                            if (q != p){
                                mo_gbarKappa_A_->add(hp, p, q, value);
                            }
                        }
                    }
                }
            }
        }
    }

    bQpqA_mo_scf_.reset();

    dcft_timer_off("DCFTSolver::Gbar<QS|PR> Kappa<R|S>");

}


/**
 * Compute the density-fitted ERI <vv||vv> tensors in G intermediates
 * and contract with lambda_ijcd.
 * Compute the density-fitted ERI <qs|pr> tensors
 * and contract with gamma<r|s>
 */
void DCFTSolver::build_DF_tensors_UHF()
{
    dcft_timer_on("DCFTSolver::build_df_tensors_UHF");

    // Form gbar<AB|CD> lambda<CD|IJ>
    build_gbarlambda_UHF_v3mem();

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

    /* Build [gbar*gamma]<q|p> */
    build_gbarGamma_UHF();

    dcft_timer_off("DCFTSolver::build_df_tensors_UHF");

}

/**
  * Compute the contraction, gbar<ab|cd> lambda<ij|cd>, using density fitting.
  * Memory required: O(V^3)
  */
void DCFTSolver::build_gbarlambda_UHF_v3mem()
{
    dcft_timer_on("DCFTSolver::DF lambda<ij|cd> gbar<ab|cd> (v3 in memory)");

    // Thread considerations
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    /********** Alpha-Alpha **********/

    // Put detailed information of b(Q|AB) block into 'block_AB'
    std::vector<std::vector<std::pair<long int, long int>>> block_AB;
    for (int hAB = 0; hAB < nirrep_; ++hAB){
        long int entrance = 0;
        std::vector<std::pair<long int, long int>> subblock;
        for (int hA = 0; hA < nirrep_; ++hA){
            int hB = hAB ^ hA;
            std::pair<long int, long int> subsubblock(entrance, navirpi_[hA] * navirpi_[hB]);
            subblock.push_back(subsubblock);
            entrance += subsubblock.second;
        }
        block_AB.push_back(subblock);
    }

    /*
     * Intermediate G <IJ|AB> = 1/2 lambda<IJ|CD> gbar<AB|CD>
     *                        = 1/2 lambda<IJ|CD> g(AC|BD) - 1/2 lambda<IJ|CD> g(AD|BC)
     *                        = lambda<IJ|CD> g(AC|BD)
     */

    dpdbuf4 Laa, Gaa;

    global_dpd_->buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    global_dpd_->buf4_init(&Gaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                           ID("[O,O]"), ID("[V,V]"), 0, "tau(temp) <OO|VV>");
    global_dpd_->buf4_scm(&Gaa, 0.0);

    for (int hAC = 0; hAC < nirrep_; ++hAC){
        for (int hA = 0; hA < nirrep_; ++hA){
            int hC = hAC ^ hA;
            int hBD = hAC;
            for (int hB = 0; hB < nirrep_; ++hB){
                int hD = hBD ^ hB;
                int hAB = hA ^ hB;
                int hCD = hC ^ hD;
                int hIJ = hCD;

                if (Laa.params->rowtot[hIJ] > 0 && Laa.params->coltot[hCD] > 0 && Gaa.params->rowtot[hIJ] > 0 && Gaa.params->coltot[hAB] > 0
                        && navirpi_[hA] > 0 && navirpi_[hC] > 0 && navirpi_[hB] > 0 && navirpi_[hD] > 0){
                    double** bQvvAp = bQabA_mo_->pointer(hAC);

                    global_dpd_->buf4_mat_irrep_init(&Laa, hIJ);
                    global_dpd_->buf4_mat_irrep_rd(&Laa, hIJ);
                    global_dpd_->buf4_mat_irrep_init(&Gaa, hIJ);
                    global_dpd_->buf4_mat_irrep_rd(&Gaa, hIJ);

                    if (hB == hD){
                        std::vector<SharedMatrix> CBD;
                        for (int i = 0; i < nthreads; ++i){
                            CBD.push_back(SharedMatrix(new Matrix("g(A'C|BD)", navirpi_[hC], navirpi_[hB]*navirpi_[hD])));
                        }
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int A = 0; A < navirpi_[hA]; ++A){
                            int thread = 0;
                            #ifdef _OPENMP
                                 thread = omp_get_thread_num();
                            #endif
                            double** CBDp = CBD[thread]->pointer();
                            // g(A'C|BD) = b(A'C|Q) b(Q|BD)
                            C_DGEMM('T', 'N', navirpi_[hC], navirpi_[hB]*navirpi_[hD], nQ_, 1.0, bQvvAp[0]+block_AB[hAC][hA].first+A*navirpi_[hC], bQabA_mo_->coldim(hAC), bQvvAp[0]+block_AB[hBD][hB].first, bQabA_mo_->coldim(hBD), 0.0, CBDp[0], navirpi_[hB]*navirpi_[hD]);
                            // G<IJ|A'B> = lambda<IJ|CD> g(A'C|DB)
                            C_DGEMM('N', 'N', Gaa.params->rowtot[hIJ], navirpi_[hB], navirpi_[hC]*navirpi_[hD], 1.0, Laa.matrix[hIJ][0]+block_AB[hCD][hC].first, Laa.params->coltot[hIJ], CBDp[0], navirpi_[hB], 1.0, Gaa.matrix[hIJ][0]+block_AB[hAB][hA].first+A*navirpi_[hB], Gaa.params->coltot[hIJ]);
                        }
                    }
                    else{
                        std::vector<SharedMatrix> CBD, CDB;
                        for (int i = 0; i < nthreads; ++i){
                            CBD.push_back(SharedMatrix(new Matrix("g(A'C|BD)", navirpi_[hC], navirpi_[hB]*navirpi_[hD])));
                            CDB.push_back(SharedMatrix(new Matrix("g(A'C|DB)", navirpi_[hC], navirpi_[hD]*navirpi_[hB])));
                        }
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int A = 0; A < navirpi_[hA]; ++A){
                            int thread = 0;
                            #ifdef _OPENMP
                                 thread = omp_get_thread_num();
                            #endif
                            double** CBDp = CBD[thread]->pointer();

                            // g(A'C|BD) = b(A'C|Q) b(Q|BD)
                            C_DGEMM('T', 'N', navirpi_[hC], navirpi_[hB]*navirpi_[hD], nQ_, 1.0, bQvvAp[0]+block_AB[hAC][hA].first+A*navirpi_[hC], bQabA_mo_->coldim(hAC), bQvvAp[0]+block_AB[hBD][hB].first, bQabA_mo_->coldim(hBD), 0.0, CBDp[0], navirpi_[hB]*navirpi_[hD]);

                            // g(A'C|BD) -> g(A'C|DB)
                            for (int B = 0; B < navirpi_[hB]; ++B){
                                for (int D = 0; D < navirpi_[hD]; ++D)
                                    CDB[thread]->set_column(0, D*navirpi_[hB]+B, CBD[thread]->get_column(0, B*navirpi_[hD]+D));
                            }

                            double** CDBp = CDB[thread]->pointer();

                            // G<IJ|A'B> = lambda<IJ|CD> g(A'C|DB)
                            C_DGEMM('N', 'N', Gaa.params->rowtot[hIJ], navirpi_[hB], navirpi_[hC]*navirpi_[hD], 1.0, Laa.matrix[hIJ][0]+block_AB[hCD][hC].first, Laa.params->coltot[hIJ], CDBp[0], navirpi_[hB], 1.0, Gaa.matrix[hIJ][0]+block_AB[hAB][hA].first+A*navirpi_[hB], Gaa.params->coltot[hIJ]);
                        }

                    }
                    global_dpd_->buf4_mat_irrep_wrt(&Gaa, hIJ);
                    global_dpd_->buf4_mat_irrep_close(&Gaa, hIJ);

                    global_dpd_->buf4_mat_irrep_close(&Laa, hIJ);
                }
            }
        }
    }

    global_dpd_->buf4_close(&Laa);
    global_dpd_->buf4_close(&Gaa);


    /********** Beta-Beta **********/

    // Put detailed information of b(Q|ab) block into 'block_ab'
    std::vector<std::vector<std::pair<long int, long int>>> block_ab;
    for (int hab = 0; hab < nirrep_; ++hab){
        long int entrance = 0;
        std::vector<std::pair<long int, long int>> subblock;
        for (int ha = 0; ha < nirrep_; ++ha){
            int hb = hab ^ ha;
            std::pair<long int, long int> subsubblock(entrance, nbvirpi_[ha] * nbvirpi_[hb]);
            subblock.push_back(subsubblock);
            entrance += subsubblock.second;
        }
        block_ab.push_back(subblock);
    }

    /*
     * Intermediate G <ij|ab> = 1/2 lambda<ij|cd> gbar<ab|cd>
     *                        = lambda<ij|cd> g(ac|bd)
     */

    dpdbuf4 Lbb, Gbb;

    global_dpd_->buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    global_dpd_->buf4_init(&Gbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                           ID("[o,o]"), ID("[v,v]"), 0, "tau(temp) <oo|vv>");
    global_dpd_->buf4_scm(&Gbb, 0.0);

    for (int hac = 0; hac < nirrep_; ++hac){
        for (int ha = 0; ha < nirrep_; ++ha){
            int hc = hac ^ ha;
            int hbd = hac;
            for (int hb = 0; hb < nirrep_; ++hb){
                int hd = hbd ^ hb;
                int hab = ha ^ hb;
                int hcd = hc ^ hd;
                int hij = hcd;

                if (Lbb.params->rowtot[hij] > 0 && Lbb.params->coltot[hcd] > 0 && Gbb.params->rowtot[hij] > 0 && Gbb.params->coltot[hab] > 0
                        && nbvirpi_[ha] > 0 && nbvirpi_[hc] > 0 && nbvirpi_[hb] > 0 && nbvirpi_[hd] > 0){
                    double** bQvvBp = bQabB_mo_->pointer(hac);

                    global_dpd_->buf4_mat_irrep_init(&Lbb, hij);
                    global_dpd_->buf4_mat_irrep_rd(&Lbb, hij);
                    global_dpd_->buf4_mat_irrep_init(&Gbb, hij);
                    global_dpd_->buf4_mat_irrep_rd(&Gbb, hij);

                    if (hb == hd){
                        std::vector<SharedMatrix> cbd;
                        for (int i = 0; i < nthreads; ++i){
                            cbd.push_back(SharedMatrix(new Matrix("g(a'c|bd)", nbvirpi_[hc], nbvirpi_[hb]*nbvirpi_[hd])));
                        }
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int a = 0; a < nbvirpi_[ha]; ++a){
                            int thread = 0;
                            #ifdef _OPENMP
                                 thread = omp_get_thread_num();
                            #endif
                            double** cbdp = cbd[thread]->pointer();
                            // g(a'c|bd) = b(a'c|Q) b(Q|bd)
                            C_DGEMM('T', 'N', nbvirpi_[hc], nbvirpi_[hb]*nbvirpi_[hd], nQ_, 1.0, bQvvBp[0]+block_ab[hac][ha].first+a*nbvirpi_[hc], bQabB_mo_->coldim(hac), bQvvBp[0]+block_ab[hbd][hb].first, bQabB_mo_->coldim(hbd), 0.0, cbdp[0], nbvirpi_[hb]*nbvirpi_[hd]);
                            // G<ij|a'b> = lambda<ij|cd> g(a'c|db)
                            C_DGEMM('N', 'N', Gbb.params->rowtot[hij], nbvirpi_[hb], nbvirpi_[hc]*nbvirpi_[hd], 1.0, Lbb.matrix[hij][0]+block_ab[hcd][hc].first, Lbb.params->coltot[hij], cbdp[0], nbvirpi_[hb], 1.0, Gbb.matrix[hij][0]+block_ab[hab][ha].first+a*nbvirpi_[hb], Gbb.params->coltot[hij]);
                        }
                    }
                    else{
                        std::vector<SharedMatrix> cbd, cdb;
                        for (int i = 0; i < nthreads; ++i){
                            cbd.push_back(SharedMatrix(new Matrix("g(a'c|bd)", nbvirpi_[hc], nbvirpi_[hb]*nbvirpi_[hd])));
                            cdb.push_back(SharedMatrix(new Matrix("g(a'c|db)", nbvirpi_[hc], nbvirpi_[hd]*nbvirpi_[hb])));
                        }
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int a = 0; a < nbvirpi_[ha]; ++a){
                            int thread = 0;
                            #ifdef _OPENMP
                                 thread = omp_get_thread_num();
                            #endif
                            double** cbdp = cbd[thread]->pointer();

                            // g(a'c|bd) = b(a'c|Q) b(Q|bd)
                            C_DGEMM('T', 'N', nbvirpi_[hc], nbvirpi_[hb]*nbvirpi_[hd], nQ_, 1.0, bQvvBp[0]+block_ab[hac][ha].first+a*nbvirpi_[hc], bQabB_mo_->coldim(hac), bQvvBp[0]+block_ab[hbd][hb].first, bQabB_mo_->coldim(hbd), 0.0, cbdp[0], nbvirpi_[hb]*nbvirpi_[hd]);

                            // g(a'c|bd) -> g(a'c|db)
                            for (int b = 0; b < nbvirpi_[hb]; ++b){
                                for (int d = 0; d < nbvirpi_[hd]; ++d)
                                    cdb[thread]->set_column(0, d*nbvirpi_[hb]+b, cbd[thread]->get_column(0, b*nbvirpi_[hd]+d));
                            }

                            double** cdbp = cdb[thread]->pointer();

                            // G<ij|a'b> = lambda<ij|cd> g(a'c|db)
                            C_DGEMM('N', 'N', Gbb.params->rowtot[hij], nbvirpi_[hb], nbvirpi_[hc]*nbvirpi_[hd], 1.0, Lbb.matrix[hij][0]+block_ab[hcd][hc].first, Lbb.params->coltot[hij], cdbp[0], nbvirpi_[hb], 1.0, Gbb.matrix[hij][0]+block_ab[hab][ha].first+a*nbvirpi_[hb], Gbb.params->coltot[hij]);
                        }

                    }
                    global_dpd_->buf4_mat_irrep_wrt(&Gbb, hij);
                    global_dpd_->buf4_mat_irrep_close(&Gbb, hij);

                    global_dpd_->buf4_mat_irrep_close(&Lbb, hij);
                }
            }
        }
    }

    global_dpd_->buf4_close(&Lbb);
    global_dpd_->buf4_close(&Gbb);

    /********** Alpha-Beta **********/

    // Put detailed information of Ab block (as in lambda<Ij|Ab>) into 'block_Ab'
    std::vector<std::vector<std::pair<long int, long int>>> block_Ab;
    for (int hAb = 0; hAb < nirrep_; ++hAb){
        long int entrance = 0;
        std::vector<std::pair<long int, long int>> subblock;
        for (int hA = 0; hA < nirrep_; ++hA){
            int hb = hAb ^ hA;
            std::pair<long int, long int> subsubblock(entrance, navirpi_[hA] * nbvirpi_[hb]);
            subblock.push_back(subsubblock);
            entrance += subsubblock.second;
        }
        block_Ab.push_back(subblock);
    }

    /*
     * Intermediate G<Ij|Ab> = lambda<Ij|Cd> gbar<Ab|Cd>
     *                       = lambda<Ij|Cd> g(AC|bd)
     */

    dpdbuf4 Lab, Gab;

    global_dpd_->buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    global_dpd_->buf4_init(&Gab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                           ID("[O,o]"), ID("[V,v]"), 0, "tau(temp) <Oo|Vv>");
    global_dpd_->buf4_scm(&Gab, 0.0);

    for (int hAC = 0; hAC < nirrep_; ++hAC){
        for (int hA = 0; hA < nirrep_; ++hA){
            int hC = hAC ^ hA;
            int hbd = hAC;
            for (int hb = 0; hb < nirrep_; ++hb){
                int hd = hbd ^ hb;
                int hAb = hA ^ hb;
                int hCd = hC ^ hd;
                int hIj = hCd;

                if (Lab.params->rowtot[hIj] > 0 && Lab.params->coltot[hCd] > 0 && Gab.params->rowtot[hIj] > 0 && Gab.params->coltot[hAb] > 0
                        && navirpi_[hA] > 0 && navirpi_[hC] > 0 && nbvirpi_[hb] > 0 && nbvirpi_[hd] > 0){
                    double** bQvvAp = bQabA_mo_->pointer(hAC);
                    double** bQvvBp = bQabB_mo_->pointer(hbd);

                    global_dpd_->buf4_mat_irrep_init(&Lab, hIj);
                    global_dpd_->buf4_mat_irrep_rd(&Lab, hIj);
                    global_dpd_->buf4_mat_irrep_init(&Gab, hIj);
                    global_dpd_->buf4_mat_irrep_rd(&Gab, hIj);

                    if (hb == hd){
                        std::vector<SharedMatrix> Cbd;
                        for (int i = 0; i < nthreads; ++i){
                            Cbd.push_back(SharedMatrix(new Matrix("g(A'C|bd)", navirpi_[hC], nbvirpi_[hb]*nbvirpi_[hd])));
                        }
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int A = 0; A < navirpi_[hA]; ++A){
                            int thread = 0;
                            #ifdef _OPENMP
                                 thread = omp_get_thread_num();
                            #endif
                            double** Cbdp = Cbd[thread]->pointer();
                            // g(A'C|bd) = b(A'C|Q) b(Q|bd)
                            C_DGEMM('T', 'N', navirpi_[hC], nbvirpi_[hb]*nbvirpi_[hd], nQ_, 1.0, bQvvAp[0]+block_AB[hAC][hA].first+A*navirpi_[hC], bQabA_mo_->coldim(hAC), bQvvBp[0]+block_ab[hbd][hb].first, bQabB_mo_->coldim(hbd), 0.0, Cbdp[0], nbvirpi_[hb]*nbvirpi_[hd]);
                            // G<Ij|A'b> = lambda<Ij|Cd> g(A'C|db)
                            C_DGEMM('N', 'N', Gab.params->rowtot[hIj], nbvirpi_[hb], navirpi_[hC]*nbvirpi_[hd], 1.0, Lab.matrix[hIj][0]+block_Ab[hCd][hC].first, Lab.params->coltot[hIj], Cbdp[0], nbvirpi_[hb], 1.0, Gab.matrix[hIj][0]+block_Ab[hAb][hA].first+A*nbvirpi_[hb], Gab.params->coltot[hIj]);
                        }
                    }
                    else{
                        std::vector<SharedMatrix> Cbd, Cdb;
                        for (int i = 0; i < nthreads; ++i){
                            Cbd.push_back(SharedMatrix(new Matrix("g(A'C|bd)", navirpi_[hC], nbvirpi_[hb]*nbvirpi_[hd])));
                            Cdb.push_back(SharedMatrix(new Matrix("g(A'C|db)", navirpi_[hC], nbvirpi_[hd]*nbvirpi_[hb])));
                        }
                        #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                        for (int A = 0; A < navirpi_[hA]; ++A){
                            int thread = 0;
                            #ifdef _OPENMP
                                 thread = omp_get_thread_num();
                            #endif
                            double** Cbdp = Cbd[thread]->pointer();

                            // g(A'C|bd) = b(A'C|Q) b(Q|bd)
                            C_DGEMM('T', 'N', navirpi_[hC], nbvirpi_[hb]*nbvirpi_[hd], nQ_, 1.0, bQvvAp[0]+block_AB[hAC][hA].first+A*navirpi_[hC], bQabA_mo_->coldim(hAC), bQvvBp[0]+block_ab[hbd][hb].first, bQabB_mo_->coldim(hbd), 0.0, Cbdp[0], nbvirpi_[hb]*nbvirpi_[hd]);

                            // g(A'C|bd) -> g(A'C|db)
                            for (int b = 0; b < nbvirpi_[hb]; ++b){
                                for (int d = 0; d < nbvirpi_[hd]; ++d)
                                    Cdb[thread]->set_column(0, d*nbvirpi_[hb]+b, Cbd[thread]->get_column(0, b*nbvirpi_[hd]+d));
                            }

                            double** Cdbp = Cdb[thread]->pointer();

                            // G<Ij|A'b> = lambda<Ij|Cd> g(A'C|db)
                            C_DGEMM('N', 'N', Gab.params->rowtot[hIj], nbvirpi_[hb], navirpi_[hC]*nbvirpi_[hd], 1.0, Lab.matrix[hIj][0]+block_Ab[hCd][hC].first, Lab.params->coltot[hIj], Cdbp[0], nbvirpi_[hb], 1.0, Gab.matrix[hIj][0]+block_Ab[hAb][hA].first+A*nbvirpi_[hb], Gab.params->coltot[hIj]);
                        }

                    }
                    global_dpd_->buf4_mat_irrep_wrt(&Gab, hIj);
                    global_dpd_->buf4_mat_irrep_close(&Gab, hIj);

                    global_dpd_->buf4_mat_irrep_close(&Lab, hIj);
                }
            }
        }
    }

    global_dpd_->buf4_close(&Lab);
    global_dpd_->buf4_close(&Gab);

    dcft_timer_off("DCFTSolver::DF lambda<ij|cd> gbar<ab|cd> (v3 in memory)");
}

/**
 * Form MO-based contraction [Gbar*Gamma]<q|p>
 * [Gbar*Gamma]<q|p> = Sum_rs Gbar<qs|pr> Gamma<r|s>
 */
void DCFTSolver::build_gbarGamma_UHF()
{
    dcft_timer_on("DCFTSolver::Gbar<QS|PR> Gamma<R|S> (FastBuilder)");

    build_gbarKappa_UHF();

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Form gamma<R|S> = kappa<R|S> + tau<R|S>
    mo_gammaA_ = SharedMatrix (new Matrix("MO-basis Gamma Alpha", nirrep_, nmopi_, nmopi_));
//    mo_gammaA_->copy(kappa_mo_a_);
//    mo_gammaA_->add(mo_tauA_);
    mo_gbarGamma_A_ = SharedMatrix (new Matrix("MO-basis Gbar_Gamma_A", nirrep_, nmopi_, nmopi_))    ;
    mo_gammaB_ = SharedMatrix (new Matrix("MO-basis Gamma Beta", nirrep_, nmopi_, nmopi_));
//    mo_gammaB_->copy(kappa_mo_b_);
//    mo_gammaB_->add(mo_tauB_);
    mo_gbarGamma_B_ = SharedMatrix (new Matrix("MO-basis Gbar_Gamma_B", nirrep_, nmopi_, nmopi_))    ;

    mo_gammaA_->copy(mo_tauA_);
    mo_gammaB_->copy(mo_tauB_);

    // Put detailed information of b(Q|pq) block into 'block'
    std::vector<std::vector<std::pair<long int, long int>>> block;
    for (int hpq = 0; hpq < nirrep_; ++hpq){
        long int entrance = 0;
        std::vector<std::pair<long int, long int>> subblock;
        for (int hp = 0; hp < nirrep_; ++hp){
            int hq = hpq ^ hp;
            std::pair<long int, long int> subsubblock(entrance, nsopi_[hp] * nsopi_[hq]);
            subblock.push_back(subsubblock);
            entrance += subsubblock.second;
        }
        block.push_back(subblock);
    }

    /*
     *  f_tilde <Q|P> = gbar<QS|PR> gamma<R|S> + gbar<Qs|Pr> gamma<r|s>
     *             = g(QP|SR) gamma<R|S> - g(QR|SP) gamma<R|S> + g(QP|sr) gamma<r|s>
     *
     *  f_tilde <q|p> = gbar<qs|pr> gamma<r|s> + gbar<qS|pR> gamma<R|S>
     *             = g(qp|sr) gamma<r|s> - g(qr|sp) gamma<r|s> + g(qp|SR) gamma<R|S>
     */

    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int hQ = 0; hQ < nirrep_; ++hQ){
        int hP = hQ;
        if (nsopi_[hQ] > 0){
            double** tFAp = mo_gbarGamma_A_->pointer(hQ);
            double** tFBp = mo_gbarGamma_B_->pointer(hQ);

            double** bQpqAp = bQpqA_mo_->pointer(0);
            double** bQpqBp = bQpqB_mo_->pointer(0);

            // (Q) = b(Q|SR)*gamma<R|S> + b(Q|sr)*gamma<r|s>
            SharedMatrix Q (new Matrix("b(Q|SR)gamma<R|S>", 1, nQ_));
            double** Qp = Q->pointer();
            for (int hR = 0; hR < nirrep_; ++hR){
                int hS = hR;
                if (nsopi_[hR] > 0){
                    double** gamma_rsAp = mo_gammaA_->pointer(hR);
                    double** gamma_rsBp = mo_gammaB_->pointer(hR);
                    // (Q) = b(Q|SR) gamma<R|S>
                    C_DGEMV('N', nQ_, nsopi_[hR]*nsopi_[hS], 1.0, bQpqAp[0]+block[0][hR].first, bQpqA_mo_->coldim(0), gamma_rsAp[0], 1, 1.0, Qp[0], 1);
                    // (Q) += b(Q|sr) gamma<r|s>
                    C_DGEMV('N', nQ_, nsopi_[hR]*nsopi_[hS], 1.0, bQpqBp[0]+block[0][hR].first, bQpqB_mo_->coldim(0), gamma_rsBp[0], 1, 1.0, Qp[0], 1);
                }
            }

            // f_tilde <Q|P> = b(QP|Aux)*(Aux) where (Aux) = (Q)
            C_DGEMV('T', nQ_, nsopi_[hP]*nsopi_[hQ], 1.0, bQpqAp[0]+block[0][hP].first, bQpqA_mo_->coldim(0), Qp[0], 1, 0.0, tFAp[0], 1);

            // f_tilde <q|p> = b(qp|Aux)*(Aux) where (Aux) = (Q)
            C_DGEMV('T', nQ_, nsopi_[hP]*nsopi_[hQ], 1.0, bQpqBp[0]+block[0][hP].first, bQpqB_mo_->coldim(0), Qp[0], 1, 0.0, tFBp[0], 1);
        }

    }


    // f_tilde <Q|P> -= b(QR|Aux) b(Aux|SP) gamma<R|S>
    for (int hQ = 0; hQ < nirrep_; ++hQ){
        int hP = hQ;
        if (nsopi_[hQ] > 0){
            for (int hR = 0; hR < nirrep_; ++hR){
                int hS = hR;
                if (nsopi_[hR] > 0){
                    double** bQpqAp = bQpqA_mo_->pointer(hQ^hR);
                    double** gamma_rsA_p = mo_gammaA_->pointer(hR);

                    std::vector<SharedMatrix> RS;
                    for (int i = 0; i < nthreads; ++i){
                        RS.push_back(SharedMatrix(new Matrix("<Q'P'|RS>", nsopi_[hR], nsopi_[hS])));
                    }

                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int Q = 0; Q < nsopi_[hQ]; ++Q){
                        for (int P = Q; P < nsopi_[hP]; ++P){
                            int thread = 0;
                            #ifdef _OPENMP
                                thread = omp_get_thread_num();
                            #endif
                            double** RSp = RS[thread]->pointer();

                            // <Q'P'|RS> = b(Q'R|Aux) b(Aux|P'S)
                            C_DGEMM('T', 'N', nsopi_[hR], nsopi_[hS], nQ_, 1.0, bQpqAp[0]+block[hQ^hR][hQ].first+Q*nsopi_[hR], bQpqA_mo_->coldim(hQ^hR), bQpqAp[0]+block[hP^hS][hP].first+P*nsopi_[hS], bQpqA_mo_->coldim(hP^hS), 0.0, RSp[0], nsopi_[hS]);
                            // - <Q'P'|RS> * gamma<R|S>
                            double value = - C_DDOT(nsopi_[hR]*nsopi_[hS], RSp[0], 1, gamma_rsA_p[0], 1);
                            mo_gbarGamma_A_->add(hP, Q, P, value);
                            if (Q != P){
                                mo_gbarGamma_A_->add(hP, P, Q, value);
                            }
                        }
                    }
                }
            }
        }
    }


    // f_tilde <q|p> -= b(qr|Aux) b(Aux|sp) gamma<r|s>
    for (int hq = 0; hq < nirrep_; ++hq){
        int hp = hq;
        if (nsopi_[hq] > 0){
            for (int hr = 0; hr < nirrep_; ++hr){
                int hs = hr;
                if (nsopi_[hr] > 0){
                    double** bQpqBp = bQpqB_mo_->pointer(hq^hr);
                    double** gamma_rsB_p = mo_gammaB_->pointer(hr);

                    std::vector<SharedMatrix> rs;
                    for (int i = 0; i < nthreads; ++i){
                        rs.push_back(SharedMatrix (new Matrix("<q'p'|rs>", nsopi_[hr], nsopi_[hs])));
                    }

                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int q = 0; q < nsopi_[hq]; ++q){
                        for (int p = q; p < nsopi_[hp]; ++p){
                            int thread = 0;
                            #ifdef _OPENMP
                                thread = omp_get_thread_num();
                            #endif
                            double** rsp = rs[thread]->pointer();

                            // <q'p'|rs> = b(q'r|Aux) b(Aux|p's)
                            C_DGEMM('T', 'N', nsopi_[hr], nsopi_[hs], nQ_, 1.0, bQpqBp[0]+block[hq^hr][hq].first+q*nsopi_[hr], bQpqB_mo_->coldim(hq^hr), bQpqBp[0]+block[hp^hs][hp].first+p*nsopi_[hs], bQpqB_mo_->coldim(hp^hs), 0.0, rsp[0], nsopi_[hs]);
                            // - <q'p'|rs> * gamma<r|s>
                            double value = - C_DDOT(nsopi_[hr]*nsopi_[hs], rsp[0], 1, gamma_rsB_p[0], 1);
                            mo_gbarGamma_B_->add(hp, q, p, value);
                            if (q != p){
                                mo_gbarGamma_B_->add(hp, p, q, value);
                            }
                        }
                    }
                }
            }
        }
    }

    mo_gbarGamma_A_->add(mo_gbarKappa_A_);
    mo_gbarGamma_B_->add(mo_gbarKappa_B_);

    dcft_timer_off("DCFTSolver::Gbar<QS|PR> Gamma<R|S> (FastBuilder)");
}

void DCFTSolver::build_gbarKappa_UHF()
{
    dcft_timer_on("DCFTSolver::Gbar<QS|PR> Kappa<R|S>");

    formb_pq_scf();

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    mo_gbarKappa_A_ = SharedMatrix (new Matrix("MO-basis Gbar_Kappa_A", nirrep_, nmopi_, nmopi_))    ;
    mo_gbarKappa_B_ = SharedMatrix (new Matrix("MO-basis Gbar_Kappa_B", nirrep_, nmopi_, nmopi_))    ;

    // Put detailed information of b(Q|pq) block into 'block'
    std::vector<std::vector<std::pair<long int, long int>>> block;
    for (int hpq = 0; hpq < nirrep_; ++hpq){
        long int entrance = 0;
        std::vector<std::pair<long int, long int>> subblock;
        for (int hp = 0; hp < nirrep_; ++hp){
            int hq = hpq ^ hp;
            std::pair<long int, long int> subsubblock(entrance, nsopi_[hp] * nsopi_[hq]);
            subblock.push_back(subsubblock);
            entrance += subsubblock.second;
        }
        block.push_back(subblock);
    }

    /*
     *  f_tilde <Q|P> = gbar<QS|PR> gamma<R|S> + gbar<Qs|Pr> gamma<r|s>
     *             = g(QP|SR) gamma<R|S> - g(QR|SP) gamma<R|S> + g(QP|sr) gamma<r|s>
     *
     *  f_tilde <q|p> = gbar<qs|pr> gamma<r|s> + gbar<qS|pR> gamma<R|S>
     *             = g(qp|sr) gamma<r|s> - g(qr|sp) gamma<r|s> + g(qp|SR) gamma<R|S>
     */

    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (int hQ = 0; hQ < nirrep_; ++hQ){
        int hP = hQ;
        if (nsopi_[hQ] > 0){
            double** tFAp = mo_gbarKappa_A_->pointer(hQ);
            double** tFBp = mo_gbarKappa_B_->pointer(hQ);

            double** bQpqAp = bQpqA_mo_scf_->pointer(0);
            double** bQpqBp = bQpqB_mo_scf_->pointer(0);

            // (Q) = b(Q|SR)*gamma<R|S> + b(Q|sr)*gamma<r|s>
            SharedMatrix Q (new Matrix("b(Q|SR)gamma<R|S>", 1, nQ_scf_));
            double** Qp = Q->pointer();
            for (int hR = 0; hR < nirrep_; ++hR){
                int hS = hR;
                if (nsopi_[hR] > 0){
                    double** gamma_rsAp = kappa_mo_a_->pointer(hR);
                    double** gamma_rsBp = kappa_mo_b_->pointer(hR);
                    // (Q) = b(Q|SR) gamma<R|S>
                    C_DGEMV('N', nQ_scf_, nsopi_[hR]*nsopi_[hS], 1.0, bQpqAp[0]+block[0][hR].first, bQpqA_mo_scf_->coldim(0), gamma_rsAp[0], 1, 1.0, Qp[0], 1);
                    // (Q) += b(Q|sr) gamma<r|s>
                    C_DGEMV('N', nQ_scf_, nsopi_[hR]*nsopi_[hS], 1.0, bQpqBp[0]+block[0][hR].first, bQpqB_mo_scf_->coldim(0), gamma_rsBp[0], 1, 1.0, Qp[0], 1);
                }
            }

            // f_tilde <Q|P> = b(QP|Aux)*(Aux) where (Aux) = (Q)
            C_DGEMV('T', nQ_scf_, nsopi_[hP]*nsopi_[hQ], 1.0, bQpqAp[0]+block[0][hP].first, bQpqA_mo_scf_->coldim(0), Qp[0], 1, 0.0, tFAp[0], 1);

            // f_tilde <q|p> = b(qp|Aux)*(Aux) where (Aux) = (Q)
            C_DGEMV('T', nQ_scf_, nsopi_[hP]*nsopi_[hQ], 1.0, bQpqBp[0]+block[0][hP].first, bQpqB_mo_scf_->coldim(0), Qp[0], 1, 0.0, tFBp[0], 1);
        }

    }


    // f_tilde <Q|P> -= b(QR|Aux) b(Aux|SP) gamma<R|S>
    for (int hQ = 0; hQ < nirrep_; ++hQ){
        int hP = hQ;
        if (nsopi_[hQ] > 0){
            for (int hR = 0; hR < nirrep_; ++hR){
                int hS = hR;
                if (nsopi_[hR] > 0){
                    double** bQpqAp = bQpqA_mo_scf_->pointer(hQ^hR);
                    double** gamma_rsA_p = kappa_mo_a_->pointer(hR);

                    std::vector<SharedMatrix> RS;
                    for (int i = 0; i < nthreads; ++i){
                        RS.push_back(SharedMatrix(new Matrix("<Q'P'|RS>", nsopi_[hR], nsopi_[hS])));
                    }

                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int Q = 0; Q < nsopi_[hQ]; ++Q){
                        for (int P = Q; P < nsopi_[hP]; ++P){
                            int thread = 0;
                            #ifdef _OPENMP
                                thread = omp_get_thread_num();
                            #endif
                            double** RSp = RS[thread]->pointer();

                            // <Q'P'|RS> = b(Q'R|Aux) b(Aux|P'S)
                            C_DGEMM('T', 'N', nsopi_[hR], nsopi_[hS], nQ_scf_, 1.0, bQpqAp[0]+block[hQ^hR][hQ].first+Q*nsopi_[hR], bQpqA_mo_scf_->coldim(hQ^hR), bQpqAp[0]+block[hP^hS][hP].first+P*nsopi_[hS], bQpqA_mo_scf_->coldim(hP^hS), 0.0, RSp[0], nsopi_[hS]);
                            // - <Q'P'|RS> * gamma<R|S>
                            double value = - C_DDOT(nsopi_[hR]*nsopi_[hS], RSp[0], 1, gamma_rsA_p[0], 1);
                            mo_gbarKappa_A_->add(hP, Q, P, value);
                            if (Q != P){
                                mo_gbarKappa_A_->add(hP, P, Q, value);
                            }
                        }
                    }
                }
            }
        }
    }


    // f_tilde <q|p> -= b(qr|Aux) b(Aux|sp) gamma<r|s>
    for (int hq = 0; hq < nirrep_; ++hq){
        int hp = hq;
        if (nsopi_[hq] > 0){
            for (int hr = 0; hr < nirrep_; ++hr){
                int hs = hr;
                if (nsopi_[hr] > 0){
                    double** bQpqBp = bQpqB_mo_scf_->pointer(hq^hr);
                    double** gamma_rsB_p = kappa_mo_b_->pointer(hr);

                    std::vector<SharedMatrix> rs;
                    for (int i = 0; i < nthreads; ++i){
                        rs.push_back(SharedMatrix (new Matrix("<q'p'|rs>", nsopi_[hr], nsopi_[hs])));
                    }

                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int q = 0; q < nsopi_[hq]; ++q){
                        for (int p = q; p < nsopi_[hp]; ++p){
                            int thread = 0;
                            #ifdef _OPENMP
                                thread = omp_get_thread_num();
                            #endif
                            double** rsp = rs[thread]->pointer();

                            // <q'p'|rs> = b(q'r|Aux) b(Aux|p's)
                            C_DGEMM('T', 'N', nsopi_[hr], nsopi_[hs], nQ_scf_, 1.0, bQpqBp[0]+block[hq^hr][hq].first+q*nsopi_[hr], bQpqB_mo_scf_->coldim(hq^hr), bQpqBp[0]+block[hp^hs][hp].first+p*nsopi_[hs], bQpqB_mo_scf_->coldim(hp^hs), 0.0, rsp[0], nsopi_[hs]);
                            // - <q'p'|rs> * gamma<r|s>
                            double value = - C_DDOT(nsopi_[hr]*nsopi_[hs], rsp[0], 1, gamma_rsB_p[0], 1);
                            mo_gbarKappa_B_->add(hp, q, p, value);
                            if (q != p){
                                mo_gbarKappa_B_->add(hp, p, q, value);
                            }
                        }
                    }
                }
            }
        }
    }

    bQpqA_mo_scf_.reset();
    bQpqB_mo_scf_.reset();

    dcft_timer_off("DCFTSolver::Gbar<QS|PR> Kappa<R|S>");

}

/**
  * Form J(P,Q)^-1/2 for SCF terms
  */
void DCFTSolver::formJm12_scf(std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<BasisSet> zero)
{
//    outfile->Printf("\tForming J(P,Q)^-1/2 ...\n\n");
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    double **J = block_matrix(nQ_scf_, nQ_scf_);
    Jm12_scf_ = block_matrix(nQ_scf_, nQ_scf_);

    // => Integrals <= //
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary, zero, auxiliary, zero));
    std::vector<std::shared_ptr<TwoBodyAOInt> > Jint;
    std::vector<const double*> buffer;
    for (int t = 0; t < nthreads; t++){
        Jint.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri()));
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
    int lwork = nQ_scf_ * 3;
    double* eigval = init_array(nQ_scf_);
    double* work = init_array(lwork);
    int status = C_DSYEV('v', 'u', nQ_scf_, J[0], nQ_scf_, eigval, work, lwork);
    if(status){
        throw PsiException("Diagonalization of J failed", __FILE__, __LINE__);
    }
    free(work);

    // Now J contains eigenvectors of the original J
    double **J_copy = block_matrix(nQ_scf_, nQ_scf_);
    C_DCOPY(nQ_scf_ * nQ_scf_, J[0], 1, J_copy[0], 1);

    // Now form J^-1/2 = U(T) * j^-1/2 * U
    // where j^-1/2 is the diagonal matrix of inverse square
    // of the eigenvalues, and U is the matrix of eigenvectors of J
    for (int i = 0; i < nQ_scf_; ++i){
        eigval[i] = (eigval[i] < 1.0E-10) ? 0.0 : 1.0 / sqrt(eigval[i]);
        // scale one set of eigenvectors by the diagonal elements j^-1/2
        C_DSCAL(nQ_scf_, eigval[i], J[i], 1);
    }
    free(eigval);

    // J^-1/2 = J_copy(T) * J
    C_DGEMM('t', 'n', nQ_scf_, nQ_scf_, nQ_scf_, 1.0, J_copy[0], nQ_scf_, J[0], nQ_scf_, 0.0, Jm12_scf_[0], nQ_scf_);
    free_block(J);
    free_block(J_copy);

}

/**
  * Form b(Q|mn) for SCF terms
  */
void DCFTSolver::formb_ao_scf(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<BasisSet> zero)
{
    bQmn_ao_scf_ = SharedMatrix(new Matrix(nQ_scf_, nso_ * nso_));
    double **Ap = bQmn_ao_scf_->pointer();
    double **Bp = block_matrix(nQ_scf_, nso_ * nso_);

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    std::shared_ptr<ERISieve> sieve = std::shared_ptr<ERISieve>(new ERISieve(primary, 1.0E-20));
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
    std::shared_ptr<IntegralFactory> rifactory2(new IntegralFactory(auxiliary, zero, primary, primary));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    std::vector<const double*> buffer;
    for (int t = 0; t < nthreads; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory2->eri()));
        buffer.push_back(eri[t]->buffer());
    }

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {

        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop  = Pstarts[block+1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary->shell(Pstart).function_index();
        int pstop  = (Pstop == auxiliary->nshell() ? nQ_scf_ : auxiliary->shell(Pstop ).function_index());
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

    C_DGEMM('N','N', nQ_scf_, nso_ * nso_, nQ_scf_, 1.0, Jm12_scf_[0], nQ_scf_, Bp[0], nso_ * nso_, 0.0, Ap[0], nso_ * nso_);
}

/**
  * Transform b(Q|mu,nu) from AO basis to SO basis for SCF terms
  */
void DCFTSolver::transform_b_ao2so_scf()
{
    dcft_timer_on("DCFTSolver::Transform b(Q|mn) AO-basis -> SO-basis");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    double** bQmn_ao_p = bQmn_ao_scf_->pointer();

    // Set up dimensions for SO-basis b(Q|mn)
    Dimension Q(nirrep_), mn(nirrep_);
    for (int hn = 0; hn < nirrep_; ++hn){
        Q[hn] = nQ_scf_;
        for (int hm = 0; hm < nirrep_; ++hm){
            mn[hm ^ hn] += nsopi_[hm] * nsopi_[hn];
        }
    }
    bQmn_so_scf_ = SharedMatrix(new Matrix("Fully-transformed b", Q, mn));

    std::vector<int> offset(nirrep_);
    for (int h = 0; h < nirrep_; ++h){
        offset.push_back(0);
    }

    // AO-basis b(Q|mn) -> SO-basis b(Q|mn)
    for (int h = 0; h < nirrep_; ++h){
        double** bQmn_so_p = bQmn_so_scf_->pointer(h);
        for (int hm = 0; hm < nirrep_; ++hm){
            int hn = h ^ hm;
            if (nsopi_[hm] > 0 && nsopi_[hn] > 0){

                SharedMatrix tmp (new Matrix("Half-transformed b", nQ_scf_, nso_ * nsopi_[hn]));
                double** tmpp = tmp->pointer();
                double** ao2so_n_p = reference_wavefunction()->aotoso()->pointer(hn);
                double** ao2so_m_p = reference_wavefunction()->aotoso()->pointer(hm);
                // First-half transformation
                C_DGEMM('N', 'N', nQ_scf_ * nso_, nsopi_[hn], nso_, 1.0, bQmn_ao_p[0], nso_, ao2so_n_p[0], nsopi_[hn], 0.0, tmpp[0], nsopi_[hn]);
                // Second-half transformation
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int Q = 0; Q < nQ_scf_; ++Q){
                    C_DGEMM('T', 'N', nsopi_[hm], nsopi_[hn], nso_, 1.0, ao2so_m_p[0], nsopi_[hm], tmpp[Q], nsopi_[hn], 0.0, bQmn_so_p[Q]+offset[h], nsopi_[hn]);
                }
            }
        offset[h] += nsopi_[hm] * nsopi_[hn];
        }
    }

    bQmn_ao_scf_.reset();

    dcft_timer_off("DCFTSolver::Transform b(Q|mn) AO-basis -> SO-basis");
}

/**
  * form b(Q,ij) for SCF terms
  */
void DCFTSolver::formb_oo_scf()
{
    dcft_timer_on("DCFTSolver::b(Q|mn) -> b(Q|ij)");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Set up dimensions for b(Q|IJ)
    Dimension OO(nirrep_), Q(nirrep_);
    for (int hI = 0; hI < nirrep_; ++hI){
        Q[hI] = nQ_scf_;
        for (int hJ = 0; hJ < nirrep_; ++hJ){
            OO[hI ^ hJ] += naoccpi_[hI] * naoccpi_[hJ];
        }
    }
    bQijA_mo_scf_ = SharedMatrix(new Matrix("b(Q|IJ)", Q, OO));

    std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
    for (int h = 0; h < nirrep_; ++h){
        offset_so.push_back(0);
        offset_mo.push_back(0);
    }

    for (int h = 0; h < nirrep_; ++h){
        double** bQmn_so_p = bQmn_so_scf_->pointer(h);
        double** bQijA_mo_p = bQijA_mo_scf_->pointer(h);
        for (int hI = 0; hI < nirrep_; ++hI){
            int hJ = h ^ hI;
            if (naoccpi_[hI] > 0 && naoccpi_[hJ] > 0){
                double** CaJp = Ca_->pointer(hJ);
                double** CaIp = Ca_->pointer(hI);
                SharedMatrix tmp (new Matrix("Half-transformed b_IJ", nQ_scf_, nsopi_[hI] * naoccpi_[hJ]));
                double** tmpp = tmp->pointer();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int Q = 0; Q < nQ_scf_; ++Q){
                    // First-half transformation
                    C_DGEMM('N', 'N', nsopi_[hI], naoccpi_[hJ], nsopi_[hJ], 1.0, bQmn_so_p[Q] + offset_so[h], nsopi_[hJ], CaJp[0], nsopi_[hJ], 0.0, tmpp[Q], naoccpi_[hJ]);
                    // Second-half transformation
                    C_DGEMM('T', 'N', naoccpi_[hI], naoccpi_[hJ], nsopi_[hI], 1.0, CaIp[0], nsopi_[hI], tmpp[Q], naoccpi_[hJ], 0.0, bQijA_mo_p[Q] + offset_mo[h], naoccpi_[hJ]);
                }
            }
            offset_so[h] += nsopi_[h ^ hI] * nsopi_[hI];
            offset_mo[h] += naoccpi_[h ^ hI] * naoccpi_[hI];
        }
    }

    if (options_.get_str("REFERENCE") != "RHF"){

        // Set up dimensions for b(Q|ij)
        Dimension oo(nirrep_), Q(nirrep_);
        for (int hi = 0; hi < nirrep_; ++hi){
            Q[hi] = nQ_scf_;
            for (int hj = 0; hj < nirrep_; ++hj){
                oo[hi ^ hj] += nboccpi_[hi] * nboccpi_[hj];
            }
        }

        bQijB_mo_scf_ = SharedMatrix(new Matrix("b(Q|ij)", Q, oo));

        std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
        for (int h = 0; h < nirrep_; ++h){
            offset_so.push_back(0);
            offset_mo.push_back(0);
        }

        for (int h = 0; h < nirrep_; ++h){
            double** bQmn_so_p = bQmn_so_scf_->pointer(h);
            double** bQijB_mo_p = bQijB_mo_scf_->pointer(h);
            for (int hi = 0; hi < nirrep_; ++hi){
                int hj = h ^ hi;
                if (nboccpi_[hi] > 0 && nboccpi_[hj] > 0){
                    double** Cbjp = Cb_->pointer(hj);
                    double** Cbip = Cb_->pointer(hi);
                    SharedMatrix tmp (new Matrix("Half-transformed b_ij", nQ_scf_, nsopi_[hi] * nboccpi_[hj]));
                    double** tmpp = tmp->pointer();
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int Q = 0; Q < nQ_scf_; ++Q){
                        // First-half transformation
                        C_DGEMM('N', 'N', nsopi_[hi], nboccpi_[hj], nsopi_[hj], 1.0, bQmn_so_p[Q] + offset_so[h], nsopi_[hj], Cbjp[0], nsopi_[hj], 0.0, tmpp[Q], nboccpi_[hj]);
                        // Second-half transformation
                        C_DGEMM('T', 'N', nboccpi_[hi], nboccpi_[hj], nsopi_[hi], 1.0, Cbip[0], nsopi_[hi], tmpp[Q], nboccpi_[hj], 0.0, bQijB_mo_p[Q] + offset_mo[h], nboccpi_[hj]);
                    }
                }
                offset_so[h] += nsopi_[h ^ hi] * nsopi_[hi];
                offset_mo[h] += nboccpi_[h ^ hi] * nboccpi_[hi];
            }
        }


    }

    dcft_timer_off("DCFTSolver::b(Q|mn) -> b(Q|ij)");
}

/**
  * form b(Q,pq) for SCF terms
  */
void DCFTSolver::formb_pq_scf()
{
    dcft_timer_on("DCFTSolver::b(Q|mn) -> b(Q|pq)");

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    // Set up dimensions for b(Aux|PQ)
    Dimension PQ(nirrep_), Aux(nirrep_);
    for (int hP = 0; hP < nirrep_; ++hP){
        Aux[hP] = nQ_scf_;
        for (int hQ = 0; hQ < nirrep_; ++hQ){
            PQ[hP ^ hQ] += nsopi_[hP] * nsopi_[hQ];
        }
    }
    bQpqA_mo_scf_ = SharedMatrix(new Matrix("b(Aux|PQ)", Aux, PQ));

    std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
    for (int h = 0; h < nirrep_; ++h){
        offset_so.push_back(0);
        offset_mo.push_back(0);
    }

    for (int h = 0; h < nirrep_; ++h){
        double** bQmn_so_p = bQmn_so_scf_->pointer(h);
        double** bQpqA_mo_p = bQpqA_mo_scf_->pointer(h);
        for (int hP = 0; hP < nirrep_; ++hP){
            int hQ = h ^ hP;
            if (nsopi_[hP] > 0 && nsopi_[hQ] > 0){
                double** Caqp = Ca_->pointer(hQ);
                double** Capp = Ca_->pointer(hP);
                SharedMatrix tmp (new Matrix("Half-transformed b_PQ", nQ_scf_, nsopi_[hP] * nsopi_[hQ]));
                double** tmpp = tmp->pointer();
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int Aux = 0; Aux < nQ_scf_; ++Aux){
                    // First-half transformation
                    C_DGEMM('N', 'N', nsopi_[hP], nsopi_[hQ], nsopi_[hQ], 1.0, bQmn_so_p[Aux] + offset_so[h], nsopi_[hQ], Caqp[0], nsopi_[hQ], 0.0, tmpp[Aux], nsopi_[hQ]);
                    // Second-half transformation
                    C_DGEMM('T', 'N', nsopi_[hP], nsopi_[hQ], nsopi_[hP], 1.0, Capp[0], nsopi_[hP], tmpp[Aux], nsopi_[hQ], 0.0, bQpqA_mo_p[Aux] + offset_mo[h], nsopi_[hQ]);
                }

            }
            offset_so[h] += nsopi_[h ^ hP] * nsopi_[hP];
            offset_mo[h] += nsopi_[h ^ hP] * nsopi_[hP];
        }
    }

    if (options_.get_str("REFERENCE") != "RHF"){
        // Set up dimensions for b(Aux|pq)
        bQpqB_mo_scf_ = SharedMatrix(new Matrix("b(Aux|pq)", Aux, PQ));

        std::vector<int> offset_so(nirrep_), offset_mo(nirrep_);
        for (int h = 0; h < nirrep_; ++h){
            offset_so.push_back(0);
            offset_mo.push_back(0);
        }

        for (int h = 0; h < nirrep_; ++h){
            double** bQmn_so_p = bQmn_so_scf_->pointer(h);
            double** bQpqB_mo_p = bQpqB_mo_scf_->pointer(h);
            for (int hp = 0; hp < nirrep_; ++hp){
                int hq = h ^ hp;
                if (nsopi_[hp] > 0 && nsopi_[hq] > 0){
                    double** Cbqp = Cb_->pointer(hq);
                    double** Cbpp = Cb_->pointer(hp);
                    SharedMatrix tmp (new Matrix("Half-transformed b_pq", nQ_scf_, nsopi_[hp] * nsopi_[hq]));
                    double** tmpp = tmp->pointer();
                    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                    for (int Aux = 0; Aux < nQ_scf_; ++Aux){
                        // First-half transformation
                        C_DGEMM('N', 'N', nsopi_[hp], nsopi_[hq], nsopi_[hq], 1.0, bQmn_so_p[Aux] + offset_so[h], nsopi_[hq], Cbqp[0], nsopi_[hq], 0.0, tmpp[Aux], nsopi_[hq]);
                        // Second-half transformation
                        C_DGEMM('T', 'N', nsopi_[hp], nsopi_[hq], nsopi_[hp], 1.0, Cbpp[0], nsopi_[hp], tmpp[Aux], nsopi_[hq], 0.0, bQpqB_mo_p[Aux] + offset_mo[h], nsopi_[hq]);
                    }
                }
                offset_so[h] += nsopi_[h ^ hp] * nsopi_[hp];
                offset_mo[h] += nsopi_[h ^ hp] * nsopi_[hp];
            }
        }
    }

    dcft_timer_off("DCFTSolver::b(Q|mn) -> b(Q|pq)");

}

}} // End namespace
