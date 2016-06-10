/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

/** Standard library includes */
#include <fstream>
#include <psifiles.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include <libmints/sieve.h>
#include "dfocc.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::tei_grad_corr()
{      
  //outfile->Printf("\ttei_grad_corr is starting... \n"); 

//===========================================================================================
//========================= Two-electron Gradient:Corr ======================================
//===========================================================================================
    // Two-electron gradients
    int df_ints_num_threads_ = 1;
    #ifdef _OPENMP
        df_ints_num_threads_ = omp_get_max_threads();
    #endif

    // Read in the basis set informations
    boost::shared_ptr<BasisSet> primary_ = BasisSet::pyconstruct_orbital(reference_wavefunction_->molecule(), 
        "BASIS", Process::environment.options.get_str("BASIS"));
    boost::shared_ptr<BasisSet> auxiliary_ = BasisSet::pyconstruct_auxiliary(reference_wavefunction_->molecule(), 
        "DF_BASIS_CC", Process::environment.options.get_str("DF_BASIS_CC"),
        "RIFIT", Process::environment.options.get_str("BASIS"));
    boost::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
    //auxiliary_->print();
    int nbasis = primary_->nbf();
    int nQ_corr = auxiliary_->nbf();
    

//===========================================================================================
//========================= Metric Gradient:Corr ============================================
//===========================================================================================
    // Read Gaux 
    Gaux = SharedTensor2d(new Tensor2d("2-Index Correlation TPDM (P|Q)", nQ_corr, nQ_corr));
    //Gaux->read(psio_, PSIF_DFOCC_DENS);
    Gaux->read_symm(psio_, PSIF_DFOCC_DENS);
 
    // JPQ_X
    timer_on("Grad: Metric:Corr");
    gradients["Metric:Corr"] = SharedMatrix(gradients["Nuclear"]->clone());
    gradients["Metric:Corr"]->set_name("Metric:Corr Gradient");
    gradients["Metric:Corr"]->zero();

    // => Sizing <= //
    int natom = primary_->molecule()->natom();
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();

    // => Integrals <= //
    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),auxiliary_,BasisSet::zero_ao_basis_set()));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > Jint;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        Jint.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));
    }

    // => Temporary Gradients <= //
    std::vector<SharedMatrix> Jtemps;
    for (int t = 0; t < df_ints_num_threads_; t++) {
         Jtemps.push_back(SharedMatrix(new Matrix("Jtemp", natom, 3)));
    }

    std::vector<std::pair<int,int> > PQ_pairs;
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int,int>(P,Q));
        }
    }

    int nthread_df = df_ints_num_threads_;
    #pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
    for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {

        int P = PQ_pairs[PQ].first;
        int Q = PQ_pairs[PQ].second;

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        Jint[thread]->compute_shell_deriv1(P,0,Q,0);
        const double* buffer = Jint[thread]->buffer();

        int nP = auxiliary_->shell(P).nfunction();
        int cP = auxiliary_->shell(P).ncartesian();
        int aP = auxiliary_->shell(P).ncenter();
        int oP = auxiliary_->shell(P).function_index();

        int nQ = auxiliary_->shell(Q).nfunction();
        int cQ = auxiliary_->shell(Q).ncartesian();
        int aQ = auxiliary_->shell(Q).ncenter();
        int oQ = auxiliary_->shell(Q).function_index();

        int ncart = cP * cQ;
        const double *Px = buffer + 0*ncart;
        const double *Py = buffer + 1*ncart;
        const double *Pz = buffer + 2*ncart;
        const double *Qx = buffer + 3*ncart;
        const double *Qy = buffer + 4*ncart;
        const double *Qz = buffer + 5*ncart;

        double perm = (P == Q ? 1.0 : 2.0);

        double** grad_Jp;
        grad_Jp = Jtemps[thread]->pointer();

        for (int p = 0; p < nP; p++) {
            for (int q = 0; q < nQ; q++) {

                    double Uval = perm * Gaux->get(p + oP, q + oQ);
                    grad_Jp[aP][0] -= Uval * (*Px);
                    grad_Jp[aP][1] -= Uval * (*Py);
                    grad_Jp[aP][2] -= Uval * (*Pz);
                    grad_Jp[aQ][0] -= Uval * (*Qx);
                    grad_Jp[aQ][1] -= Uval * (*Qy);
                    grad_Jp[aQ][2] -= Uval * (*Qz);

                Px++;
                Py++;
                Pz++;
                Qx++;
                Qy++;
                Qz++;
            }
        }
    }
    Gaux.reset();

    // => Temporary Gradient Reduction <= //
    for (int t = 0; t < df_ints_num_threads_; t++) {
         gradients["Metric:Corr"]->add(Jtemps[t]);
    }

    //gradients["Metric:Corr"]->print_atom_vector();
    timer_off("Grad: Metric:Corr");

//===========================================================================================
//========================= 3-Index Gradient:Corr ===========================================
//===========================================================================================
    // Read gQso
    gQso = SharedTensor2d(new Tensor2d("Correlation 3-Index TPDM (Q|nn)", nQ_corr, nso_, nso_));
    gQso->read(psio_, PSIF_DFOCC_DENS, true, true);

    // (Q | mu nu)^X
    timer_on("Grad: 3-Index:Corr");    

    //int naux = nQ;
    gradients["3-Index:Corr"] = SharedMatrix(gradients["Nuclear"]->clone());
    gradients["3-Index:Corr"]->set_name("3-Index:Corr Gradient");
    gradients["3-Index:Corr"]->zero();


    // => Sizing <= //
    //int natom = primary_->molecule()->natom();
    //int nso = primary_->nbf();
    //int naux = auxiliary_->nbf();

    boost::shared_ptr<ERISieve> sieve_ = boost::shared_ptr<ERISieve>(new ERISieve(primary_, 0.0));
    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //
    int max_rows;
    max_rows = auxiliary_->nshell();

    // => Block Sizing <= //
    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int nP = auxiliary_->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(auxiliary_->nshell());

    // => Integrals <= //
    boost::shared_ptr<IntegralFactory> rifactory2(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory2->eri(1)));
    }

    // => Temporary Gradients <= //
    std::vector<SharedMatrix> Jtemps2;
    for (int t = 0; t < df_ints_num_threads_; t++) {
         Jtemps2.push_back(SharedMatrix(new Matrix("Jtemp2", natom, 3)));
    }


    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {

        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop  = Pstarts[block+1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop  = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop ).function_index());
        int np = pstop - pstart;

        // > Integrals < //
        int nthread_df = df_ints_num_threads_;
        #pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int P =  PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell_deriv1(P,0,M,N);

            const double* buffer = eri[thread]->buffer();

            int nP = auxiliary_->shell(P).nfunction();
            int cP = auxiliary_->shell(P).ncartesian();
            int aP = auxiliary_->shell(P).ncenter();
            //int oP = auxiliary_->shell(P).function_index() - pstart;
            int oP = auxiliary_->shell(P).function_index();

            int nM = primary_->shell(M).nfunction();
            int cM = primary_->shell(M).ncartesian();
            int aM = primary_->shell(M).ncenter();
            int oM = primary_->shell(M).function_index();

            int nN = primary_->shell(N).nfunction();
            int cN = primary_->shell(N).ncartesian();
            int aN = primary_->shell(N).ncenter();
            int oN = primary_->shell(N).function_index();

            int ncart = cP * cM * cN;
            const double *Px = buffer + 0*ncart;
            const double *Py = buffer + 1*ncart;
            const double *Pz = buffer + 2*ncart;
            const double *Mx = buffer + 3*ncart;
            const double *My = buffer + 4*ncart;
            const double *Mz = buffer + 5*ncart;
            const double *Nx = buffer + 6*ncart;
            const double *Ny = buffer + 7*ncart;
            const double *Nz = buffer + 8*ncart;

            double perm = (M == N ? 1.0 : 2.0);

            double** grad_Jp;
            grad_Jp = Jtemps2[thread]->pointer();

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {

                            //double Ival = 1.0 * perm * gQso->get(p + oP + pstart, (m + oM) * nso + (n + oN));
                            double Ival = 1.0 * perm * gQso->get(p + oP, (m + oM) * nso + (n + oN));
                            grad_Jp[aP][0] += Ival * (*Px);
                            grad_Jp[aP][1] += Ival * (*Py);
                            grad_Jp[aP][2] += Ival * (*Pz);
                            grad_Jp[aM][0] += Ival * (*Mx);
                            grad_Jp[aM][1] += Ival * (*My);
                            grad_Jp[aM][2] += Ival * (*Mz);
                            grad_Jp[aN][0] += Ival * (*Nx);
                            grad_Jp[aN][1] += Ival * (*Ny);
                            grad_Jp[aN][2] += Ival * (*Nz);

                        Px++;
                        Py++;
                        Pz++;
                        Mx++;
                        My++;
                        Mz++;
                        Nx++;
                        Ny++;
                        Nz++;
                    }
                }
            }
        }
    }
    gQso.reset();

    // => Temporary Gradient Reduction <= //
    for (int t = 0; t < df_ints_num_threads_; t++) {
         gradients["3-Index:Corr"]->add(Jtemps2[t]);
    }

    //gradients["3-Index:Corr"]->print_atom_vector();
    timer_off("Grad: 3-Index:Corr");    

//outfile->Printf("\tref_grad is done. \n"); 
}// end 


}} // End Namespaces


