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


#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/sieve.h"
#include "psi4/psi4-dec.h"

#include "defines.h"
#include "dfocc.h"
#include "tensors.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::trans_corr()
{
    // Read SO integrals
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);

    trans_ab = 1;
    if (orb_opt_ == "TRUE" || dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") {
        // Form B(Q,ij)
        timer_on("Form B(Q,ij)");
        b_oo();
        timer_off("Form B(Q,ij)");

        // Form B(Q,ia)
        timer_on("Form B(Q,ia)");
        b_ov();
        timer_off("Form B(Q,ia)");

        // Form B(Q,ab)
        timer_on("Form B(Q,ab)");
        b_vv();
        timer_off("Form B(Q,ab)");
    }

    else {
        // Form B(Q,ij)
        timer_on("Form B(Q,ij)");
        b_ij();
        timer_off("Form B(Q,ij)");

        // Form B(Q,ia)
        timer_on("Form B(Q,ia)");
        b_ia();
        timer_off("Form B(Q,ia)");

        // Form B(Q,ab)
        timer_on("Form B(Q,ab)");
        b_ab();
        timer_off("Form B(Q,ab)");
    }
    bQso.reset();

    // Trans OEI
    timer_on("Trans OEI");
    trans_oei();
    timer_off("Trans OEI");
}

//=======================================================
//          trans for mp2 energy
//=======================================================
void DFOCC::trans_mp2()
{
    // Read SO integrals
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);

    // Form B(Q,ia)
    trans_ab = 0;
    timer_on("Form B(Q,ia)");
    b_ia();
    timer_off("Form B(Q,ia)");
    bQso.reset();
}

//=======================================================
//          DF CC
//=======================================================
void DFOCC::df_corr()
{
    //outfile->Printf("\tComputing DF-BASIS-CC integrals... \n");

    // Read in the basis set informations
    std::shared_ptr<BasisSet> auxiliary_ = get_basisset("DF_BASIS_CC");
    std::shared_ptr<BasisSet> primary_ = get_basisset("ORBITAL");
    
    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
    //auxiliary_->print();

    // Read number of auxilary basis
    nQ = auxiliary_->nbf();

    // Form J^-1/2
    timer_on("Form J");
    formJ(auxiliary_, zero);
    timer_off("Form J");

    // Form B(Q,mu nu)
    timer_on("Form B(Q,munu)");
    b_so(primary_, auxiliary_, zero);
    timer_off("Form B(Q,munu)");

} // end df_corr

//=======================================================
//          form J(P,Q)^-1/2
//=======================================================
void DFOCC::formJ(std::shared_ptr<BasisSet> auxiliary_, std::shared_ptr<BasisSet> zero)
{

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    double **J = block_matrix(nQ, nQ);
    J_mhalf = block_matrix(nQ, nQ);

    // => Integrals <= //
    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_,zero,auxiliary_,zero));
    std::vector<std::shared_ptr<TwoBodyAOInt> > Jint;
    std::vector<const double*> buffer;
    for (int t = 0; t < nthreads; t++) {
        Jint.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri()));
        buffer.push_back(Jint[t]->buffer());
    }

    std::vector<std::pair<int,int> > PQ_pairs;
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int,int>(P,Q));
        }
    }

    #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {

        int P = PQ_pairs[PQ].first;
        int Q = PQ_pairs[PQ].second;

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        Jint[thread]->compute_shell(P,0,Q,0);
        //const double* buffer = Jint[thread]->buffer();

        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();

        int nQ = auxiliary_->shell(Q).nfunction();
        int oQ = auxiliary_->shell(Q).function_index();

        int index = 0;
        for (int p = 0; p < nP; p++) {
            for (int q = 0; q < nQ; q++, ++index) {
                 J[p + oP][q + oQ] = buffer[thread][index];
            }
        }
    }

    /*
    // Create integral factories for the RI basis
    std::shared_ptr<IntegralFactory> rifactory_J(new IntegralFactory(auxiliary_, zero, auxiliary_, zero));
    std::shared_ptr<TwoBodyAOInt> Jint(rifactory_J->eri());

    double **J = block_matrix(nQ, nQ);
    J_mhalf = block_matrix(nQ, nQ);
    const double *Jbuffer = Jint->buffer();

    // Compute J_PQ metric
    int index = 0;
    for (int MU=0; MU < auxiliary_->nshell(); ++MU) {
        int nummu = auxiliary_->shell(MU).nfunction();
        for (int NU=0; NU < auxiliary_->nshell(); ++NU) {
            int numnu = auxiliary_->shell(NU).nfunction();
            Jint->compute_shell(MU, 0, NU, 0);
            index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = auxiliary_->shell(MU).function_index() + mu;
                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = auxiliary_->shell(NU).function_index() + nu;
                    J[omu][onu] = Jbuffer[index];
                }
            }
        }
    }
    */

    // First, diagonalize J
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    int lwork = nQ * 3;
    double* eigval = init_array(nQ);
    double* work = init_array(lwork);
    int status = C_DSYEV('v', 'u', nQ, J[0], nQ, eigval, work, lwork);
    if(status){
        throw PsiException("Diagonalization of J failed", __FILE__, __LINE__);
    }
    free(work);

    // Now J contains the eigenvectors of the original J
    // Copy J to J_copy
    double **J_copy = block_matrix(nQ, nQ);
    C_DCOPY(nQ*nQ, J[0], 1, J_copy[0], 1);

    // Now form J^{-1/2} = U(T)*j^{-1/2}*U,
    // where j^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of J
    for(int i=0; i<nQ; ++i){
        eigval[i] = (eigval[i] < 1.0E-10) ? 0.0 : 1.0 / sqrt(eigval[i]);
        // scale one set of eigenvectors by the diagonal elements j^{-1/2}
        C_DSCAL(nQ, eigval[i], J[i], 1);
    }
    free(eigval);

    // J_mhalf = J_copy(T) * J
    C_DGEMM('t','n', nQ, nQ, nQ, 1.0, J_copy[0], nQ, J[0], nQ, 0.0, J_mhalf[0], nQ);
    free_block(J);
    free_block(J_copy);

    // write J
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_CC Jmhalf <P|Q>", nQ, nQ));
    Jmhalf->set(J_mhalf);
    Jmhalf->write(psio_, PSIF_DFOCC_INTS);
    Jmhalf.reset();

} // end formJ

//=======================================================
//          form b(Q, mu nu)
//=======================================================
void DFOCC::b_so(std::shared_ptr<BasisSet> primary_, std::shared_ptr<BasisSet> auxiliary_, std::shared_ptr<BasisSet> zero)
{
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    double** Ap = block_matrix(nQ, nso2_);
    double** Bp = block_matrix(nQ, nso2_);

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    std::shared_ptr<ERISieve> sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff));
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
    std::shared_ptr<IntegralFactory> rifactory2(new IntegralFactory(auxiliary_, zero, primary_, primary_));
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

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop  = (Pstop == auxiliary_->nshell() ? nQ : auxiliary_->shell(Pstop ).function_index());
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

            int nP = auxiliary_->shell(P).nfunction();
            int oP = auxiliary_->shell(P).function_index();

            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();

            int nN = primary_->shell(N).nfunction();
            int oN = primary_->shell(N).function_index();

            int index = 0;
            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++, index++) {
                         //Bp[p + oP][(m + oM) * nso_ + (n + oN)] = buffer[thread][p * nM * nN + m * nN + n];
                         //Bp[p + oP][(n + oN) * nso_ + (m + oM)] = buffer[thread][p * nM * nN + m * nN + n];
                         Bp[p + oP][(m + oM) * nso_ + (n + oN)] = buffer[thread][index];
                         Bp[p + oP][(n + oN) * nso_ + (m + oM)] = buffer[thread][index];
                    }
                }
            }
        }
    }


    /*
    std::shared_ptr<IntegralFactory> fact(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    std::shared_ptr<TwoBodyAOInt> eri(fact->eri());
    const double* buffer = eri->buffer();

    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        for (int M = 0; M < primary_->nshell(); M++) {
            int nm = primary_->shell(M).nfunction();
            int mstart = primary_->shell(M).function_index();
            for (int N = 0; N < primary_->nshell(); N++) {
                int nn = primary_->shell(N).nfunction();
                int nstart = primary_->shell(N).function_index();

                eri->compute_shell(P,0,M,N);

                for (int p = 0, index = 0; p < np; p++) {
                    for (int m = 0; m < nm; m++) {
                        for (int n = 0; n < nn; n++, index++) {
                            Bp[p + pstart][(m + mstart) * nso_ + (n + nstart)] = buffer[index];
                        }
                    }
                }
            }
        }
    }
    */

    C_DGEMM('N','N', nQ, nso2_, nQ, 1.0, J_mhalf[0], nQ, Bp[0], nso2_, 0.0, Ap[0], nso2_);
    bQso->set(Ap);
    bQso->write(psio_, PSIF_DFOCC_INTS, true, true);
    if (print_ > 3) bQso->print();
    free_block(Bp);
    free_block(J_mhalf);
    free_block(Ap);
    bQso.reset();

    /*
    // Build C(Q, mu nu)
    double** Cp = block_matrix(nQ, nso2_);
    C_DGEMM('N','N', nQ, nso2_, nQ, 1.0, J_mhalf[0], nQ, Ap[0], nso2_, 0.0, Cp[0], nso2_);
    free_block(J_mhalf);
    free_block(Ap);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mn)", nQ, nso2_));
    cQso->set(Cp);
    free_block(Cp);
    cQso->write(psio_, PSIF_DFOCC_INTS);
    cQso.reset();
    */

} // end b_so

//=======================================================
//          form b(Q,ij) : active
//=======================================================
void DFOCC::b_ij()
{
    bQnoA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mI)", nQ, nso_ * naoccA));
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA * naoccA));
    bQnoA->contract(false, false, nQ * nso_, naoccA, nso_, bQso, CaoccA, 1.0, 0.0);
    bQijA->contract233(true, false, naoccA, naoccA, CaoccA, bQnoA, 1.0, 0.0);
    bQnoA.reset();
    bQijA->write(psio_, PSIF_DFOCC_INTS);
    bQijA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQnoB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mi)", nQ, nso_ * naoccB));
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB * naoccB));
    bQnoB->contract(false, false, nQ * nso_, naoccB, nso_, bQso, CaoccB, 1.0, 0.0);
    bQijB->contract233(true, false, naoccB, naoccB, CaoccB, bQnoB, 1.0, 0.0);
    bQnoB.reset();
    bQijB->write(psio_, PSIF_DFOCC_INTS);
    bQijB.reset();
 }

} // end b_ij

//=======================================================
//          form b(Q,ij) : all
//=======================================================
void DFOCC::b_oo()
{
    bQnoA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mO)", nQ, nso_ * noccA));
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OO)", nQ, noccA, noccA));
    bQnoA->contract(false, false, nQ * nso_, noccA, nso_, bQso, CoccA, 1.0, 0.0);
    bQooA->contract233(true, false, noccA, noccA, CoccA, bQnoA, 1.0, 0.0);
    bQnoA.reset();
    bQooA->write(psio_, PSIF_DFOCC_INTS);

    // Form active b(Q,ij)
    bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
    bQijA->form_b_ij(nfrzc, bQooA);
    bQooA.reset();
    bQijA->write(psio_, PSIF_DFOCC_INTS);
    bQijA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQnoB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mo)", nQ, nso_ * noccB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|oo)", nQ, noccB, noccB));
    bQnoB->contract(false, false, nQ * nso_, noccB, nso_, bQso, CoccB, 1.0, 0.0);
    bQooB->contract233(true, false, noccB, noccB, CoccB, bQnoB, 1.0, 0.0);
    bQnoB.reset();
    bQooB->write(psio_, PSIF_DFOCC_INTS);

    // Form active b(Q,ij)
    bQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB));
    bQijB->form_b_ij(nfrzc, bQooB);
    bQooB.reset();
    bQijB->write(psio_, PSIF_DFOCC_INTS);
    bQijB.reset();
 }

} // end b_oo

//=======================================================
//          form b(Q,ia) : active
//=======================================================
void DFOCC::b_ia()
{
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mA)", nQ, nso_ * navirA));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA * navirA));
    bQnvA->contract(false, false, nQ * nso_, navirA, nso_, bQso, CavirA, 1.0, 0.0);
    bQiaA->contract233(true, false, naoccA, navirA, CaoccA, bQnvA, 1.0, 0.0);
    bQiaA->write(psio_, PSIF_DFOCC_INTS);
    bQnvA->write(psio_, PSIF_DFOCC_INTS);
    bQiaA.reset();
    bQnvA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ma)", nQ, nso_ * navirB));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB * navirB));
    bQnvB->contract(false, false, nQ * nso_, navirB, nso_, bQso, CavirB, 1.0, 0.0);
    bQiaB->contract233(true, false, naoccB, navirB, CaoccB, bQnvB, 1.0, 0.0);
    bQiaB->write(psio_, PSIF_DFOCC_INTS);
    bQnvB->write(psio_, PSIF_DFOCC_INTS);
    bQiaB.reset();
    bQnvB.reset();
 }

} // end b_ia

//=======================================================
//          form b(Q,ov) : all
//=======================================================
void DFOCC::b_ov()
{
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mV)", nQ, nso_ * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|OV)", nQ, noccA, nvirA));
    bQnvA->contract(false, false, nQ * nso_, nvirA, nso_, bQso, CvirA, 1.0, 0.0);
    bQovA->contract233(true, false, noccA, nvirA, CoccA, bQnvA, 1.0, 0.0);
    bQovA->write(psio_, PSIF_DFOCC_INTS);
    bQnvA->write(psio_, PSIF_DFOCC_INTS);
    bQnvA.reset();

    // Form active b(Q,ia)
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
    bQiaA->form_b_ia(nfrzc, bQovA);
    bQovA.reset();
    bQiaA->write(psio_, PSIF_DFOCC_INTS);
    bQiaA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mv)", nQ, nso_ * nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ov)", nQ, noccB, nvirB));
    bQnvB->contract(false, false, nQ * nso_, nvirB, nso_, bQso, CvirB, 1.0, 0.0);
    bQovB->contract233(true, false, noccB, nvirB, CoccB, bQnvB, 1.0, 0.0);
    bQovB->write(psio_, PSIF_DFOCC_INTS);
    bQnvB->write(psio_, PSIF_DFOCC_INTS);
    bQnvB.reset();

    // Form active b(Q,ia)
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ia)", nQ, naoccB, navirB));
    bQiaB->form_b_ia(nfrzc, bQovB);
    bQovB.reset();
    bQiaB->write(psio_, PSIF_DFOCC_INTS);
    bQiaB.reset();
 }

} // end b_ov

//=======================================================
//          form b(Q,ab) : active
//=======================================================
void DFOCC::b_ab()
{
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mA)", nQ, nso_ * navirA));
    bQnvA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->contract233(true, false, navirA, navirA, CavirA, bQnvA, 1.0, 0.0);
    bQnvA.reset();
    bQabA->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQabA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ma)", nQ, nso_ * navirB));
    bQnvB->read(psio_, PSIF_DFOCC_INTS);
    bQabB->contract233(true, false, navirB, navirB, CavirB, bQnvB, 1.0, 0.0);
    bQnvB.reset();
    bQabB->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQabB.reset();
 }

} // end b_ab

//=======================================================
//          form b(Q,vv) : all
//=======================================================
void DFOCC::b_vv()
{
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|VV)", nQ, nvirA, nvirA));
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mV)", nQ, nso_ * nvirA));
    bQnvA->read(psio_, PSIF_DFOCC_INTS);
    bQvvA->contract233(true, false, nvirA, nvirA, CvirA, bQnvA, 1.0, 0.0);
    bQnvA.reset();
    bQvvA->write(psio_, PSIF_DFOCC_INTS, true, true);

    // Form active b(Q,ab)
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->form_b_ab(bQvvA);
    bQvvA.reset();
    bQabA->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQabA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mv)", nQ, nso_ * nvirB));
    bQnvB->read(psio_, PSIF_DFOCC_INTS);
    bQvvB->contract233(true, false, nvirB, nvirB, CvirB, bQnvB, 1.0, 0.0);
    bQnvB.reset();
    bQvvB->write(psio_, PSIF_DFOCC_INTS, true, true);

    // Form active b(Q,ab)
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQabB->form_b_ab(bQvvB);
    bQvvB.reset();
    bQabB->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQabB.reset();
 }

} // end b_vv

//=======================================================
//          form c(Q,ij): active
//=======================================================
void DFOCC::c_ij()
{
    cQnoA = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mI)", nQ, nso_ * naoccA));
    cQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|IJ)", nQ, naoccA * naoccA));
    cQnoA->contract(false, false, nQ * nso_, naoccA, nso_, cQso, CaoccA, 1.0, 0.0);
    cQijA->contract233(true, false, naoccA, naoccA, CaoccA, cQnoA, 1.0, 0.0);
    cQnoA.reset();
    cQijA->write(psio_, PSIF_DFOCC_INTS);
    cQijA.reset();

 if (reference_ == "UNRESTRICTED") {
    cQnoB = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mi)", nQ, nso_ * naoccB));
    cQijB = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|ij)", nQ, naoccB * naoccB));
    cQnoB->contract(false, false, nQ * nso_, naoccB, nso_, cQso, CaoccB, 1.0, 0.0);
    cQijB->contract233(true, false, naoccB, naoccB, CaoccB, cQnoB, 1.0, 0.0);
    cQnoB.reset();
    cQijB->write(psio_, PSIF_DFOCC_INTS);
    cQijB.reset();
 }

} // end c_ij

//=======================================================
//          form c(Q,ij): all
//=======================================================
void DFOCC::c_oo()
{
    cQnoA = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mO)", nQ, nso_ * noccA));
    cQooA = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|OO)", nQ, noccA * noccA));
    cQnoA->contract(false, false, nQ * nso_, noccA, nso_, cQso, CoccA, 1.0, 0.0);
    cQooA->contract233(true, false, noccA, noccA, CoccA, cQnoA, 1.0, 0.0);
    cQnoA.reset();
    cQooA->write(psio_, PSIF_DFOCC_INTS);
    cQooA.reset();

 if (reference_ == "UNRESTRICTED") {
    cQnoB = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mo)", nQ, nso_ * noccB));
    cQooB = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|oo)", nQ, noccB * noccB));
    cQnoB->contract(false, false, nQ * nso_, noccB, nso_, cQso, CoccB, 1.0, 0.0);
    cQooB->contract233(true, false, noccB, noccB, CoccB, cQnoB, 1.0, 0.0);
    cQnoB.reset();
    cQooB->write(psio_, PSIF_DFOCC_INTS);
    cQooB.reset();
 }

} // end c_oo

//=======================================================
//          form c(Q,ia) : active
//=======================================================
void DFOCC::c_ia()
{
    cQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mA)", nQ, nso_ * navirA));
    cQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|IA)", nQ, naoccA * navirA));
    cQnvA->contract(false, false, nQ * nso_, navirA, nso_, cQso, CavirA, 1.0, 0.0);
    cQiaA->contract233(true, false, naoccA, navirA, CaoccA, cQnvA, 1.0, 0.0);
    if (trans_ab == 0) cQnvA.reset();
    cQiaA->write(psio_, PSIF_DFOCC_INTS);
    cQiaA.reset();

 if (reference_ == "UNRESTRICTED") {
    cQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|ma)", nQ, nso_ * navirB));
    cQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|ia)", nQ, naoccB * navirB));
    cQnvB->contract(false, false, nQ * nso_, navirB, nso_, cQso, CavirB, 1.0, 0.0);
    cQiaB->contract233(true, false, naoccB, navirB, CaoccB, cQnvB, 1.0, 0.0);
    if (trans_ab == 0) cQnvB.reset();
    cQiaB->write(psio_, PSIF_DFOCC_INTS);
    cQiaB.reset();
 }

} // end c_ia

//=======================================================
//          form c(Q,ia) : all
//=======================================================
void DFOCC::c_ov()
{
    cQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mV)", nQ, nso_ * nvirA));
    cQovA = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|OV)", nQ, noccA * nvirA));
    cQnvA->contract(false, false, nQ * nso_, nvirA, nso_, cQso, CvirA, 1.0, 0.0);
    cQovA->contract233(true, false, noccA, nvirA, CoccA, cQnvA, 1.0, 0.0);
    if (trans_ab == 0) cQnvA.reset();
    cQovA->write(psio_, PSIF_DFOCC_INTS);
    cQovA.reset();

 if (reference_ == "UNRESTRICTED") {
    cQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|mv)", nQ, nso_ * nvirB));
    cQovB = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|ov)", nQ, noccB * nvirB));
    cQnvB->contract(false, false, nQ * nso_, nvirB, nso_, cQso, CvirB, 1.0, 0.0);
    cQovB->contract233(true, false, noccB, nvirB, CoccB, cQnvB, 1.0, 0.0);
    if (trans_ab == 0) cQnvB.reset();
    cQovB->write(psio_, PSIF_DFOCC_INTS);
    cQovB.reset();
 }

} // end c_ov

//=======================================================
//          form c(Q,ab) : active
//=======================================================
void DFOCC::c_ab()
{
    cQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|AB)", nQ, navirA * navirA));
    cQabA->contract233(true, false, navirA, navirA, CavirA, cQnvA, 1.0, 0.0);
    cQnvA.reset();
    cQabA->write(psio_, PSIF_DFOCC_INTS);
    cQabA.reset();

 if (reference_ == "UNRESTRICTED") {
    cQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|ab)", nQ, navirB * navirB));
    cQabB->contract233(true, false, navirB, navirB, CavirB, cQnvB, 1.0, 0.0);
    cQnvB.reset();
    cQabB->write(psio_, PSIF_DFOCC_INTS);
    cQabB.reset();
 }

} // end c_ab

//=======================================================
//          form c(Q,ab) : all
//=======================================================
void DFOCC::c_vv()
{
    cQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|VV)", nQ, nvirA * nvirA));
    cQvvA->contract233(true, false, nvirA, nvirA, CvirA, cQnvA, 1.0, 0.0);
    cQnvA.reset();
    cQvvA->write(psio_, PSIF_DFOCC_INTS);
    cQvvA.reset();

 if (reference_ == "UNRESTRICTED") {
    cQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC C (Q|VV)", nQ, nvirB * nvirB));
    cQvvB->contract233(true, false, nvirB, nvirB, CvirB, cQnvB, 1.0, 0.0);
    cQnvB.reset();
    cQvvB->write(psio_, PSIF_DFOCC_INTS);
    cQvvB.reset();
 }

} // end c_vv

//=======================================================
//          Trans OEI ints
//=======================================================
void DFOCC::trans_oei()
{
    // Alpha
    HmoA->transform(Hso, CmoA);
    if (print_ > 2) HmoA->print();
    // Blocks
    HooA->form_oo(HmoA);
    HvoA->form_vo(HmoA);
    HovA->form_ov(HmoA);
    //HovA = HvoA->transpose();
    HvvA->form_vv(noccA, HmoA);

 if (reference_ == "UNRESTRICTED") {
    HmoB->transform(Hso, CmoB);
    if (print_ > 2) HmoB->print();
    // Blocks
    HooB->form_oo(HmoB);
    HvoB->form_vo(HmoB);
    HovB->form_ov(HmoB);
    //HovB = HvoB->transpose();
    HvvB->form_vv(noccB, HmoB);
 } // uhf

} // end trans_oei

//=======================================================
//       Form non-zero DF ints: Experimental
//=======================================================
void DFOCC::b_so_non_zero()
{

    // defs
    SharedTensor2d K, L, J, G;
    int nmn_nz, syc;

    // Read SO integrals
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);


    ndf_nz = 0;
    #pragma omp parallel for
    for(int Q = 0 ; Q < nQ; ++Q){
        for(int m = 0 ; m < nso_; ++m){
            for(int n = 0 ; n < nso_; ++n){
                int mn = n + (m * nso_);
		if (fabs(bQso->get(Q,mn)) > int_cutoff_) ndf_nz++;
            }
        }
    }

    int ndf_ao = 0;
    ndf_ao = nQ*nso_*nso_;
    double perct_ = 0.0;
    perct_ = (double)ndf_nz / (double)ndf_ao;
    perct_ *= 100;
    //outfile->Printf("\tNumber of AO-basis DF-CC integrals          : %3d\n", ndf_ao);
    //outfile->Printf("\tNumber of non-zero AO-basis DF-CC integrals : %3d\n", ndf_nz);
    //outfile->Printf("\tPercent of non-zero DF-CC integrals         : %2.2f\n", perct_);

    K = SharedTensor2d(new Tensor2d("DF_BASIS_CC NONZERO B (Q|mn)", ndf_nz, 1));
    ndf_nz = 0;
    #pragma omp parallel for
    for(int Q = 0 ; Q < nQ; ++Q){
        for(int m = 0 ; m < nso_; ++m){
            for(int n = 0 ; n < nso_; ++n){
                int mn = n + (m * nso_);

		if (fabs(bQso->get(Q,mn)) > int_cutoff_) {
		    K->set(ndf_nz, 0, bQso->get(Q,mn));
		    ndf_nz++;
		    //if (m >= n) outfile->Printf("\tQ, m, n: %3d %3d %3d \n", Q, m, n);
		}

            }
        }
    }
    //K->write(psio_, PSIF_DFOCC_INTS);


    L = SharedTensor2d(new Tensor2d("DF_BASIS_CC AO-Basis (mn|ls)", nso_, nso_, nso_, nso_));
    L->gemm(true, false, bQso, bQso, 1.0, 0.0);
    //L->print();

    /*
    ndf_nz = 0;
    #pragma omp parallel for
    for(int m = 0 ; m < nso_; ++m){
        for(int n = 0 ; n < nso_; ++n){
            int mn = n + (m * nso_);
            for(int l = 0 ; l < nso_; ++l){
                for(int s = 0 ; s < nso_; ++s){
                    int ls = s + (l * nso_);
		    if (fabs(L->get(mn,ls)) > int_cutoff_) ndf_nz++;
		}
            }
        }
    }
    L.reset();

    ndf_ao = nso_*nso_*nso_*nso_;
    perct_ = 0.0;
    perct_ = (double)ndf_nz / (double)ndf_ao;
    perct_ *= 100;
    outfile->Printf("\tNumber of (mn|ls) integrals                 : %3d\n", ndf_ao);
    outfile->Printf("\tNumber of non-zero (mn|ls) integrals        : %3d\n", ndf_nz);
    outfile->Printf("\tPercent of non-zero (mn|ls) integrals       : %2.2f\n", perct_);
    */

    ndf_nz = 0;
    #pragma omp parallel for
    for(int m = 0 ; m < nso_; ++m){
        for(int n = 0 ; n <=m; ++n){
            int mn2 = n + (m * nso_);
	    int mn = index2(m,n);
            for(int l = 0 ; l < nso_; ++l){
                for(int s = 0 ; s <=l; ++s){
                    int ls2 = s + (l * nso_);
	            int ls = index2(l,s);

		    if (mn >= ls ) {
		        if (fabs(L->get(mn2,ls2)) > int_cutoff_) ndf_nz++;
		    }

		}
            }
        }
    }
    L.reset();

    ndf_ao = ntri_so*(ntri_so+1) / 2;
    perct_ = 0.0;
    perct_ = (double)ndf_nz / (double)ndf_ao;
    perct_ *= 100;
    outfile->Printf("\tNumber of (mn|ls) integrals                 : %3d\n", ndf_ao);
    outfile->Printf("\tNumber of non-zero (mn|ls) integrals        : %3d\n", ndf_nz);
    outfile->Printf("\tPercent of non-zero (mn|ls) integrals       : %2.2f\n", perct_);


    // screening
    G = SharedTensor2d(new Tensor2d("Presecreening (mn|mn)", nso_, nso_));

    #pragma omp parallel for
    for(int m = 0 ; m < nso_; ++m){
        for(int n = 0 ; n < nso_; ++n){
            int mn = n + (m * nso_);
	    double sum = 0.0;
            for(int Q = 0 ; Q < nQ; ++Q){
                sum += bQso->get(Q,mn) * bQso->get(Q,mn);
            }
	    double value = sqrt(sum);
	    G->set(m, n, value);
        }
    }
    //G->print();


    ndf_nz = 0;
    #pragma omp parallel for
    for(int m = 0 ; m < nso_; ++m){
        for(int n = 0 ; n <=m; ++n){
	    int mn = index2(m,n);
            for(int l = 0 ; l < nso_; ++l){
                for(int s = 0 ; s <=l; ++s){
	            int ls = index2(l,s);

		    if (mn >= ls ) {
		        if ( fabs(G->get(m,n)*G->get(l,s)) > int_cutoff_) ndf_nz++;
		    }

		}
            }
        }
    }
    G.reset();

    perct_ = 0.0;
    perct_ = (double)ndf_nz / (double)ndf_ao;
    perct_ *= 100;
    outfile->Printf("\tNumber of (mn|ls) integrals                 : %3d\n", ndf_ao);
    outfile->Printf("\tNumber of prescreened (mn|ls) integrals     : %3d\n", ndf_nz);
    outfile->Printf("\tPercent of non-zero (mn|ls) integrals       : %2.2f\n", perct_);

    ndf_nz = 0;
    #pragma omp parallel for
    for(int m = 0 ; m < nso_; ++m){
        for(int n = 0 ; n <=m; ++n){
	    int mn = index2(m,n);
            for(int l = 0 ; l < nso_; ++l){
                for(int s = 0 ; s <=l; ++s){
	            int ls = index2(l,s);

		    if (mn >= ls ) {
		        if ( fabs(Sso->get(m,n)*Sso->get(l,s)) > int_cutoff_) ndf_nz++;
		    }

		}
            }
        }
    }

    perct_ = 0.0;
    perct_ = (double)ndf_nz / (double)ndf_ao;
    perct_ *= 100;
    outfile->Printf("\tNumber of (mn|ls) integrals                 : %3d\n", ndf_ao);
    outfile->Printf("\tNumber of overlap-prescreened (mn|ls)       : %3d\n", ndf_nz);
    outfile->Printf("\tPercent of non-zero (mn|ls) integrals       : %2.2f\n", perct_);

    // free
    bQso.reset();
    K.reset();

} // end b_so_non_zero

//=======================================================
//       Form LDL ABCD ints: Experimental
//=======================================================
void DFOCC::ldl_abcd_ints()
{
    timer_on("LDL <AB|CD>");

    // Variables
    nQ_cd = 1;
    int n = navirA * navirA;
    double dmax = 0.0;
    double dmin = 0.0;
    int maxQ = n;
    int Q = 0;
    int ndiag = 0;

    SharedTensor1d D, D2, L1, U1, R1;
    SharedTensor2d R, L, L2, Dmat, U, LU, J, K;
    SharedTensor1i o2n, n2o;
    SharedTensor1i pair_to_idx1, pair_to_idx2;

    // Title
    outfile->Printf("\n\tGenerating LDL factors ...\n");
    outfile->Printf("\tLDL decomposition threshold: %8.2le\n", tol_ldl);

    outfile->Printf( "\n\t          LDL for <AB|CD> \n");
    outfile->Printf( "\t   ------------------------------ \n");
    outfile->Printf( "\tIter      max(|D_Q|)        # of LDL vectors  \n");
    outfile->Printf( "\t----    ---------------    ------------------ \n");

    // Pair to idx mapping
    pair_to_idx1 = SharedTensor1i(new Tensor1i("AB -> A", n));
    pair_to_idx2 = SharedTensor1i(new Tensor1i("AB -> B", n));
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        for(int b = 0 ; b < navirA; ++b){
            int ab = ab_idxAA->get(a,b);
	    pair_to_idx1->set(ab,a);
	    pair_to_idx2->set(ab,b);
	}
    }

    /*
    // compute diagonals (AB|AB)
    D = SharedTensor1d(new Tensor1d("D", n));
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        for(int b = 0 ; b < navirA; ++b){
            int ab = ab_idxAA->get(a,b);
	    double sum = 0.0;
            for(int P = 0 ; P < nQ; ++P){
		sum += bQabA->get(P,ab) * bQabA->get(P,ab);
	    }
	    D->set(ab,sum);
	}
    }
    */

    // compute diagonals <AB|AB>
    D = SharedTensor1d(new Tensor1d("D", n));
    #pragma omp parallel for
    for(int a = 0 ; a < navirA; ++a){
        int aa = ab_idxAA->get(a,a);
        for(int b = 0 ; b < navirA; ++b){
            int bb = ab_idxAA->get(b,b);
            int ab = ab_idxAA->get(a,b);
	    double sum = 0.0;
            for(int P = 0 ; P < nQ; ++P){
		sum += bQabA->get(P,aa) * bQabA->get(P,bb);
	    }
	    D->set(ab,sum);
	}
    }

    // Count number of non-zero diagonals
    for(int i=0; i< n; i++) {
        if (fabs(D->get(i)) > tol_ldl) ndiag++;
    }

    // Initialize mapping arrays
    o2n = SharedTensor1i(new Tensor1i("Old -> New", n));
    n2o = SharedTensor1i(new Tensor1i("New -> Old", n));
    #pragma omp parallel for
    for(int i=0; i < n; i++) {
	o2n->set(i,i);
	n2o->set(i,i);
    }

    // Descending ordering for D
    for(int i=0; i< (n-1); i++) {
        for(int j = (i+1); j < n; j++) {
            if (fabs(D->get(i)) < fabs(D->get(j))) {
                double temp= D->get(i);
         	D->set(i,D->get(j));
		D->set(j,temp);
                int i_org = n2o->get(i);
                int j_org = n2o->get(j);
		o2n->set(i_org,j);
		o2n->set(j_org,i);
		n2o->set(i,j_org);
		n2o->set(j,i_org);
	    }
	}
    }

    // RMS
    dmax = fabs(D->get(Q));
    dmin = fabs(D->get(ndiag-1));
    outfile->Printf("\tNumber of complete LDL factors:   %5li\n",n);
    outfile->Printf("\tEstimated number of LDL factors:   %5li\n",ndiag);
    outfile->Printf("\tmax(|D_Q|) =%12.8f\n",dmax);
    outfile->Printf("\tmin(|D_Q|) =%12.8f\n",dmin);
    outfile->Printf("\t%3d     %12.8f          %3d\n",Q,dmax,nQ_cd);

    // compute the off-diagonal
    R1 = SharedTensor1d(new Tensor1d("R1", n));
    int j0 = n2o->get(Q);
    int c0 = pair_to_idx1->get(j0);
    int d0 = pair_to_idx2->get(j0);
    #pragma omp parallel for
    for(int i = Q+1; i < n; i++) {
        int ab = n2o->get(i);
        int a0 = pair_to_idx1->get(ab);
        int b0 = pair_to_idx2->get(ab);
        int ac = ab_idxAA->get(a0,c0);
        int bd = ab_idxAA->get(b0,d0);
	double sum = 0.0;
        for(int P = 0 ; P < nQ; ++P){
            sum += bQabA->get(P,ac) * bQabA->get(P,bd);
	}
	R1->set(i,sum);
    }

    // Compute the first L vector
    L1 = SharedTensor1d(new Tensor1d("L1", n));
    L1->set(0,1.0);
    #pragma omp parallel for
    for(int i = Q+1; i < n; i++) {
	double value = R1->get(i)/D->get(Q);
	L1->set(i,value);
    }

    // Form the global L matrix
    L = SharedTensor2d(new Tensor2d("L <AB|Q>", n, nQ_cd));
    #pragma omp parallel for
    for(int i = 0; i < n; i++) {
	L->set(i,0,L1->get(i));
    }

//==========================================================================================
//========================= Head of the Loop ===============================================
//==========================================================================================
    // Start Iterations
    do
    {
        // increment Q
        Q++;

	// Update the diagonal vector
        #pragma omp parallel for
        for(int j = Q; j < n; j++) {
	    double value = L->get(j,Q-1) * L->get(j,Q-1) * D->get(Q-1);
	    D->subtract(j,value);
	}

	/*
    	// Revert L to the original ordering
    	L2 = SharedTensor2d(new Tensor2d("L-copy", n, nQ_cd));
    	L2->copy(L);
        #pragma omp parallel for
    	for(int i = 0; i < n; i++) {
	    int i_new = o2n->get(i);
            for(int j = 0; j < nQ_cd; j++) {
	    	L->set(i,j,L2->get(i_new,j));
	     }
    	}
    	L2.reset();

        // Descending ordering for D
        for(int i=Q; i < (n-1); i++) {
            for(int j = (i+1); j < n; j++) {
                if (fabs(D->get(i)) < fabs(D->get(j))) {
                    double temp = D->get(i);
         	    D->set(i,D->get(j));
		    D->set(j,temp);
                    int i_org = n2o->get(i);
                    int j_org = n2o->get(j);
		    o2n->set(i_org,j);
		    o2n->set(j_org,i);
		    n2o->set(i,j_org);
		    n2o->set(j,i_org);
	        }
	    }
        }
	*/

        // Choose the dmax
        dmax = fabs(D->get(Q));
        #pragma omp parallel for
        for(int i = Q+1; i < n; i++) {
	    if (fabs(D->get(i)) > dmax) {
                dmax = fabs(D->get(i));
	    }
        }

        // Print
        outfile->Printf("\t%3d     %12.8f          %3d\n",Q,dmax,nQ_cd);

	/*
    	// Reorder L to the new ordering
        L2 = SharedTensor2d(new Tensor2d("L-copy", n, nQ_cd));
        L2->copy(L);
        #pragma omp parallel for
        for(int i = 0; i < n; i++) {
	    int i_old = n2o->get(i);
            for(int j = 0; j < nQ_cd; j++) {
	        L->set(i,j,L2->get(i_old,j));
	    }
        }
        L2.reset();
	*/

	// check dmax
	if (dmax <= tol_ldl) break;

        // Form U1
	// U[Q](P) = L(P,Q)*D(P)
        U1 = SharedTensor1d(new Tensor1d("U1", nQ_cd));
        #pragma omp parallel for
        for(int P = 0; P < nQ_cd; P++) {
	    U1->set(P,L->get(Q,P)*D->get(P));
        }

        // Compute the off-diagonal: Part-1
	// R(ab,Q) = <ab|Q>
	R1->zero();
        int j0 = n2o->get(Q);
        int c0 = pair_to_idx1->get(j0);
        int d0 = pair_to_idx2->get(j0);
        #pragma omp parallel for
        for(int i = Q+1; i < n; i++) {
	    if (fabs(D->get(i)) * fabs(D->get(Q)) > tol_ldl) {
            int ab = n2o->get(i);
            int a0 = pair_to_idx1->get(ab);
            int b0 = pair_to_idx2->get(ab);
            int ac = ab_idxAA->get(a0,c0);
            int bd = ab_idxAA->get(b0,d0);
	    double sum = 0.0;
            for(int P = 0 ; P < nQ; ++P){
                sum += bQabA->get(P,ac) * bQabA->get(P,bd);
	    }
	    R1->set(i,sum);
	    }
        }

        // Compute the off-diagonal: Part-2
	// R(ab,Q) -= \sum_{P=0 to Q-1} L(ab,P) * L(Q,P) * D(P)
	R1->gemv(false, L, U1, -1.0, 1.0);
        U1.reset();

        // Compute the next L vector
	L1->zero();
	L1->set(Q,1.0); // set diagonal to 1
        #pragma omp parallel for
        for(int i = Q+1; i < n; i++) {
	    if (fabs(D->get(Q)) > tol_ldl) {
	        double value = R1->get(i)/D->get(Q);
	        L1->set(i,value);
	    }
        }

        // Form the global L matrix
        nQ_cd++;
        L2 = SharedTensor2d(new Tensor2d("New L", n, nQ_cd));
	// copy previous L
        #pragma omp parallel for
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < nQ_cd-1; j++) {
		L2->set(i,j,L->get(i,j));
	    }
	}
        L.reset();
        L = SharedTensor2d(new Tensor2d("L <AB|Q>", n, nQ_cd));
        L->copy(L2);
        L2.reset();
	// add L1
        #pragma omp parallel for
        for(int i = 0; i < n; i++) {
	    L->set(i,Q,L1->get(i));
	}

    }
    while(Q < (maxQ-1));
    outfile->Printf("\tIteratons were done.\n");
    outfile->Printf("\tmax(|D_Q|) =%12.8f\n",dmax);
    outfile->Printf("\tNumber of computed LDL factors:   %5li\n",nQ_cd);
//==========================================================================================
//========================= End of the Loop ================================================
//==========================================================================================
    // Form U
    U = SharedTensor2d(new Tensor2d("U <Q|CD>", nQ_cd, n));
    #pragma omp parallel for
    for(int Q = 0; Q < nQ_cd; Q++) {
        for(int i = 0; i < n; i++) {
	    U->set(Q,i,L->get(i,Q)*D->get(Q));
	}
    }

    // Revert L to the original ordering
    L2 = SharedTensor2d(new Tensor2d("L-copy", n, nQ_cd));
    L2->copy(L);
    #pragma omp parallel for
    for(int i = 0; i < n; i++) {
	int i_new = o2n->get(i);
        for(int j = 0; j < nQ_cd; j++) {
	    L->set(i,j,L2->get(i_new,j));
	}
    }
    L2.reset();

    // Revert U to the original ordering
    L2 = SharedTensor2d(new Tensor2d("L-copy", nQ_cd, n));
    L2->copy(U);
    #pragma omp parallel for
    for(int i = 0; i < nQ_cd; i++) {
        for(int j = 0; j < n; j++) {
	    int j_new = o2n->get(j);
	    U->set(i,j,L2->get(i,j_new));
	}
    }
    L2.reset();

    // Write L matrix
    L->write(psio_, PSIF_DFOCC_INTS);
    U->write(psio_, PSIF_DFOCC_INTS);

    /*
    // Verification
    LU = SharedTensor2d(new Tensor2d("LU", n, n));
    LU->gemm(false,false,L,U,1.0,0.0);

    // Form exact <AB|CD>
    J = SharedTensor2d(new Tensor2d("J (AC|BD)", navirA, navirA, navirA, navirA));
    J->gemm(true, false, bQabA, bQabA, 1.0, 0.0);
    K = SharedTensor2d(new Tensor2d("K <AB|CD>", navirA, navirA, navirA, navirA));
    K->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Final Results
    L->print();
    U->print();
    LU->print();
    LU.reset();
    K->print();
    K.reset();
    */

    // reset
    D.reset();
    U.reset();
    L.reset();
    L1.reset();
    R1.reset();
    o2n.reset();
    n2o.reset();
    pair_to_idx1.reset();
    pair_to_idx2.reset();

    timer_off("LDL <AB|CD>");

} // end ldl_abcd_ints

//=======================================================
//       Form LDL PQRS ints: Experimental
//=======================================================
void DFOCC::ldl_pqrs_ints(int dim1, int dim2, SharedTensor2d &bQ)
{
    timer_on("LDL <PQ|RS>");

    // Variables
    nQ_cd = 1;
    int n = dim1 * dim2;
    double dmax = 0.0;
    double dmin = 0.0;
    int maxQ = n;
    int Q = 0;
    int ndiag = 0;

    SharedTensor1d D, D2, L1, U1, R1;
    SharedTensor2d R, L, L2, Dmat, U, LU, J, K;
    SharedTensor1i o2n, n2o;
    SharedTensor1i pair_to_idx1, pair_to_idx2;

    // Title
    outfile->Printf("\n\tGenerating LDL factors ...\n");
    outfile->Printf("\tLDL decomposition threshold: %8.2le\n", tol_ldl);

    // Pair to idx mapping
    pair_to_idx1 = SharedTensor1i(new Tensor1i("AB -> A", n));
    pair_to_idx2 = SharedTensor1i(new Tensor1i("AB -> B", n));
    #pragma omp parallel for
    for(int a = 0 ; a < dim1; ++a){
        for(int b = 0 ; b < dim2; ++b){
            int ab = (a*dim2) + b;
	    pair_to_idx1->set(ab,a);
	    pair_to_idx2->set(ab,b);
	}
    }

    /*
    // compute diagonals (PQ|PQ)
    D = SharedTensor1d(new Tensor1d("D", n));
    #pragma omp parallel for
    for(int a = 0 ; a < dim1; ++a){
        for(int b = 0 ; b < dim2; ++b){
            int ab = (a*dim2) + b;
	    double sum = 0.0;
            for(int P = 0 ; P < nQ; ++P){
		sum += bQ->get(P,ab) * bQ->get(P,ab);
	    }
	    D->set(ab,sum);
	}
    }
    */

    // compute diagonals <PQ|PQ>
    D = SharedTensor1d(new Tensor1d("D", n));
    #pragma omp parallel for
    for(int a = 0 ; a < dim1; ++a){
        int aa = (a*dim2) + a;
        for(int b = 0 ; b < dim2; ++b){
            int bb = (b*dim2) + b;
            int ab = (a*dim2) + b;
	    double sum = 0.0;
            for(int P = 0 ; P < nQ; ++P){
		sum += bQ->get(P,aa) * bQ->get(P,bb);
	    }
	    D->set(ab,sum);
	}
    }

    // Count number of non-zero diagonals
    for(int i=0; i< n; i++) {
        if (fabs(D->get(i)) > tol_ldl) ndiag++;
    }

    // Initialize mapping arrays
    o2n = SharedTensor1i(new Tensor1i("Old -> New", n));
    n2o = SharedTensor1i(new Tensor1i("New -> Old", n));
    #pragma omp parallel for
    for(int i=0; i < n; i++) {
	o2n->set(i,i);
	n2o->set(i,i);
    }

    // Descending ordering for D
    for(int i=0; i< (n-1); i++) {
        for(int j = (i+1); j < n; j++) {
            if (fabs(D->get(i)) < fabs(D->get(j))) {
                double temp= D->get(i);
         	D->set(i,D->get(j));
		D->set(j,temp);
                int i_org = n2o->get(i);
                int j_org = n2o->get(j);
		o2n->set(i_org,j);
		o2n->set(j_org,i);
		n2o->set(i,j_org);
		n2o->set(j,i_org);
	    }
	}
    }

    // RMS
    dmax = fabs(D->get(Q));
    dmin = fabs(D->get(ndiag-1));
    outfile->Printf("\tNumber of complete LDL factors:   %5li\n",n);
    outfile->Printf("\tEstimated number of LDL factors:   %5li\n",ndiag);
    outfile->Printf("\tmax(|D_Q|) =%12.8f\n",dmax);
    outfile->Printf("\tmin(|D_Q|) =%12.8f\n",dmin);
    //outfile->Printf("\t%3d     %12.8f     %3d\n",Q,dmax,nQ_cd);

    // compute the off-diagonal
    R1 = SharedTensor1d(new Tensor1d("R1", n));
    int j0 = n2o->get(Q);
    int c0 = pair_to_idx1->get(j0);
    int d0 = pair_to_idx2->get(j0);
    #pragma omp parallel for
    for(int i = Q+1; i < n; i++) {
        int ab = n2o->get(i);
        int a0 = pair_to_idx1->get(ab);
        int b0 = pair_to_idx2->get(ab);
        int ac = (a0*dim2) + c0;
        int bd = (b0*dim2) + d0;
	double sum = 0.0;
        for(int P = 0 ; P < nQ; ++P){
            sum += bQ->get(P,ac) * bQ->get(P,bd);
	}
	R1->set(i,sum);
    }

    // Compute the first L vector
    L1 = SharedTensor1d(new Tensor1d("L1", n));
    L1->set(0,1.0);
    #pragma omp parallel for
    for(int i = Q+1; i < n; i++) {
	double value = R1->get(i)/D->get(Q);
	L1->set(i,value);
    }

    // Form the global L matrix
    L = SharedTensor2d(new Tensor2d("L <AB|Q>", n, nQ_cd));
    #pragma omp parallel for
    for(int i = 0; i < n; i++) {
	L->set(i,0,L1->get(i));
    }

//==========================================================================================
//========================= Head of the Loop ===============================================
//==========================================================================================
    // Start Iterations
    do
    {
        // increment Q
        Q++;

	// Update the diagonal vector
        #pragma omp parallel for
        for(int j = Q; j < n; j++) {
	    double value = L->get(j,Q-1) * L->get(j,Q-1) * D->get(Q-1);
	    D->subtract(j,value);
	}

    	// Revert L to the original ordering
    	L2 = SharedTensor2d(new Tensor2d("L-copy", n, nQ_cd));
    	L2->copy(L);
        #pragma omp parallel for
    	for(int i = 0; i < n; i++) {
	    int i_new = o2n->get(i);
            for(int j = 0; j < nQ_cd; j++) {
	    	L->set(i,j,L2->get(i_new,j));
	     }
    	}
    	L2.reset();

        // Descending ordering for D
        for(int i=Q; i < (n-1); i++) {
            for(int j = (i+1); j < n; j++) {
                if (fabs(D->get(i)) < fabs(D->get(j))) {
                    double temp = D->get(i);
         	    D->set(i,D->get(j));
		    D->set(j,temp);
                    int i_org = n2o->get(i);
                    int j_org = n2o->get(j);
		    o2n->set(i_org,j);
		    o2n->set(j_org,i);
		    n2o->set(i,j_org);
		    n2o->set(j,i_org);
	        }
	    }
        }

        // Choose the dmax
        dmax = fabs(D->get(Q));
        #pragma omp parallel for
        for(int i = Q+1; i < n; i++) {
	    if (fabs(D->get(i)) > dmax) {
                dmax = fabs(D->get(i));
	    }
        }

        // Print
        //outfile->Printf("\t%3d     %12.8f          %3d\n",Q,dmax,nQ_cd);

    	// Reorder L to the new ordering
        L2 = SharedTensor2d(new Tensor2d("L-copy", n, nQ_cd));
        L2->copy(L);
        #pragma omp parallel for
        for(int i = 0; i < n; i++) {
	    int i_old = n2o->get(i);
            for(int j = 0; j < nQ_cd; j++) {
	        L->set(i,j,L2->get(i_old,j));
	    }
        }
        L2.reset();

	// check dmax
	if (dmax <= tol_ldl) break;

        // Form U1
	// U[Q](P) = L(P,Q)*D(P)
        U1 = SharedTensor1d(new Tensor1d("U1", nQ_cd));
        for(int P = 0; P < nQ_cd; P++) {
	    U1->set(P,L->get(Q,P)*D->get(P));
        }

        // Compute the off-diagonal: Part-1
	// R(ab,Q) = <ab|Q>
        int j0 = n2o->get(Q);
        int c0 = pair_to_idx1->get(j0);
        int d0 = pair_to_idx2->get(j0);
        #pragma omp parallel for
        for(int i = Q+1; i < n; i++) {
            int ab = n2o->get(i);
            int a0 = pair_to_idx1->get(ab);
            int b0 = pair_to_idx2->get(ab);
            int ac = (a0*dim2) + c0;
            int bd = (b0*dim2) + d0;
	    double sum = 0.0;
            for(int P = 0 ; P < nQ; ++P){
                sum += bQ->get(P,ac) * bQ->get(P,bd);
	    }
	    R1->set(i,sum);
        }

        // Compute the off-diagonal: Part-2
	// R(ab,Q) -= \sum_{P=0 to Q-1} L(ab,P) * L(Q,P) * D(P)
	R1->gemv(false, L, U1, -1.0, 1.0);
        U1.reset();

        // Compute the next L vector
	L1->zero();
	L1->set(Q,1.0); // set diagonal to 1
        #pragma omp parallel for
        for(int i = Q+1; i < n; i++) {
	    if (fabs(D->get(Q)) > tol_ldl) {
	        double value = R1->get(i)/D->get(Q);
	        L1->set(i,value);
	    }
        }

        // Form the global L matrix
        nQ_cd++;
        L2 = SharedTensor2d(new Tensor2d("New L", n, nQ_cd));
	// copy previous L
        #pragma omp parallel for
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < nQ_cd-1; j++) {
		L2->set(i,j,L->get(i,j));
	    }
	}
        L.reset();
	// add L1
        #pragma omp parallel for
        for(int i = 0; i < n; i++) {
	    L2->set(i,Q,L1->get(i));
	}
        L = SharedTensor2d(new Tensor2d("L <AB|Q>", n, nQ_cd));
        L->copy(L2);
        L2.reset();

    }
    while(Q < (maxQ-1));
    outfile->Printf("\tIteratons were done.\n");
    outfile->Printf("\tmax(|D_Q|) =%12.8f\n",dmax);
    outfile->Printf("\tNumber of computed LDL factors:   %5li\n",nQ_cd);
//==========================================================================================
//========================= End of the Loop ================================================
//==========================================================================================
    // Form U
    U = SharedTensor2d(new Tensor2d("U <Q|CD>", nQ_cd, n));
    #pragma omp parallel for
    for(int Q = 0; Q < nQ_cd; Q++) {
        for(int i = 0; i < n; i++) {
	    U->set(Q,i,L->get(i,Q)*D->get(Q));
	}
    }

    // Revert L to the original ordering
    L2 = SharedTensor2d(new Tensor2d("L-copy", n, nQ_cd));
    L2->copy(L);
    #pragma omp parallel for
    for(int i = 0; i < n; i++) {
	int i_new = o2n->get(i);
        for(int j = 0; j < nQ_cd; j++) {
	    L->set(i,j,L2->get(i_new,j));
	}
    }
    L2.reset();

    // Revert U to the original ordering
    L2 = SharedTensor2d(new Tensor2d("L-copy", nQ_cd, n));
    L2->copy(U);
    #pragma omp parallel for
    for(int i = 0; i < nQ_cd; i++) {
        for(int j = 0; j < n; j++) {
	    int j_new = o2n->get(j);
	    U->set(i,j,L2->get(i,j_new));
	}
    }
    L2.reset();

    // Write L matrix
    L->write(psio_, PSIF_DFOCC_INTS);
    U->write(psio_, PSIF_DFOCC_INTS);

    /*
    // Verification
    LU = SharedTensor2d(new Tensor2d("LU", n, n));
    LU->gemm(false,false,L,U,1.0,0.0);

    // Form exact <AB|CD>
    J = SharedTensor2d(new Tensor2d("J (AC|BD)", navirA, navirA, navirA, navirA));
    J->gemm(true, false, bQabA, bQabA, 1.0, 0.0);
    K = SharedTensor2d(new Tensor2d("K <AB|CD>", navirA, navirA, navirA, navirA));
    K->sort(1324, J, 1.0, 0.0);
    J.reset();

    // Final Results
    L->print();
    U->print();
    LU->print();
    LU.reset();
    K->print();
    K.reset();
    */

    // reset
    D.reset();
    U.reset();
    L.reset();
    L1.reset();
    R1.reset();
    o2n.reset();
    n2o.reset();
    pair_to_idx1.reset();
    pair_to_idx2.reset();

    timer_off("LDL <PQ|RS>");

} // end ldl_pqrs_ints

//=======================================================
//       Form CD (MN|LS) ints
//=======================================================
void DFOCC::cd_aob_cints()
{
    timer_on("CD (MN|LS)");

    // Variables
    int dim1 = nso_;
    int dim2 = nso_;
    int naux = nQ;

    SharedTensor1d D, D2, L1, U1, R1;
    SharedTensor2d bQ, R, L_, L2, Dmat, U, LU, J, K;
    SharedTensor1i pair_to_idx1, pair_to_idx2;

    // Title
    outfile->Printf("\n\tGenerating CD factors ...\n");
    outfile->Printf("\tCD decomposition threshold: %8.2le\n", tol_ldl);

    // SO basis
    bQ = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso_, nso_));
    bQ->read(psio_, PSIF_DFOCC_INTS, true, true);

    // Initial dimension
    size_t n = dim1 * dim2;
    size_t Q = 0;

    outfile->Printf("\tNumber of complete CD factors:   %5li\n",n);

    // Pair to idx mapping
    pair_to_idx1 = SharedTensor1i(new Tensor1i("AB -> A", n));
    pair_to_idx2 = SharedTensor1i(new Tensor1i("AB -> B", n));
    #pragma omp parallel for
    for(int a = 0 ; a < dim1; ++a){
        for(int b = 0 ; b < dim2; ++b){
            int ab = (a*dim2) + b;
	    pair_to_idx1->set(ab,a);
	    pair_to_idx2->set(ab,b);
	}
    }

    // Memory constrasize_t on rows
    size_t max_size_t = std::numeric_limits<int>::max();

    ULI max_rows_ULI = ((memory - n) / (2L * n));
    size_t max_rows = (max_rows_ULI > max_size_t ? max_size_t : max_rows_ULI);

    // Get the diagonal (Q|Q)^(0)
    double* diag = new double[n];
    #pragma omp parallel for
    for(int a = 0 ; a < dim1; ++a){
        for(int b = 0 ; b < dim2; ++b){
            int ab = (a*dim2) + b;
	    double sum = 0.0;
            for(int P = 0 ; P < naux; ++P){
		sum += bQ->get(P,ab) * bQ->get(P,ab);
	    }
	    diag[ab] = sum;
	}
    }

    // Temporary cholesky factor
    std::vector<double*> L;

    // List of selected pivots
    std::vector<int> pivots;

    // Cholesky procedure
    while (Q < n) {

        // Select the pivot
        size_t pivot = 0;
        double dmax = diag[0];
        for (size_t P = 0; P < n; P++) {
            if (dmax < diag[P]) {
                dmax = diag[P];
                pivot = P;
            }
        }

        // Check to see if convergence reached
        if (dmax < tol_ldl || dmax < 0.0) break;

	// Print
        //outfile->Printf("\t%3d     %12.8f\n",Q,dmax);

        // If here, we're trying to add this row
        pivots.push_back(pivot);
        double L_QQ = sqrt(dmax);

        // Check to see if memory constraints are OK
        if (Q > max_rows) {
            throw PSIEXCEPTION("Cholesky: Memory constraints exceeded.");
        }

        // If here, we're really going to add this row
        L.push_back(new double[n]);

        // Compute (ab|Q)
        int c0 = pair_to_idx1->get(pivot);
        int d0 = pair_to_idx2->get(pivot);
        #pragma omp parallel for
        for(size_t i = 0; i < n; i++) {
            int a0 = pair_to_idx1->get(i);
            int b0 = pair_to_idx2->get(i);
            int ab = (a0*dim2) + b0;
            int cd = (c0*dim2) + d0;
	    double sum = 0.0;
            for(int P = 0 ; P < naux; ++P){
                sum += bQ->get(P,ab) * bQ->get(P,cd);
	    }
	    L[Q][i] = sum;
        }

        // [(ab|Q) - L_ab^P L_Q^P]
        for (size_t P = 0; P < Q; P++) {
            C_DAXPY(n,-L[P][pivots[Q]],L[P],1,L[Q],1);
        }

        // 1/L_QQ [(ab|Q) - L_ab^P L_Q^P]
        C_DSCAL(n, 1.0 / L_QQ, L[Q], 1);

        // Zero the upper triangle
        for (size_t P = 0; P < pivots.size(); P++) {
            L[Q][pivots[P]] = 0.0;
        }

        // Set the pivot factor
        L[Q][pivot] = L_QQ;

        // Update the Schur complement diagonal
        for (size_t P = 0; P < n; P++) {
            diag[P] -= L[Q][P] * L[Q][P];
        }

        // Force truly zero elements to zero
        for (size_t P = 0; P < pivots.size(); P++) {
            diag[pivots[P]] = 0.0;
        }

        Q++;
    }
    nQ_cd = static_cast<int>(Q);
    int n_ = static_cast<int>(n);
    outfile->Printf("\tIteratons were done.\n");
    outfile->Printf("\tNumber of computed CD factors:   %5li\n",nQ_cd);

    // Form L
    U = SharedTensor2d(new Tensor2d("L <Q|AB>", nQ_cd, n_));
    #pragma omp parallel for
    for(size_t P = 0; P < Q; P++) {
        int PP = static_cast<int>(P);
        for(size_t i = 0; i < n; i++) {
            int ii = static_cast<int>(i);
	    U->set(PP,ii,L[P][i]);
	}
    }

    // Write L matrix
    U->write(psio_, PSIF_DFOCC_INTS);
    U.reset();

    /*
    // Verification
    LU = SharedTensor2d(new Tensor2d("LU", n_, n_));
    LU->gemm(true,false,U,U,1.0,0.0);
    LU->print();
    LU.reset();

    // Form exact (AB|CD)
    J = SharedTensor2d(new Tensor2d("J (AC|BD)", dim1, dim2, dim1, dim2));
    J->gemm(true, false, bQ, bQ, 1.0, 0.0);
    J->print();
    J.reset();
    */

    // reset
    bQ.reset();
    pair_to_idx1.reset();
    pair_to_idx2.reset();

    timer_off("CD (MN|LS)");

} // end cd_aob_cints

//=======================================================
//       Form CD (AB|CD) ints
//=======================================================
void DFOCC::cd_abcd_cints()
{
    timer_on("CD (AB|CD)");

    // Variables
    int dim1 = navirA;
    int dim2 = navirA;
    int naux = nQ;

    SharedTensor1d D, D2, L1, U1, R1;
    SharedTensor2d bQ, R, L_, L2, Dmat, U, LU, J, K;
    SharedTensor1i pair_to_idx1, pair_to_idx2;

    // Title
    outfile->Printf("\n\tGenerating CD factors ...\n");
    outfile->Printf("\tCD decomposition threshold: %8.2le\n", tol_ldl);

    // SO basis
    //bQ = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    //bQ->read(psio_, PSIF_DFOCC_INTS, true, true);
    bQ = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, ntri_abAA));
    bQ->read(psio_, PSIF_DFOCC_INTS);

    // Initial dimension
    //size_t n = dim1 * dim2;
    size_t n = ntri_abAA;
    size_t Q = 0;

    outfile->Printf("\tNumber of complete CD factors:   %5li\n",n);

    /*
    // Pair to idx mapping
    pair_to_idx1 = SharedTensor1i(new Tensor1i("AB -> A", n));
    pair_to_idx2 = SharedTensor1i(new Tensor1i("AB -> B", n));
    #pragma omp parallel for
    for(int a = 0 ; a < dim1; ++a){
        for(int b = 0 ; b < dim2; ++b){
            int ab = (a*dim2) + b;
	    pair_to_idx1->set(ab,a);
	    pair_to_idx2->set(ab,b);
	}
    }
    */

    // Memory constrasize_t on rows
    size_t max_size_t = std::numeric_limits<int>::max();

    ULI max_rows_ULI = ((memory - n) / (2L * n));
    size_t max_rows = (max_rows_ULI > max_size_t ? max_size_t : max_rows_ULI);

    // Get the diagonal (Q|Q)^(0)
    double* diag = new double[n];
    /*
    #pragma omp parallel for
    for(int a = 0 ; a < dim1; ++a){
        for(int b = 0 ; b < dim2; ++b){
            int ab = (a*dim2) + b;
	    double sum = 0.0;
            for(int P = 0 ; P < naux; ++P){
		sum += bQ->get(P,ab) * bQ->get(P,ab);
	    }
	    diag[ab] = sum;
	}
    }
    */
    #pragma omp parallel for
    for(int ab = 0; ab < n; ++ab){
	 double sum = 0.0;
         for(int P = 0 ; P < naux; ++P){
	     sum += bQ->get(P,ab) * bQ->get(P,ab);
	 }
	 diag[ab] = sum;
    }

    // Temporary cholesky factor
    std::vector<double*> L;

    // List of selected pivots
    std::vector<int> pivots;

    // Cholesky procedure
    while (Q < n) {

        // Select the pivot
        size_t pivot = 0;
        double dmax = diag[0];
        for (size_t P = 0; P < n; P++) {
            if (dmax < diag[P]) {
                dmax = diag[P];
                pivot = P;
            }
        }

        // Check to see if convergence reached
        if (dmax < tol_ldl || dmax < 0.0) break;

	// Print
        //outfile->Printf("\t%3d     %12.8f\n",Q,dmax);

        // If here, we're trying to add this row
        pivots.push_back(pivot);
        double L_QQ = sqrt(dmax);

        // Check to see if memory constraints are OK
        if (Q > max_rows) {
            throw PSIEXCEPTION("Cholesky: Memory constraints exceeded.");
        }

        // If here, we're really going to add this row
        L.push_back(new double[n]);

        // Compute (ab|Q)
	/*
        int c0 = pair_to_idx1->get(pivot);
        int d0 = pair_to_idx2->get(pivot);
        #pragma omp parallel for
        for(size_t i = 0; i < n; i++) {
            int a0 = pair_to_idx1->get(i);
            int b0 = pair_to_idx2->get(i);
            int ab = (a0*dim2) + b0;
            int cd = (c0*dim2) + d0;
	    double sum = 0.0;
            for(int P = 0 ; P < naux; ++P){
                sum += bQ->get(P,ab) * bQ->get(P,cd);
	    }
	    L[Q][i] = sum;
        }
	*/
        #pragma omp parallel for
        for(size_t i = 0; i < n; i++) {
	    double sum = 0.0;
            for(int P = 0 ; P < naux; ++P){
                sum += bQ->get(P,i) * bQ->get(P,pivot);
	    }
	    L[Q][i] = sum;
        }

        // [(ab|Q) - L_ab^P L_Q^P]
        for (size_t P = 0; P < Q; P++) {
            C_DAXPY(n,-L[P][pivots[Q]],L[P],1,L[Q],1);
        }

        // 1/L_QQ [(ab|Q) - L_ab^P L_Q^P]
        C_DSCAL(n, 1.0 / L_QQ, L[Q], 1);

        // Zero the upper triangle
        for (size_t P = 0; P < pivots.size(); P++) {
            L[Q][pivots[P]] = 0.0;
        }

        // Set the pivot factor
        L[Q][pivot] = L_QQ;

        // Update the Schur complement diagonal
        for (size_t P = 0; P < n; P++) {
            diag[P] -= L[Q][P] * L[Q][P];
        }

        // Force truly zero elements to zero
        for (size_t P = 0; P < pivots.size(); P++) {
            diag[pivots[P]] = 0.0;
        }

        Q++;
    }
    nQ_cd = static_cast<int>(Q);
    int n_ = static_cast<int>(n);
    outfile->Printf("\tIteratons were done.\n");
    outfile->Printf("\tNumber of computed CD factors:   %5li\n",nQ_cd);

    // Form L
    U = SharedTensor2d(new Tensor2d("L <Q|AB>", nQ_cd, n_));
    #pragma omp parallel for
    for(size_t P = 0; P < Q; P++) {
        int PP = static_cast<int>(P);
        for(size_t i = 0; i < n; i++) {
            int ii = static_cast<int>(i);
	    U->set(PP,ii,L[P][i]);
	}
    }

    // Write L matrix
    U->write(psio_, PSIF_DFOCC_INTS);
    U.reset();

    /*
    // Verification
    LU = SharedTensor2d(new Tensor2d("LU", n_, n_));
    LU->gemm(true,false,U,U,1.0,0.0);
    LU->print();
    LU.reset();

    // Form exact (AB|CD)
    //J = SharedTensor2d(new Tensor2d("J (AB|CD)", dim1, dim2, dim1, dim2));
    J = SharedTensor2d(new Tensor2d("J (A>=B|C>=D)", n_, n_));
    J->gemm(true, false, bQ, bQ, 1.0, 0.0);
    J->print();
    J.reset();
    */

    // reset
    bQ.reset();
    //pair_to_idx1.reset();
    //pair_to_idx2.reset();

    timer_off("CD (AB|CD)");

} // end cd_abcd_cints

//=======================================================
//       Form CD <AB|CD> ints
//=======================================================
void DFOCC::cd_abcd_xints()
{
    timer_on("CD <AB|CD>");

    // Variables
    int dim1 = navirA;
    int dim2 = navirA;
    int naux = nQ;

    SharedTensor1d D, D2, L1, U1, R1;
    SharedTensor2d bQ, R, L_, L2, Dmat, U, LU, J, K;
    SharedTensor1i pair_to_idx1, pair_to_idx2;

    // Title
    outfile->Printf("\n\tGenerating CD factors ...\n");
    outfile->Printf("\tCD decomposition threshold: %8.2le\n", tol_ldl);

    // SO basis
    bQ = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQ->read(psio_, PSIF_DFOCC_INTS, true, true);

    // Initial dimension
    size_t n = dim1 * dim2;
    size_t Q = 0;

    outfile->Printf("\tNumber of complete CD factors:   %5li\n",n);

    // Pair to idx mapping
    pair_to_idx1 = SharedTensor1i(new Tensor1i("AB -> A", n));
    pair_to_idx2 = SharedTensor1i(new Tensor1i("AB -> B", n));
    #pragma omp parallel for
    for(int a = 0 ; a < dim1; ++a){
        for(int b = 0 ; b < dim2; ++b){
            int ab = (a*dim2) + b;
	    pair_to_idx1->set(ab,a);
	    pair_to_idx2->set(ab,b);
	}
    }

    // Memory constrasize_t on rows
    size_t max_size_t = std::numeric_limits<int>::max();

    ULI max_rows_ULI = ((memory - n) / (2L * n));
    size_t max_rows = (max_rows_ULI > max_size_t ? max_size_t : max_rows_ULI);

    // Get the diagonal (Q|Q)^(0)
    double* diag = new double[n];
    #pragma omp parallel for
    for(int a = 0 ; a < dim1; ++a){
        int aa = (a*dim2) + a;
        for(int b = 0 ; b < dim2; ++b){
            int bb = (b*dim2) + b;
	    double sum = 0.0;
            for(int P = 0 ; P < naux; ++P){
		sum += bQ->get(P,aa) * bQ->get(P,bb);
	    }
            int ab = (a*dim2) + b;
	    diag[ab] = sum;
	}
    }

    // Temporary cholesky factor
    std::vector<double*> L;

    // List of selected pivots
    std::vector<int> pivots;

    // Cholesky procedure
    while (Q < n) {

        // Select the pivot
        size_t pivot = 0;
        double dmax = diag[0];
        for (size_t P = 0; P < n; P++) {
            if (dmax < diag[P]) {
                dmax = diag[P];
                pivot = P;
            }
        }

        // Check to see if convergence reached
        if (dmax < tol_ldl || dmax < 0.0) break;

	// Print
        //outfile->Printf("\t%3d     %12.8f\n",Q,dmax);

        // If here, we're trying to add this row
        pivots.push_back(pivot);
        double L_QQ = sqrt(dmax);

        // Check to see if memory constraints are OK
        if (Q > max_rows) {
            throw PSIEXCEPTION("Cholesky: Memory constraints exceeded.");
        }

        // If here, we're really going to add this row
        L.push_back(new double[n]);

        // Compute (ab|Q)
        int c0 = pair_to_idx1->get(pivot);
        int d0 = pair_to_idx2->get(pivot);
        #pragma omp parallel for
        for(size_t i = 0; i < n; i++) {
            int a0 = pair_to_idx1->get(i);
            int b0 = pair_to_idx2->get(i);
            int ac = (a0*dim2) + c0;
            int bd = (b0*dim2) + d0;
	    double sum = 0.0;
            for(int P = 0 ; P < naux; ++P){
                sum += bQ->get(P,ac) * bQ->get(P,bd);
	    }
	    L[Q][i] = sum;
        }

        // [(ab|Q) - L_ab^P L_Q^P]
        for (size_t P = 0; P < Q; P++) {
            C_DAXPY(n,-L[P][pivots[Q]],L[P],1,L[Q],1);
        }

        // 1/L_QQ [(ab|Q) - L_ab^P L_Q^P]
        C_DSCAL(n, 1.0 / L_QQ, L[Q], 1);

        // Zero the upper triangle
        for (size_t P = 0; P < pivots.size(); P++) {
            L[Q][pivots[P]] = 0.0;
        }

        // Set the pivot factor
        L[Q][pivot] = L_QQ;

        // Update the Schur complement diagonal
        for (size_t P = 0; P < n; P++) {
            diag[P] -= L[Q][P] * L[Q][P];
        }

        // Force truly zero elements to zero
        for (size_t P = 0; P < pivots.size(); P++) {
            diag[pivots[P]] = 0.0;
        }

        Q++;
    }
    nQ_cd = static_cast<int>(Q);
    int n_ = static_cast<int>(n);
    outfile->Printf("\tIteratons were done.\n");
    outfile->Printf("\tNumber of computed CD factors:   %5li\n",nQ_cd);

    // Form L
    U = SharedTensor2d(new Tensor2d("L <Q|AB>", nQ_cd, n_));
    #pragma omp parallel for
    for(size_t P = 0; P < Q; P++) {
        int PP = static_cast<int>(P);
        for(size_t i = 0; i < n; i++) {
            int ii = static_cast<int>(i);
	    U->set(PP,ii,L[P][i]);
	}
    }

    // Write L matrix
    U->write(psio_, PSIF_DFOCC_INTS);
    U.reset();

    /*
    // Verification
    LU = SharedTensor2d(new Tensor2d("LU", n_, n_));
    LU->gemm(true,false,U,U,1.0,0.0);
    LU->print();
    LU.reset();

    // Form exact <AB|CD>
    J = SharedTensor2d(new Tensor2d("J (AC|BD)", dim1, dim2, dim1, dim2));
    J->gemm(true, false, bQ, bQ, 1.0, 0.0);
    K = SharedTensor2d(new Tensor2d("K <AB|CD>", navirA, navirA, navirA, navirA));
    K->sort(1324, J, 1.0, 0.0);
    J.reset();
    K->print();
    K.reset();
    */

    // reset
    bQ.reset();
    pair_to_idx1.reset();
    pair_to_idx2.reset();

    timer_off("CD <AB|CD>");

} // end cd_abcd_xints





}} // Namespaces






