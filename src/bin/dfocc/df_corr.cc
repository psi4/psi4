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

#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libmints/mints.h>
#include <libmints/sieve.h>
#include "psi4-dec.h"

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
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso2_));
    bQso->read(psio_, PSIF_DFOCC_INTS);

    trans_ab = 1;
    if (orb_opt_ == "TRUE" || dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
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
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso2_));
    bQso->read(psio_, PSIF_DFOCC_INTS);

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
    //fprintf(outfile,"\tComputing DF-BASIS-CC integrals... \n"); fflush(outfile);

    // Read in the basis set informations
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> auxiliary_ = BasisSet::construct(parser, reference_wavefunction_->molecule(), "DF_BASIS_CC");
    boost::shared_ptr<BasisSet> primary_ = BasisSet::construct(parser, reference_wavefunction_->molecule(), "BASIS");
    boost::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
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
void DFOCC::formJ(boost::shared_ptr<BasisSet> auxiliary_, boost::shared_ptr<BasisSet> zero)
{

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    double **J = block_matrix(nQ, nQ);
    J_mhalf = block_matrix(nQ, nQ);

    // => Integrals <= //
    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_,zero,auxiliary_,zero));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > Jint;
    std::vector<const double*> buffer;
    for (int t = 0; t < nthreads; t++) {
        Jint.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri()));
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
    boost::shared_ptr<IntegralFactory> rifactory_J(new IntegralFactory(auxiliary_, zero, auxiliary_, zero));
    boost::shared_ptr<TwoBodyAOInt> Jint(rifactory_J->eri());

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
void DFOCC::b_so(boost::shared_ptr<BasisSet> primary_, boost::shared_ptr<BasisSet> auxiliary_, boost::shared_ptr<BasisSet> zero)
{
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mn)", nQ, nso2_));
    double** Ap = block_matrix(nQ, nso2_); 
    double** Bp = block_matrix(nQ, nso2_); 

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    boost::shared_ptr<ERISieve> sieve_ = boost::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff));
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
    boost::shared_ptr<IntegralFactory> rifactory2(new IntegralFactory(auxiliary_, zero, primary_, primary_));
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
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    boost::shared_ptr<TwoBodyAOInt> eri(fact->eri());
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
    bQso->write(psio_, PSIF_DFOCC_INTS);
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
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA * navirA));
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mA)", nQ, nso_ * navirA));
    bQnvA->read(psio_, PSIF_DFOCC_INTS);
    bQabA->contract233(true, false, navirA, navirA, CavirA, bQnvA, 1.0, 0.0);
    bQnvA.reset();
    bQabA->write(psio_, PSIF_DFOCC_INTS);
    bQabA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB * navirB));
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ma)", nQ, nso_ * navirB));
    bQnvB->read(psio_, PSIF_DFOCC_INTS);
    bQabB->contract233(true, false, navirB, navirB, CavirB, bQnvB, 1.0, 0.0);
    bQnvB.reset();
    bQabB->write(psio_, PSIF_DFOCC_INTS);
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
    bQvvA->write(psio_, PSIF_DFOCC_INTS);

    // Form active b(Q,ab)
    bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
    bQabA->form_b_ab(bQvvA);
    bQvvA.reset();
    bQabA->write(psio_, PSIF_DFOCC_INTS);
    bQabA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|vv)", nQ, nvirB, nvirB));
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|mv)", nQ, nso_ * nvirB));
    bQnvB->read(psio_, PSIF_DFOCC_INTS);
    bQvvB->contract233(true, false, nvirB, nvirB, CvirB, bQnvB, 1.0, 0.0);
    bQnvB.reset();
    bQvvB->write(psio_, PSIF_DFOCC_INTS);

    // Form active b(Q,ab)
    bQabB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|ab)", nQ, navirB, navirB));
    bQabB->form_b_ab(bQvvB);
    bQvvB.reset();
    bQabB->write(psio_, PSIF_DFOCC_INTS);
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
    HovA = HvoA->transpose();
    HvvA->form_vv(noccA, HmoA);

 if (reference_ == "UNRESTRICTED") {
    HmoB->transform(Hso, CmoB);
    if (print_ > 2) HmoB->print();
    // Blocks
    HooB->form_oo(HmoB);
    HvoB->form_vo(HmoB);
    HovB = HvoB->transpose();
    HvvB->form_vv(noccB, HmoB);
 } // uhf

} // end trans_oei

}} // Namespaces



