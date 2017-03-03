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
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/sieve.h"
#include "psi4/psifiles.h"
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

void DFOCC::trans_ref()
{
    // Read SO integrals
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", nQ_ref, nso_, nso_));
    bQso->read(psio_, PSIF_DFOCC_INTS, true, true);

    // Form B(Q,ij)
    timer_on("Form B(Q,ij)");
    b_oo_ref();
    timer_off("Form B(Q,ij)");

    // Form B(Q,ia)
    timer_on("Form B(Q,ia)");
    b_ov_ref();
    timer_off("Form B(Q,ia)");

    // Form B(Q,ab)
    timer_on("Form B(Q,ab)");
    b_vv_ref();
    timer_off("Form B(Q,ab)");
    bQso.reset();

    /*
    if (time4grad == 1) {
        cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", nQ_ref, nso2_));
        c_oo_ref();
        cQso.reset();
    }
    */
}

//=======================================================
//          DF SCF
//=======================================================
void DFOCC::df_ref()
{
    //outfile->Printf("\tComputing DF-BASIS-SCF integrals... \n");

  //if (read_scf_3index == "TRUE" && dertype == "NONE") {
  // 1.  read scf 3-index integrals from disk

  // get ntri from sieve
  std::shared_ptr<ERISieve> sieve (new ERISieve(basisset_, options_.get_double("INTS_TOLERANCE")));
  const std::vector<std::pair<int, int> >& function_pairs = sieve->function_pairs();
  long int ntri_cd = function_pairs.size();

      // read integrals from disk if they were generated in the SCF
      if ( options_.get_str("SCF_TYPE") == "DF") {
          outfile->Printf("\tReading DF integrals from disk ...\n");
          std::shared_ptr<BasisSet> primary = get_basisset("ORBITAL");
          std::shared_ptr<BasisSet> auxiliary = get_basisset("DF_BASIS_SCF");
          std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
          nQ_ref = auxiliary->nbf();

          // ntri comes from sieve above
          std::shared_ptr<Matrix> Qmn = SharedMatrix(new Matrix("Qmn Integrals",nQ_ref,ntri_cd));
          double** Qmnp = Qmn->pointer();
          psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
          psio_->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri_cd * nQ_ref);
          psio_->close(PSIF_DFSCF_BJ,1);

          bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", nQ_ref, nso_, nso_));
          for (long int mn = 0; mn < ntri_cd; mn++) {
              long int m = function_pairs[mn].first;
              long int n = function_pairs[mn].second;
              for (long int P = 0; P < nQ_ref; P++) {
                  bQso->set(P, (m*nso_) + n, Qmnp[P][mn]);
                  bQso->set(P, (n*nso_) + m, Qmnp[P][mn]);
              }
          }
          bQso->write(psio_, PSIF_DFOCC_INTS, true, true);

          if (dertype == "FIRST") {
              // Form J^-1/2
              timer_on("Form J");
              formJ_ref(auxiliary, zero);
              timer_off("Form J");
          }// end if (dertype == "FIRST")

      }// end if ( options_.get_str("SCF_TYPE") == "DF" )

      // read integrals from disk if they were generated in the SCF
      else if ( options_.get_str("SCF_TYPE") == "CD") {
          outfile->Printf("\tReading Cholesky vectors from disk ...\n");
          //nQ_ref = Process::environment.globals["NAUX (SCF)"];
          outfile->Printf("\tCholesky decomposition threshold: %8.2le\n", options_.get_double("CHOLESKY_TOLERANCE"));
          //outfile->Printf("\tNumber of Cholesky vectors:   %5li\n",nQ_ref);

          // ntri comes from sieve above
          psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
          psio_->read_entry(PSIF_DFSCF_BJ, "length", (char*)&nQ_ref, sizeof(long int));
          std::shared_ptr<Matrix> Qmn = SharedMatrix(new Matrix("Qmn Integrals",nQ_ref,ntri_cd));
          double** Qmnp = Qmn->pointer();
          psio_->read_entry(PSIF_DFSCF_BJ, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri_cd * nQ_ref);
          psio_->close(PSIF_DFSCF_BJ,1);

          bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", nQ_ref, nso_, nso_));
          for (long int mn = 0; mn < ntri_cd; mn++) {
              long int m = function_pairs[mn].first;
              long int n = function_pairs[mn].second;
              for (long int P = 0; P < nQ_ref; P++) {
                  bQso->set(P, (m*nso_) + n, Qmnp[P][mn]);
                  bQso->set(P, (n*nso_) + m, Qmnp[P][mn]);
              }
          }
          bQso->write(psio_, PSIF_DFOCC_INTS, true, true);
      }// end else if ( options_.get_str("SCF_TYPE") == "CD" )

      //else throw PSIEXCEPTION("SCF_TYPE should be DF or CD");
  //}// end if (read_scf_3index == "TRUE")


  //else if (read_scf_3index == "FALSE") {
  else {
    // Read in the basis set informations
    std::shared_ptr<BasisSet> primary_ = get_basisset("ORBITAL");
    std::shared_ptr<BasisSet> auxiliary_ = get_basisset("DF_BASIS_SCF");
    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());

    // Read number of auxilary basis
    nQ_ref = auxiliary_->nbf();

    // Form J^-1/2
    timer_on("Form J");
    formJ_ref(auxiliary_, zero);
    timer_off("Form J");

    // Form B(Q,mu nu)
    timer_on("Form B(Q,munu)");
    b_so_ref(primary_, auxiliary_, zero);
    timer_off("Form B(Q,munu)");
  }// end if (read_scf_3index == "FALSE")

  //outfile->Printf("\tDF-BASIS-SCF integrals were done. \n");
} // end df_ref


//=======================================================
//          form J(P,Q)^-1/2
//=======================================================
void DFOCC::formJ_ref(std::shared_ptr<BasisSet> auxiliary_, std::shared_ptr<BasisSet> zero)
{
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    double **J = block_matrix(nQ_ref, nQ_ref);
    J_mhalf = block_matrix(nQ_ref, nQ_ref);

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

    double **J = block_matrix(nQ_ref, nQ_ref);
    J_mhalf = block_matrix(nQ_ref, nQ_ref);
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
    int lwork = nQ_ref * 3;
    double* eigval = init_array(nQ_ref);
    double* work = init_array(lwork);
    int status = C_DSYEV('v', 'u', nQ_ref, J[0], nQ_ref, eigval, work, lwork);
    if(status){
        throw PsiException("Diagonalization of J failed", __FILE__, __LINE__);
    }
    free(work);

    // Now J contains the eigenvectors of the original J
    // Copy J to J_copy
    double **J_copy = block_matrix(nQ_ref, nQ_ref);
    C_DCOPY(nQ_ref*nQ_ref, J[0], 1, J_copy[0], 1);

    // Now form J^{-1/2} = U(T)*j^{-1/2}*U,
    // where j^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of J
    for(int i=0; i<nQ_ref; ++i){
        eigval[i] = (eigval[i] < 1.0E-10) ? 0.0 : 1.0 / sqrt(eigval[i]);
        // scale one set of eigenvectors by the diagonal elements j^{-1/2}
        C_DSCAL(nQ_ref, eigval[i], J[i], 1);
    }
    free(eigval);

    // J_mhalf = J_copy(T) * J
    C_DGEMM('t','n', nQ_ref, nQ_ref, nQ_ref, 1.0, J_copy[0], nQ_ref, J[0], nQ_ref, 0.0, J_mhalf[0], nQ_ref);
    free_block(J);
    free_block(J_copy);

    // write J
    Jmhalf = SharedTensor2d(new Tensor2d("DF_BASIS_SCF Jmhalf <P|Q>", nQ_ref, nQ_ref));
    Jmhalf->set(J_mhalf);
    Jmhalf->write(psio_, PSIF_DFOCC_INTS);
    Jmhalf.reset();

} // end formJ

//=======================================================
//          form b(Q, mu nu)
//=======================================================
void DFOCC::b_so_ref(std::shared_ptr<BasisSet> primary_, std::shared_ptr<BasisSet> auxiliary_, std::shared_ptr<BasisSet> zero)
{
    bQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mn)", nQ_ref, nso_, nso_));
    double** Ap = block_matrix(nQ_ref, nso2_);
    double** Bp = block_matrix(nQ_ref, nso2_);

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = Process::environment.get_n_threads();
    #endif

    std::shared_ptr<ERISieve> sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, 0.0));
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
        int pstop  = (Pstop == auxiliary_->nshell() ? nQ_ref : auxiliary_->shell(Pstop ).function_index());
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

    C_DGEMM('N','N', nQ_ref, nso2_, nQ_ref, 1.0, J_mhalf[0], nQ_ref, Bp[0], nso2_, 0.0, Ap[0], nso2_);
    bQso->set(Ap);
    bQso->write(psio_, PSIF_DFOCC_INTS, true, true);
    if (print_ > 3) bQso->print();
    free_block(Bp);
    free_block(J_mhalf);
    free_block(Ap);
    bQso.reset();

    /*
    // Build C(Q, mu nu)
    double** Cp = block_matrix(nQ_ref, nso2_);
    C_DGEMM('N','N', nQ_ref, nso2_, nQ_ref, 1.0, J_mhalf[0], nQ_ref, Ap[0], nso2_, 0.0, Cp[0], nso2_);
    free_block(J_mhalf);
    free_block(Ap);
    cQso = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mn)", nQ_ref, nso2_));
    cQso->set(Cp);
    free_block(Cp);
    cQso->write(psio_, PSIF_DFOCC_INTS);
    cQso.reset();
    */

} // end b_so

//=======================================================
//          form b(Q,ij) : all
//=======================================================
void DFOCC::b_oo_ref()
{
    bQnoA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mO)", nQ_ref, nso_ * noccA));
    bQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OO)", nQ_ref, noccA * noccA));
    bQnoA->contract(false, false, nQ_ref * nso_, noccA, nso_, bQso, CoccA, 1.0, 0.0);
    bQooA->contract233(true, false, noccA, noccA, CoccA, bQnoA, 1.0, 0.0);
    bQnoA.reset();
    bQooA->write(psio_, PSIF_DFOCC_INTS);
    bQooA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQnoB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mo)", nQ_ref, nso_ * noccB));
    bQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|oo)", nQ_ref, noccB * noccB));
    bQnoB->contract(false, false, nQ_ref * nso_, noccB, nso_, bQso, CoccB, 1.0, 0.0);
    bQooB->contract233(true, false, noccB, noccB, CoccB, bQnoB, 1.0, 0.0);
    bQnoB.reset();
    bQooB->write(psio_, PSIF_DFOCC_INTS);
    bQooB.reset();
 }
} // end b_oo

//=======================================================
//          form b(Q,ia) : all
//=======================================================
void DFOCC::b_ov_ref()
{
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mV)", nQ_ref, nso_ * nvirA));
    bQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|OV)", nQ_ref, noccA * nvirA));
    bQnvA->contract(false, false, nQ_ref * nso_, nvirA, nso_, bQso, CvirA, 1.0, 0.0);
    bQovA->contract233(true, false, noccA, nvirA, CoccA, bQnvA, 1.0, 0.0);
    bQovA->write(psio_, PSIF_DFOCC_INTS);
    bQnvA->write(psio_, PSIF_DFOCC_INTS);
    bQovA.reset();
    bQnvA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mv)", nQ_ref, nso_ * nvirB));
    bQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|ov)", nQ_ref, noccB * nvirB));
    bQnvB->contract(false, false, nQ_ref * nso_, nvirB, nso_, bQso, CvirB, 1.0, 0.0);
    bQovB->contract233(true, false, noccB, nvirB, CoccB, bQnvB, 1.0, 0.0);
    bQovB->write(psio_, PSIF_DFOCC_INTS);
    bQnvB->write(psio_, PSIF_DFOCC_INTS);
    bQovB.reset();
    bQnvB.reset();
 }
} // end b_ov

//=======================================================
//          form b(Q,ab) : all
//=======================================================
void DFOCC::b_vv_ref()
{
    bQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|VV)", nQ_ref, nvirA, nvirA));
    bQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mV)", nQ_ref, nso_ * nvirA));
    bQnvA->read(psio_, PSIF_DFOCC_INTS);
    bQvvA->contract233(true, false, nvirA, nvirA, CvirA, bQnvA, 1.0, 0.0);
    bQnvA.reset();
    bQvvA->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQvvA.reset();

 if (reference_ == "UNRESTRICTED") {
    bQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|vv)", nQ_ref, nvirB, nvirB));
    bQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF B (Q|mv)", nQ_ref, nso_ * nvirB));
    bQnvB->read(psio_, PSIF_DFOCC_INTS);
    bQvvB->contract233(true, false, nvirB, nvirB, CvirB, bQnvB, 1.0, 0.0);
    bQnvB.reset();
    bQvvB->write(psio_, PSIF_DFOCC_INTS, true, true);
    bQvvB.reset();
 }
} // end b_vv

//=======================================================
//          form c(Q,ij) : all
//=======================================================
void DFOCC::c_oo_ref()
{
    cQnoA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mO)", nQ_ref, nso_ * noccA));
    cQooA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|OO)", nQ_ref, noccA * noccA));
    cQnoA->contract(false, false, nQ_ref * nso_, noccA, nso_, cQso, CoccA, 1.0, 0.0);
    cQooA->contract233(true, false, noccA, noccA, CoccA, cQnoA, 1.0, 0.0);
    cQnoA.reset();
    cQooA->write(psio_, PSIF_DFOCC_INTS);
    cQooA.reset();

 if (reference_ == "UNRESTRICTED") {
    cQnoB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mo)", nQ_ref, nso_ * noccB));
    cQooB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|oo)", nQ_ref, noccB * noccB));
    cQnoB->contract(false, false, nQ_ref * nso_, noccB, nso_, cQso, CoccB, 1.0, 0.0);
    cQooB->contract233(true, false, noccB, noccB, CoccB, cQnoB, 1.0, 0.0);
    cQnoB.reset();
    cQooB->write(psio_, PSIF_DFOCC_INTS);
    cQooB.reset();
 }
} // end c_oo

//=======================================================
//          form c(Q,ia) : all
//=======================================================
void DFOCC::c_ov_ref()
{
    cQnvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mV)", nQ_ref, nso_ * nvirA));
    cQovA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|OV)", nQ_ref, noccA * nvirA));
    cQnvA->contract(false, false, nQ_ref * nso_, nvirA, nso_, cQso, CvirA, 1.0, 0.0);
    cQovA->contract233(true, false, noccA, nvirA, CoccA, cQnvA, 1.0, 0.0);
    //if (trans_ab == 0) cQnvA.reset();
    cQovA->write(psio_, PSIF_DFOCC_INTS);
    cQovA.reset();

 if (reference_ == "UNRESTRICTED") {
    cQnvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|mv)", nQ_ref, nso_ * nvirB));
    cQovB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|ov)", nQ_ref, noccB * nvirB));
    cQnvB->contract(false, false, nQ_ref * nso_, nvirB, nso_, cQso, CvirB, 1.0, 0.0);
    cQovB->contract233(true, false, noccB, nvirB, CoccB, cQnvB, 1.0, 0.0);
    //if (trans_ab == 0) cQnvB.reset();
    cQovB->write(psio_, PSIF_DFOCC_INTS);
    cQovB.reset();
 }
} // end c_ov

//=======================================================
//          form c(Q,ab) : all
//=======================================================
void DFOCC::c_vv_ref()
{
    cQvvA = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|VV)", nQ_ref, nvirA * nvirA));
    cQvvA->contract233(true, false, nvirA, nvirA, CvirA, cQnvA, 1.0, 0.0);
    cQnvA.reset();
    cQvvA->write(psio_, PSIF_DFOCC_INTS);
    cQvvA.reset();

 if (reference_ == "UNRESTRICTED") {
    cQvvB = SharedTensor2d(new Tensor2d("DF_BASIS_SCF C (Q|vv)", nQ_ref, nvirB * nvirB));
    cQvvB->contract233(true, false, nvirB, nvirB, CvirB, cQnvB, 1.0, 0.0);
    cQnvB.reset();
    cQvvB->write(psio_, PSIF_DFOCC_INTS);
    cQvvB.reset();
 }
} // end c_vv

}} // Namespaces

