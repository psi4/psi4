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

#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>
#include <libqt/qt.h>
#include <math.h>
#include <limits>
#include <vector>
#include "cholesky.h"
#include <psifiles.h>
#include <psi4-dec.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.h>
#ifdef _OPENMP
    #include <omp.h>
#endif

namespace psi {

Cholesky::Cholesky(double delta, unsigned long int memory)
    : delta_(delta), memory_(memory), Q_(0)
{
}
Cholesky::~Cholesky()
{
}
void Cholesky::choleskify()
{
    if(read_previous_cholesky_vector_)
    {
        boost::shared_ptr<PSIO> psio (new PSIO());
        psio_address addr = PSIO_ZERO;
        int file_unit = PSIF_DFSCF_BJ;
        psio->open(file_unit, PSIO_OPEN_OLD);
        psio->read_entry(file_unit, "length", (char*) &Q_, sizeof(long int));

        size_t n = N();
        L_ = SharedMatrix(new Matrix("Partial Cholesky", Q_, n));

        double** Lp = L_->pointer();
        psio->read_entry(file_unit,"(Q|mn) Integrals",(char*) Lp[0], sizeof(double) * Q_ * n);
        psio->close(file_unit, 1);
    }
    else {

        // Initial dimensions
        size_t n = N();
        Q_ = 0;

        // Memory constrasize_t on rows
        size_t max_size_t = std::numeric_limits<int>::max();

        ULI max_rows_ULI = ((memory_ - n) / (2L * n));
        size_t max_rows = (max_rows_ULI > max_size_t ? max_size_t : max_rows_ULI);

        // Get the diagonal (Q|Q)^(0)
        double* diag = new double[n];
        //outfile->Printf("\n Compute diagonal:");
        Timer diagonal_time;
        compute_diagonal(diag);
        //outfile->Printf(" %8.6f s\n", diagonal_time.get());

        // Temporary cholesky factor
        std::vector<double*> L;

        // List of selected pivots
        std::vector<int> pivots;

        // Cholesky procedure
        Timer cholesky_procedure;
        while (Q_ < n) {

            // Select the pivot
            size_t pivot = 0;
            double Dmax = diag[0];
            for (size_t P = 0; P < n; P++) {
                if (Dmax < diag[P]) {
                    Dmax = diag[P];
                    pivot = P;
                }
            }

            // Check to see if convergence reached
            if (Dmax < delta_ || Dmax < 0.0) break;

            // If here, we're trying to add this row
            pivots.push_back(pivot);
            double L_QQ = sqrt(Dmax);

            // Check to see if memory constraints are OK
            if (Q_ > max_rows) {
            }

            // If here, we're really going to add this row
            L.push_back(new double[n]);

            // (m|Q)
            compute_row(pivot, L[Q_]);

            // [(m|Q) - L_m^P L_Q^P]
            Timer daxpy_time;
            for (size_t P = 0; P < Q_; P++) {
                C_DAXPY(n,-L[P][pivots[Q_]],L[P],1,L[Q_],1);
            }

            // 1/L_QQ [(m|Q) - L_m^P L_Q^P]
            C_DSCAL(n, 1.0 / L_QQ, L[Q_], 1);

            // Zero the upper triangle
            for (size_t P = 0; P < pivots.size(); P++) {
                L[Q_][pivots[P]] = 0.0;
            }

            // Set the pivot factor
            L[Q_][pivot] = L_QQ;

            // Update the Schur complement diagonal
            for (size_t P = 0; P < n; P++) {
                diag[P] -= L[Q_][P] * L[Q_][P];
            }

            // Force truly zero elements to zero
            for (size_t P = 0; P < pivots.size(); P++) {
                diag[pivots[P]] = 0.0;
            }

            Q_++;
        }
        //outfile->Printf("\n Cholesky Procedure takes %8.8f s", cholesky_procedure.get());

        // Copy into a more permanant Matrix object
        L_ = SharedMatrix(new Matrix("Partial Cholesky", Q_, n));
        double** Lp = L_->pointer();

        for (size_t Q = 0; Q < Q_; Q++) {
            ::memcpy(static_cast<void*>(Lp[Q]), static_cast<void*>(L[Q]), n * sizeof(double));
            delete[] L[Q];
        }
    }
}

CholeskyMatrix::CholeskyMatrix(SharedMatrix A, double delta, unsigned long int memory) :
    A_(A), Cholesky(delta, memory)
{
    if (A_->nirrep() != 1)
        throw PSIEXCEPTION("CholeskyMatrix only supports C1 matrices");
    if (A_->rowspi()[0] != A_->colspi()[0])
        throw PSIEXCEPTION("CholeskyMatrix only supports square matrices");
}
CholeskyMatrix::~CholeskyMatrix()
{
}
size_t CholeskyMatrix::N()
{
    return A_->rowspi()[0];
}
void CholeskyMatrix::compute_diagonal(double* target)
{
    size_t n = N();
    double** Ap = A_->pointer();
    for (size_t i = 0; i < n; i++) {
        target[i] = Ap[i][i];
    }
}
void CholeskyMatrix::compute_row(int row, double* target)
{
    ::memcpy(static_cast<void*>(target),static_cast<void*>(A_->pointer()[row]),N() * sizeof(double));
}

CholeskyERI::CholeskyERI(boost::shared_ptr<TwoBodyAOInt> integral, double schwarz,
    double delta, unsigned long int memory) :
    integral_(integral), schwarz_(schwarz), Cholesky(delta, memory)
{
    basisset_ = integral_->basis();
    //{
    //    integral_threads_.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->eri()));
    //}
}
CholeskyERI::~CholeskyERI()
{
}
size_t CholeskyERI::N()
{
    return basisset_->nbf() * basisset_->nbf();
}
void CholeskyERI::compute_diagonal(double* target)
{
    //outfile->Printf("\n Compute Diagonal integrals: ");
    Timer diagonal_ints;
    const double* buffer = integral_->buffer();
    //int nthread = 1;
    //#ifdef _OPENMP
    //    nthread = omp_get_max_threads();
    //#endif
    //std::vector<const double*> buffer;
    //for(int thread = 0; thread < nthread; thread++)
    //{
    //    buffer.push_back(integral_threads_[thread]->buffer());
    //}

    //#pragma omp parallel for schedule(dynamic) num_threads(nthread)
    for (size_t M = 0; M < basisset_->nshell(); M++) {
        for (size_t N = 0; N < basisset_->nshell(); N++) {
     //       int thread = 0;
     //       #ifdef _OPENMP
     //           thread = omp_get_thread_num();
     //       #endif

            integral_->compute_shell(M,N,M,N);

            size_t nM = basisset_->shell(M).nfunction();
            size_t nN = basisset_->shell(N).nfunction();
            size_t mstart = basisset_->shell(M).function_index();
            size_t nstart = basisset_->shell(N).function_index();

            for (size_t om = 0; om < nM; om++) {
                for (size_t on = 0; on < nN; on++) {
                    target[(om + mstart) * basisset_->nbf() + (on + nstart)] =
                        buffer[om * nN * nM * nN + on * nM * nN + om * nN + on];
                }
            }
        }
    }
//    outfile->Printf("\n Diagonal Done in %8.8f s.", diagonal_ints.get());
}
void CholeskyERI::compute_row(int row, double* target)
{
//    outfile->Printf("\n Row computed takes ");
    Timer chol_row;
    size_t r = row / basisset_->nbf();
    size_t s = row % basisset_->nbf();
    size_t R = basisset_->function_to_shell(r);
    size_t S = basisset_->function_to_shell(s);

    size_t nR = basisset_->shell(R).nfunction();
    size_t nS = basisset_->shell(S).nfunction();
    size_t rstart = basisset_->shell(R).function_index();
    size_t sstart = basisset_->shell(S).function_index();

    size_t oR = r - rstart;
    size_t os = s - sstart;
    int nshell = basisset_->nshell();

    //int nthread = 1;
    //#ifdef _OPENMP
    //    nthread = omp_get_max_threads();
    //#endif
    //std::vector<const double*> buffer;
    //for(int thread = 0; thread < nthread; thread++)
    //{
    //    buffer.push_back(integral_threads_[thread]->buffer());
    //}


    //#pragma omp parallel for schedule(dynamic) num_threads(nthread)
    const double* buffer = integral_->buffer();
    for (size_t M = 0; M < basisset_->nshell(); M++) {
        for (size_t N = M; N < basisset_->nshell(); N++) {
 //           int thread = 0;
 //           #ifdef _OPENMP
 //               thread = omp_get_thread_num();
 //           #endif
            integral_->compute_shell(M,N,R,S);

            size_t nM = basisset_->shell(M).nfunction();
            size_t nN = basisset_->shell(N).nfunction();
            size_t mstart = basisset_->shell(M).function_index();
            size_t nstart = basisset_->shell(N).function_index();

            for (size_t om = 0; om < nM; om++) {
                for (size_t on = 0; on < nN; on++) {
                    target[(om + mstart) * basisset_->nbf() + (on + nstart)] =
                    target[(on + nstart) * basisset_->nbf() + (om + mstart)] =
                        buffer[om * nN * nR * nS + on * nR * nS + oR * nS + os];
                }
            }
        }
    }
    //outfile->Printf(" %8.8f s. ", chol_row.get());
}

CholeskyMP2::CholeskyMP2(SharedMatrix Qia,
    boost::shared_ptr<Vector> eps_aocc,
    boost::shared_ptr<Vector> eps_avir,
    bool symmetric,
    double delta, unsigned long int memory) :
    Qia_(Qia), eps_aocc_(eps_aocc), eps_avir_(eps_avir),
    symmetric_(symmetric), Cholesky(delta, memory)
{
}
CholeskyMP2::~CholeskyMP2()
{
}
size_t CholeskyMP2::N()
{
    return Qia_->colspi()[0];
}
void CholeskyMP2::compute_diagonal(double* target)
{
    size_t naocc = eps_aocc_->dimpi()[0];
    size_t navir = eps_avir_->dimpi()[0];
    size_t nQ = Qia_->rowspi()[0];

    double** Qp = Qia_->pointer();
    double* eop = eps_aocc_->pointer();
    double* evp = eps_avir_->pointer();

    for (size_t i = 0, ia = 0; i < naocc; i++) {
        for (size_t a = 0; a < navir; a++, ia++) {
            target[ia] = C_DDOT(nQ,&Qp[0][ia], naocc * (ULI) navir, &Qp[0][ia], naocc * (ULI) navir) /
                (symmetric_ ? sqrt(2.0 * (evp[a] - eop[i])) : (2.0 * (evp[a] - eop[i])));
        }
    }
}
void CholeskyMP2::compute_row(int row, double* target)
{
    size_t naocc = eps_aocc_->dimpi()[0];
    size_t navir = eps_avir_->dimpi()[0];
    size_t nQ = Qia_->rowspi()[0];

    size_t j = row / navir;
    size_t b = row % navir;

    double** Qp = Qia_->pointer();
    double* eop = eps_aocc_->pointer();
    double* evp = eps_avir_->pointer();

    for (size_t i = 0, ia = 0; i < naocc; i++) {
        for (size_t a = 0; a < navir; a++, ia++) {
            target[ia] = C_DDOT(nQ,&Qp[0][ia], naocc * (ULI) navir, &Qp[0][row], naocc * (ULI) navir) /
                (symmetric_ ? sqrt(evp[a] + evp[b] - eop[i] - eop[j]) : (evp[a] + evp[b] - eop[i] - eop[j]));
        }
    }
}

CholeskyDelta::CholeskyDelta(
    boost::shared_ptr<Vector> eps_aocc,
    boost::shared_ptr<Vector> eps_avir,
    double delta, unsigned long int memory) :
    eps_aocc_(eps_aocc), eps_avir_(eps_avir),
    Cholesky(delta, memory)
{
}
CholeskyDelta::~CholeskyDelta()
{
}
size_t CholeskyDelta::N()
{
    return eps_aocc_->dimpi()[0] * eps_avir_->dimpi()[0];
}
void CholeskyDelta::compute_diagonal(double* target)
{
    size_t naocc = eps_aocc_->dimpi()[0];
    size_t navir = eps_avir_->dimpi()[0];

    double* eop = eps_aocc_->pointer();
    double* evp = eps_avir_->pointer();

    for (size_t i = 0, ia = 0; i < naocc; i++) {
        for (size_t a = 0; a < navir; a++, ia++) {
            target[ia] = 1.0 / (2.0 * (evp[a] - eop[i]));
        }
    }
}
void CholeskyDelta::compute_row(int row, double* target)
{
    size_t naocc = eps_aocc_->dimpi()[0];
    size_t navir = eps_avir_->dimpi()[0];

    size_t j = row / navir;
    size_t b = row % navir;

    double* eop = eps_aocc_->pointer();
    double* evp = eps_avir_->pointer();

    for (size_t i = 0, ia = 0; i < naocc; i++) {
        for (size_t a = 0; a < navir; a++, ia++) {
            target[ia] = 1.0 / (evp[a] + evp[b] - eop[i] - eop[j]);
        }
    }
}

CholeskyLocal::CholeskyLocal(
    SharedMatrix C,
    double delta, unsigned long int memory) :
    C_(C), Cholesky(delta, memory)
{
}
CholeskyLocal::~CholeskyLocal()
{
}
size_t CholeskyLocal::N()
{
    return C_->rowspi()[0];
}
void CholeskyLocal::compute_diagonal(double* target)
{
    size_t n = C_->rowspi()[0];
    size_t nocc = C_->colspi()[0];

    double** Cp = C_->pointer();

    for (size_t m = 0; m < n; m++) {
        target[m] = C_DDOT(nocc, Cp[m], 1, Cp[m], 1);
    }
}
void CholeskyLocal::compute_row(int row, double* target)
{
    size_t n = C_->rowspi()[0];
    size_t nocc = C_->colspi()[0];

    double** Cp = C_->pointer();

    for (size_t m = 0; m < n; m++) {
        target[m] = C_DDOT(nocc, Cp[m], 1, Cp[row], 1);
    }
}

} // Namespace psi
