#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>
#include <libqt/qt.h>
#include <math.h>
#include <limits>
#include <vector>
#include "cholesky.h"

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
    // Initial dimensions
    int n = N();
    Q_ = 0;

    // Memory constraint on rows
    int max_int = std::numeric_limits<int>::max();

    ULI max_rows_ULI = ((memory_ - n) / (2L * n));
    int max_rows = (max_rows_ULI > max_int ? max_int : max_rows_ULI);

    // Get the diagonal (Q|Q)^(0)
    double* diag = new double[n];
    compute_diagonal(diag);

    // Temporary cholesky factor
    std::vector<double*> L;

    // List of selected pivots
    std::vector<int> pivots;

    // Cholesky procedure
    while (Q_ < n) {

        // Select the pivot
        int pivot = 0;
        double Dmax = diag[0];
        for (int P = 0; P < n; P++) {
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
            throw PSIEXCEPTION("Cholesky: Memory constraints exceeded. Fire your theorist.");
        }

        // If here, we're really going to add this row
        L.push_back(new double[n]);

        // (m|Q)
        compute_row(pivot, L[Q_]);

        // [(m|Q) - L_m^P L_Q^P]
        for (int P = 0; P < Q_; P++) {
            C_DAXPY(n,-L[P][pivots[Q_]],L[P],1,L[Q_],1);
        }

        // 1/L_QQ [(m|Q) - L_m^P L_Q^P]
        C_DSCAL(n, 1.0 / L_QQ, L[Q_], 1);

        // Zero the upper triangle
        for (int P = 0; P < pivots.size(); P++) {
            L[Q_][pivots[P]] = 0.0;
        }

        // Set the pivot factor
        L[Q_][pivot] = L_QQ;

        // Update the Schur complement diagonal
        for (int P = 0; P < n; P++) {
            diag[P] -= L[Q_][P] * L[Q_][P];
        }

        // Force truly zero elements to zero
        for (int P = 0; P < pivots.size(); P++) {
            diag[pivots[P]] = 0.0;
        }

        Q_++;
    }

    // Copy into a more permanant Matrix object
    L_ = boost::shared_ptr<Matrix>(new Matrix("Partial Cholesky", Q_, n));
    double** Lp = L_->pointer();

    for (int Q = 0; Q < Q_; Q++) {
        ::memcpy(static_cast<void*>(Lp[Q]), static_cast<void*>(L[Q]), n * sizeof(double));
        delete[] L[Q];
    }
}

CholeskyMatrix::CholeskyMatrix(boost::shared_ptr<Matrix> A, double delta, unsigned long int memory) :
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
int CholeskyMatrix::N()
{
    return A_->rowspi()[0];
}
void CholeskyMatrix::compute_diagonal(double* target)
{
    int n = N();
    double** Ap = A_->pointer();
    for (int i = 0; i < n; i++) {
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
}
CholeskyERI::~CholeskyERI()
{
}
int CholeskyERI::N()
{
    return basisset_->nbf() * basisset_->nbf();
}
void CholeskyERI::compute_diagonal(double* target)
{
    const double* buffer = integral_->buffer();

    for (int M = 0; M < basisset_->nshell(); M++) {
        for (int N = 0; N < basisset_->nshell(); N++) {

            integral_->compute_shell(M,N,M,N);

            int nM = basisset_->shell(M)->nfunction();
            int nN = basisset_->shell(N)->nfunction();
            int mstart = basisset_->shell(M)->function_index();
            int nstart = basisset_->shell(N)->function_index();

            for (int om = 0; om < nM; om++) {
                for (int on = 0; on < nN; on++) {
                    target[(om + mstart) * basisset_->nbf() + (on + nstart)] =
                        buffer[om * nN * nM * nN + on * nM * nN + om * nN + on];
                }
            }
        }
    }
}
void CholeskyERI::compute_row(int row, double* target)
{
    const double* buffer = integral_->buffer();

    int r = row / basisset_->nbf();
    int s = row % basisset_->nbf();
    int R = basisset_->function_to_shell(r);
    int S = basisset_->function_to_shell(s);

    int nR = basisset_->shell(R)->nfunction();
    int nS = basisset_->shell(S)->nfunction();
    int rstart = basisset_->shell(R)->function_index();
    int sstart = basisset_->shell(S)->function_index();

    int oR = r - rstart;
    int os = s - sstart;

    for (int M = 0; M < basisset_->nshell(); M++) {
        for (int N = 0; N < basisset_->nshell(); N++) {

            integral_->compute_shell(M,N,R,S);

            int nM = basisset_->shell(M)->nfunction();
            int nN = basisset_->shell(N)->nfunction();
            int mstart = basisset_->shell(M)->function_index();
            int nstart = basisset_->shell(N)->function_index();

            for (int om = 0; om < nM; om++) {
                for (int on = 0; on < nN; on++) {
                    target[(om + mstart) * basisset_->nbf() + (on + nstart)] =
                        buffer[om * nN * nR * nS + on * nR * nS + oR * nS + os];
                }
            }
        }
    }
}

CholeskyMP2::CholeskyMP2(boost::shared_ptr<Matrix> Qia,
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
int CholeskyMP2::N()
{
    return Qia_->colspi()[0];
}
void CholeskyMP2::compute_diagonal(double* target)
{
    int naocc = eps_aocc_->dimpi()[0];
    int navir = eps_avir_->dimpi()[0];
    int nQ = Qia_->rowspi()[0];

    double** Qp = Qia_->pointer();
    double* eop = eps_aocc_->pointer();
    double* evp = eps_avir_->pointer();

    for (int i = 0, ia = 0; i < naocc; i++) {
        for (int a = 0; a < navir; a++, ia++) {
            target[ia] = C_DDOT(nQ,&Qp[0][ia], naocc * (ULI) navir, &Qp[0][ia], naocc * (ULI) navir) /
                (symmetric_ ? sqrt(2.0 * (evp[a] - eop[i])) : (2.0 * (evp[a] - eop[i])));
        }
    }
}
void CholeskyMP2::compute_row(int row, double* target)
{
    int naocc = eps_aocc_->dimpi()[0];
    int navir = eps_avir_->dimpi()[0];
    int nQ = Qia_->rowspi()[0];

    int j = row / navir;
    int b = row % navir;

    double** Qp = Qia_->pointer();
    double* eop = eps_aocc_->pointer();
    double* evp = eps_avir_->pointer();

    for (int i = 0, ia = 0; i < naocc; i++) {
        for (int a = 0; a < navir; a++, ia++) {
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
int CholeskyDelta::N()
{
    return eps_aocc_->dimpi()[0] * eps_avir_->dimpi()[0];
}
void CholeskyDelta::compute_diagonal(double* target)
{
    int naocc = eps_aocc_->dimpi()[0];
    int navir = eps_avir_->dimpi()[0];

    double* eop = eps_aocc_->pointer();
    double* evp = eps_avir_->pointer();

    for (int i = 0, ia = 0; i < naocc; i++) {
        for (int a = 0; a < navir; a++, ia++) {
            target[ia] = 1.0 / (2.0 * (evp[a] - eop[i]));
        }
    }
}
void CholeskyDelta::compute_row(int row, double* target)
{
    int naocc = eps_aocc_->dimpi()[0];
    int navir = eps_avir_->dimpi()[0];

    int j = row / navir;
    int b = row % navir;

    double* eop = eps_aocc_->pointer();
    double* evp = eps_avir_->pointer();

    for (int i = 0, ia = 0; i < naocc; i++) {
        for (int a = 0; a < navir; a++, ia++) {
            target[ia] = 1.0 / (evp[a] + evp[b] - eop[i] - eop[j]);
        }
    }
}

CholeskyLocal::CholeskyLocal(
    boost::shared_ptr<Matrix> C,
    double delta, unsigned long int memory) :
    C_(C), Cholesky(delta, memory)
{
}
CholeskyLocal::~CholeskyLocal()
{
}
int CholeskyLocal::N()
{
    return C_->rowspi()[0];
}
void CholeskyLocal::compute_diagonal(double* target)
{
    int n = C_->rowspi()[0];
    int nocc = C_->colspi()[0];

    double** Cp = C_->pointer();

    for (int m = 0; m < n; m++) {
        target[m] = C_DDOT(nocc, Cp[m], 1, Cp[m], 1);
    }
}
void CholeskyLocal::compute_row(int row, double* target)
{
    int n = C_->rowspi()[0];
    int nocc = C_->colspi()[0];

    double** Cp = C_->pointer();

    for (int m = 0; m < n; m++) {
        target[m] = C_DDOT(nocc, Cp[m], 1, Cp[row], 1);
    }
}

} // Namespace psi
