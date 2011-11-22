#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>
#include <libqt/qt.h>
#include <math.h>
#include "qr.h"

namespace psi {

QR::QR(SharedMatrix A, double delta) :
    A_(A), delta_(delta), print_(0), debug_(0)
{
    if (A_->nirrep() != 1)
        throw PSIEXCEPTION("QR: QR does not support symmetry");
    if (A_->rowspi()[0] != A_->colspi()[0])
        throw PSIEXCEPTION("QR: QR only supported for Hermitian matrices");
}
QR::~QR()
{
}
void QR::decompose()
{
    form_QR();
    form_PN();
}
void QR::form_QR()
{
    int n = A_->rowspi()[0];

    if (debug_) {

        SharedMatrix R2(new Matrix("R",n,n));
        R2->copy(A_);
        double** R2p = R2->pointer();

        boost::shared_ptr<Vector> T2(new Vector("T",n));
        double* T2p = T2->pointer();

        boost::shared_ptr<IntVector> P2(new IntVector("Pivots",n));
        int* P2p = P2->pointer();

        double work_val2 = 0;

        C_DGEQP3(n,n,R2p[0],n,P2p,T2p,&work_val2,-1);

        int lwork2 = (int) work_val2;

        double* work2 = new double[lwork2];

        C_DGEQP3(n,n,R2p[0],n,P2p,T2p,work2,lwork2);

        delete[] work2;

        R2->print();
        T2->print();
        P2->print(outfile);
    }

    // Combination of R and reflector generators
    SharedMatrix Rtemp(new Matrix("R",n,n));
    Rtemp->copy(A_);
    double** Rp = Rtemp->pointer();

    // Pivot values of reflector generators
    boost::shared_ptr<Vector> T(new Vector("T",n));
    double* Tp = T->pointer();

    // Initial pivots
    for (int l = 0; l < n; l++) {
        pivots_.push_back(l);
    }

    // Temp for pivoting columns
    double* temp = new double[n];
    // Temp for Householder vector
    double* v = new double[n];
    // Temp for Householder cumulant
    double* q = new double[n];

    int nQ = 0;

    // ==> R <== //

    // n - 1 steps naively (usually exit early due to low l2)
    for (int l = 0; l < n; l++) {

        // Dimension of the residual
        int n2 = n - l;

        // Determine pivot by column of residual with max l2 norm
        int pivot = l;
        double max_pivot = 0.0;
        for (int m = l; m < n; m++) {
            double l2 = C_DNRM2(n2,&Rp[l][m],n);
            if (max_pivot <= l2) {
                max_pivot = l2;
                pivot = m;
            }
        }

        if (debug_) {
            fprintf(outfile, "Step l = %d: pivoting to column %d\n\n", l, pivot);
        }

        // Perform pivot
        C_DCOPY(n,&Rp[0][pivot],n,temp,1);
        C_DCOPY(n,&Rp[0][l],n,&Rp[0][pivot],n);
        C_DCOPY(n,temp,1,&Rp[0][l],n);

        int temp_pivot = pivots_[pivot];
        pivots_[pivot] = pivots_[l];
        pivots_[l] = temp_pivot;

        // Termination based on l2 norm of pivot residual, not
        // pivot itself
        if (max_pivot < delta_)
            break;

        // We're adding a column
        nQ++;

        // Last column is a void reflector
        if (l == n - 1)
            break;

        // Alpha takes the opposite sign as the pivot, for stability
        double alpha = (Rp[l][l] < 0.0 ? -max_pivot : max_pivot);

        // Build the Householder vector
        C_DCOPY(n2,&Rp[l][l],n,v,1);
        v[0] += alpha;
        C_DSCAL(n2,1.0/C_DNRM2(n2,v,1),v,1);

        // Apply the Householder reflection to the rest of R
        C_DGEMV('T',n2,n2-1,1.0,&Rp[l][l+1],n,v,1,0.0,q,1);
        C_DGER(n2,n2-1,-2.0,v,1,&q[0],1,&Rp[l][l+1],n);

        // Explicitly reflect the l-th column
        Rp[l][l] = -alpha;

        // Save the Householder reflector generator for later use
        // The body of the reflector goes in the lower triangle of R,
        // the value of the pivot of the reflector goes in T
        C_DCOPY(n2-1,&v[1],1,&Rp[l+1][l],n);
        C_DSCAL(n2-1,1.0/v[0],&Rp[l+1][l],n);
        Tp[l] = 2.0 * v[0] * v[0];

        if (debug_) {
            Rtemp->print();
            T->print();
        }
    }

    // ==> Q <== //
    Q_ = SharedMatrix(new Matrix("Q",nQ,n));
    double** Qp = Q_->pointer();

    for (int i = 0; i < nQ; i++) {
        // Definitional, but not required
        Qp[i][i] = 1.0;

        // No elements below this point (null reflector)
        if (i == n - 1)
            break;

        // How many elements remain to be copied?
        int nl = n - i - 1;

        // FORTRAN order
        C_DCOPY(nl,&Rp[i+1][i],n,&Qp[i][i+1],1);
    }

    if (debug_)
        Q_->print();

    double work_query = 0.0;

    C_DORGQR(n,nQ,nQ,Qp[0],n,Tp,&work_query,-1);

    int lwork = (int) work_query;
    double* work = new double[lwork];

    C_DORGQR(n,nQ,nQ,Qp[0],n,Tp,work,lwork);

    delete[] work;

    // Repackage R

    R_ = SharedMatrix(new Matrix("R",nQ,n));
    double** R3p = R_->pointer();

    for (int i = 0; i < nQ; i++) {
        ::memcpy(static_cast<void*>(&R3p[i][i]),static_cast<void*>(&Rp[i][i]),(n - i) * sizeof(double));
    }

    if (debug_) {
        Q_->print();
        R_->print();
        fflush(outfile);
    }
}
void QR::form_PN()
{
    int nQ = Q_->rowspi()[0];
    int n = Q_->colspi()[0];

    SharedMatrix T(new Matrix("T",nQ,n));
    SharedMatrix D(new Matrix("D",nQ,nQ));
    SharedMatrix V(new Matrix("V",nQ,nQ));
    boost::shared_ptr<Vector> d(new Vector("d",nQ));

    double** Ap = A_->pointer();
    double** Qp = Q_->pointer();
    double** Tp = T->pointer();
    double** Dp = D->pointer();
    double** Vp = V->pointer();
    double*  dp = d->pointer();

    // Transform to reduced basis
    C_DGEMM('N','N',nQ,n,n,1.0,Qp[0],n,Ap[0],n,0.0,Tp[0],n);
    C_DGEMM('N','T',nQ,nQ,n,1.0,Tp[0],n,Qp[0],n,0.0,Dp[0],nQ);

    T.reset();

    // Eigendecompose in the reduced basis
    D->diagonalize(V,d);

    // Determine the number of negative/positive roots
    int nN = 0;
    int nP = 0;
    for (int i = 0; i < nQ; i++) {
        if (dp[i] < 0)
            nN++;
        else
            nP++;
    }

    // Build the P/N vectors
    N_ = SharedMatrix(new Matrix("N",nN,n));
    if (nN > 0) {
        double** Np = N_->pointer();
        for (int i = 0; i < nN; i++) {
            C_DSCAL(nQ,sqrt(fabs(dp[i])),&Vp[0][i],nQ);
        }
        C_DGEMM('T','N',nN,n,nQ,1.0,&Vp[0][0],nQ,Qp[0],n,0.0,Np[0],n);
    }

    P_ = SharedMatrix(new Matrix("P",nP,n));
    if (nP > 0) {
        double** Pp = P_->pointer();
        for (int i = nN; i < nQ; i++) {
            C_DSCAL(nQ,sqrt(dp[i]),&Vp[0][i],nQ);
        }
        C_DGEMM('T','N',nP,n,nQ,1.0,&Vp[0][nN],nQ,Qp[0],n,0.0,Pp[0],n);
    }
}

} // Namespace psi
