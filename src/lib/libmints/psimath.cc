#include <libqt/qt.h>
#include <libmints/matrix.h>
#include <libmints/vector.h>
#include <psi4-dec.h>
#include <libmints/psimath.h>

using namespace boost;

namespace psi {

/// PSI_DGBMV, a wrapper to C_DGBMV using objects
void PSI_DGBMV(int irrep, char trans, int m, int n, int kl, int ku, double alpha, SharedMatrix a, int lda, boost::shared_ptr<Vector> x, int incx, double beta, boost::shared_ptr<Vector> y, int incy)
{
    C_DGBMV(trans, m, n, kl, ku, alpha, a->pointer(irrep)[0], lda, x->pointer(irrep), incx, beta, y->pointer(irrep), incy);
}
/// PSI_DGEMM, a wrapper to C_DGEMM using objects
void PSI_DGEMM(int irrep, char transa, char transb, int m, int n, int k, double alpha, SharedMatrix a, int lda, SharedMatrix b, int ldb, double beta, SharedMatrix c, int ldc)
{
    C_DGEMM(transa, transb, m, n, k, alpha, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, beta, c->pointer(irrep)[0], ldc);
}
/// PSI_DGEMV, a wrapper to C_DGEMV using objects
void PSI_DGEMV(int irrep, char trans, int m, int n, double alpha, SharedMatrix a, int lda, boost::shared_ptr<Vector> x, int incx, double beta, boost::shared_ptr<Vector> y, int incy)
{
    C_DGEMV(trans, m, n, alpha, a->pointer(irrep)[0], lda, x->pointer(irrep), incx, beta, y->pointer(irrep), incy);
}
/// PSI_DGER, a wrapper to C_DGER using objects
void PSI_DGER(int irrep, int m, int n, double alpha, boost::shared_ptr<Vector> x, int incx, boost::shared_ptr<Vector> y, int incy, SharedMatrix a, int lda)
{
    C_DGER(m, n, alpha, x->pointer(irrep), incx, y->pointer(irrep), incy, a->pointer(irrep)[0], lda);
}
/// PSI_DSBMV, a wrapper to C_DSBMV using objects
void PSI_DSBMV(int irrep, char uplo, int n, int k, double alpha, SharedMatrix a, int lda, boost::shared_ptr<Vector> x, int incx, double beta, boost::shared_ptr<Vector> y, int incy)
{
    C_DSBMV(uplo, n, k, alpha, a->pointer(irrep)[0], lda, x->pointer(irrep), incx, beta, y->pointer(irrep), incy);
}
/// PSI_DSYMM, a wrapper to C_DSYMM using objects
void PSI_DSYMM(int irrep, char side, char uplo, int m, int n, double alpha, SharedMatrix a, int lda, SharedMatrix b, int ldb, double beta, SharedMatrix c, int ldc)
{
    C_DSYMM(side, uplo, m, n, alpha, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, beta, c->pointer(irrep)[0], ldc);
}
/// PSI_DSYMV, a wrapper to C_DSYMV using objects
void PSI_DSYMV(int irrep, char uplo, int n, double alpha, SharedMatrix a, int lda, boost::shared_ptr<Vector> x, int incx, double beta, boost::shared_ptr<Vector> y, int incy)
{
    C_DSYMV(uplo, n, alpha, a->pointer(irrep)[0], lda, x->pointer(irrep), incx, beta, y->pointer(irrep), incy);
}
/// PSI_DSYR, a wrapper to C_DSYR using objects
void PSI_DSYR(int irrep, char uplo, int n, double alpha, boost::shared_ptr<Vector> x, int incx, SharedMatrix a, int lda)
{
    C_DSYR(uplo, n, alpha, x->pointer(irrep), incx, a->pointer(irrep)[0], lda);
}
/// PSI_DSYR2, a wrapper to C_DSYR2 using objects
void PSI_DSYR2(int irrep, char uplo, int n, double alpha, boost::shared_ptr<Vector> x, int incx, boost::shared_ptr<Vector> y, int incy, SharedMatrix a, int lda)
{
    C_DSYR2(uplo, n, alpha, x->pointer(irrep), incx, y->pointer(irrep), incy, a->pointer(irrep)[0], lda);
}
/// PSI_DSYR2K, a wrapper to C_DSYR2K using objects
void PSI_DSYR2K(int irrep, char uplo, char trans, int n, int k, double alpha, SharedMatrix a, int lda, SharedMatrix b, int ldb, double beta, SharedMatrix c, int ldc)
{
    C_DSYR2K(uplo, trans, n, k, alpha, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, beta, c->pointer(irrep)[0], ldc);
}
/// PSI_DSYRK, a wrapper to C_DSYRK using objects
void PSI_DSYRK(int irrep, char uplo, char trans, int n, int k, double alpha, SharedMatrix a, int lda, double beta, SharedMatrix c, int ldc)
{
    C_DSYRK(uplo, trans, n, k, alpha, a->pointer(irrep)[0], lda, beta, c->pointer(irrep)[0], ldc);
}
/// PSI_DTBMV, a wrapper to C_DTBMV using objects
void PSI_DTBMV(int irrep, char uplo, char trans, char diag, int n, int k, SharedMatrix a, int lda, boost::shared_ptr<Vector> x, int incx)
{
    C_DTBMV(uplo, trans, diag, n, k, a->pointer(irrep)[0], lda, x->pointer(irrep), incx);
}
/// PSI_DTBSV, a wrapper to C_DTBSV using objects
void PSI_DTBSV(int irrep, char uplo, char trans, char diag, int n, int k, SharedMatrix a, int lda, boost::shared_ptr<Vector> x, int incx)
{
    C_DTBSV(uplo, trans, diag, n, k, a->pointer(irrep)[0], lda, x->pointer(irrep), incx);
}
/// PSI_DTRMM, a wrapper to C_DTRMM using objects
void PSI_DTRMM(int irrep, char side, char uplo, char transa, char diag, int m, int n, double alpha, SharedMatrix a, int lda, SharedMatrix b, int ldb)
{
    C_DTRMM(side, uplo, transa, diag, m, n, alpha, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb);
}
/// PSI_DTRMV, a wrapper to C_DTRMV using objects
void PSI_DTRMV(int irrep, char uplo, char trans, char diag, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> x, int incx)
{
    C_DTRMV(uplo, trans, diag, n, a->pointer(irrep)[0], lda, x->pointer(irrep), incx);
}
/// PSI_DTRSM, a wrapper to C_DTRSM using objects
void PSI_DTRSM(int irrep, char side, char uplo, char transa, char diag, int m, int n, double alpha, SharedMatrix a, int lda, SharedMatrix b, int ldb)
{
    C_DTRSM(side, uplo, transa, diag, m, n, alpha, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb);
}
/// PSI_DTRSV, a wrapper to C_DTRSV using objects
void PSI_DTRSV(int irrep, char uplo, char trans, char diag, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> x, int incx)
{
    C_DTRSV(uplo, trans, diag, n, a->pointer(irrep)[0], lda, x->pointer(irrep), incx);
}

/// PSI_DROT, a wrapper to C_DROT using objects
void PSI_DROT(int irrep, unsigned long int n, boost::shared_ptr<Vector> x, int incx, boost::shared_ptr<Vector> y, int incy, double c, double s)
{
    C_DROT(n, x->pointer(irrep), incx, y->pointer(irrep), incy, c, s);
}
/// PSI_DSWAP, a wrapper to C_DSWAP using objects
void PSI_DSWAP(int irrep, unsigned long int n, boost::shared_ptr<Vector> x, int incx, boost::shared_ptr<Vector> y, int incy)
{
    C_DSWAP(n, x->pointer(irrep), incx, y->pointer(irrep), incy);
}
/// PSI_DCOPY, a wrapper to C_DCOPY using objects
void PSI_DCOPY(int irrep, unsigned long int n, boost::shared_ptr<Vector> x, int incx, boost::shared_ptr<Vector> y, int incy)
{
    C_DCOPY(n, x->pointer(irrep), incx, y->pointer(irrep), incy);
}
/// PSI_DSCAL, a wrapper to C_DSCAL using objects
void PSI_DSCAL(int irrep, unsigned long int n, double alpha, boost::shared_ptr<Vector> x, int incx)
{
    C_DSCAL(n, alpha, x->pointer(irrep), incx);
}
/// PSI_DAXPY, a wrapper to C_DAXPY using objects
void PSI_DAXPY(int irrep, unsigned long int n, double alpha, boost::shared_ptr<Vector> x, int incx, boost::shared_ptr<Vector> y, int incy)
{
    C_DAXPY(n, alpha, x->pointer(irrep), incx, y->pointer(irrep), incy);
}
/// PSI_DDOT, a wrapper to C_DDOT using objects
double PSI_DDOT(int irrep, unsigned long int n, boost::shared_ptr<Vector> x, int incx, boost::shared_ptr<Vector> y, int incy)
{
    return C_DDOT(n, x->pointer(irrep), incx, y->pointer(irrep), incy);
}
/// PSI_DNRM2, a wrapper to C_DNRM2 using objects
double PSI_DNRM2(int irrep, unsigned long int n, boost::shared_ptr<Vector> x, int incx)
{
    return C_DNRM2(n, x->pointer(irrep), incx);
}
/// PSI_DASUM, a wrapper to C_DASUM using objects
double PSI_DASUM(int irrep, unsigned long int n, boost::shared_ptr<Vector> x, int incx)
{
    return C_DASUM(n, x->pointer(irrep), incx);
}
/// PSI_IDAMAX, a wrapper to C_IDAMAX using objects
unsigned long int PSI_IDAMAX(int irrep, unsigned long int n, boost::shared_ptr<Vector> x, int incx)
{
    return C_IDAMAX(n, x->pointer(irrep), incx);
}


/// LAPACK
/// PSI_DBDSDC, a wrapper to return C_DBDSDC using objects
int PSI_DBDSDC(int irrep, char uplo, char compq, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix u, int ldu, SharedMatrix vt, int ldvt, boost::shared_ptr<Vector> q, boost::shared_ptr<IntVector> iq, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DBDSDC(uplo, compq, n, d->pointer(irrep), e->pointer(irrep), u->pointer(irrep)[0], ldu, vt->pointer(irrep)[0], ldvt, q->pointer(irrep), iq->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DBDSQR, a wrapper to return C_DBDSQR using objects
int PSI_DBDSQR(int irrep, char uplo, int n, int ncvt, int nru, int ncc, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix vt, int ldvt, SharedMatrix u, int ldu, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work)
{
    return C_DBDSQR(uplo, n, ncvt, nru, ncc, d->pointer(irrep), e->pointer(irrep), vt->pointer(irrep)[0], ldvt, u->pointer(irrep)[0], ldu, c->pointer(irrep)[0], ldc, work->pointer(irrep));
}
/// PSI_DDISNA, a wrapper to return C_DDISNA using objects
int PSI_DDISNA(int irrep, char job, int m, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> sep)
{
    return C_DDISNA(job, m, n, d->pointer(irrep), sep->pointer(irrep));
}
/// PSI_DGBBRD, a wrapper to return C_DGBBRD using objects
int PSI_DGBBRD(int irrep, char vect, int m, int n, int ncc, int kl, int ku, SharedMatrix ab, int ldab, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix q, int ldq, SharedMatrix pt, int ldpt, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work)
{
    return C_DGBBRD(vect, m, n, ncc, kl, ku, ab->pointer(irrep)[0], ldab, d->pointer(irrep), e->pointer(irrep), q->pointer(irrep)[0], ldq, pt->pointer(irrep)[0], ldpt, c->pointer(irrep)[0], ldc, work->pointer(irrep));
}
/// PSI_DGBCON, a wrapper to return C_DGBCON using objects
int PSI_DGBCON(int irrep, char norm, int n, int kl, int ku, SharedMatrix ab, int ldab, boost::shared_ptr<IntVector> ipiv, double anorm, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DGBCON(norm, n, kl, ku, ab->pointer(irrep)[0], ldab, ipiv->pointer(irrep), anorm, rcond->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DGBEQU, a wrapper to return C_DGBEQU using objects
int PSI_DGBEQU(int irrep, int m, int n, int kl, int ku, SharedMatrix ab, int ldab, boost::shared_ptr<Vector> r, boost::shared_ptr<Vector> c, boost::shared_ptr<Vector> rowcnd, boost::shared_ptr<Vector> colcnd, boost::shared_ptr<Vector> amax)
{
    return C_DGBEQU(m, n, kl, ku, ab->pointer(irrep)[0], ldab, r->pointer(irrep), c->pointer(irrep), rowcnd->pointer(irrep), colcnd->pointer(irrep), amax->pointer(irrep));
}
/// PSI_DGBRFS, a wrapper to return C_DGBRFS using objects
int PSI_DGBRFS(int irrep, char trans, int n, int kl, int ku, int nrhs, SharedMatrix ab, int ldab, SharedMatrix afb, int ldafb, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DGBRFS(trans, n, kl, ku, nrhs, ab->pointer(irrep)[0], ldab, afb->pointer(irrep)[0], ldafb, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DGBSV, a wrapper to return C_DGBSV using objects
int PSI_DGBSV(int irrep, int n, int kl, int ku, int nrhs, SharedMatrix ab, int ldab, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb)
{
    return C_DGBSV(n, kl, ku, nrhs, ab->pointer(irrep)[0], ldab, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DGBSVX, a wrapper to return C_DGBSVX using objects
int PSI_DGBSVX(int irrep, char fact, char trans, int n, int kl, int ku, int nrhs, SharedMatrix ab, int ldab, SharedMatrix afb, int ldafb, boost::shared_ptr<IntVector> ipiv, char equed, boost::shared_ptr<Vector> r, boost::shared_ptr<Vector> c, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DGBSVX(fact, trans, n, kl, ku, nrhs, ab->pointer(irrep)[0], ldab, afb->pointer(irrep)[0], ldafb, ipiv->pointer(irrep), equed, r->pointer(irrep), c->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, rcond->pointer(irrep), ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DGBTRF, a wrapper to return C_DGBTRF using objects
int PSI_DGBTRF(int irrep, int m, int n, int kl, int ku, SharedMatrix ab, int ldab, boost::shared_ptr<IntVector> ipiv)
{
    return C_DGBTRF(m, n, kl, ku, ab->pointer(irrep)[0], ldab, ipiv->pointer(irrep));
}
/// PSI_DGBTRS, a wrapper to return C_DGBTRS using objects
int PSI_DGBTRS(int irrep, char trans, int n, int kl, int ku, int nrhs, SharedMatrix ab, int ldab, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb)
{
    return C_DGBTRS(trans, n, kl, ku, nrhs, ab->pointer(irrep)[0], ldab, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DGEBAK, a wrapper to return C_DGEBAK using objects
int PSI_DGEBAK(int irrep, char job, char side, int n, int ilo, int ihi, boost::shared_ptr<Vector> scale, int m, SharedMatrix v, int ldv)
{
    return C_DGEBAK(job, side, n, ilo, ihi, scale->pointer(irrep), m, v->pointer(irrep)[0], ldv);
}
/// PSI_DGEBAL, a wrapper to return C_DGEBAL using objects
int PSI_DGEBAL(int irrep, char job, int n, SharedMatrix a, int lda, boost::shared_ptr<IntVector> ilo, boost::shared_ptr<IntVector> ihi, boost::shared_ptr<Vector> scale)
{
    return C_DGEBAL(job, n, a->pointer(irrep)[0], lda, ilo->pointer(irrep), ihi->pointer(irrep), scale->pointer(irrep));
}
/// PSI_DGEBRD, a wrapper to return C_DGEBRD using objects
int PSI_DGEBRD(int irrep, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, boost::shared_ptr<Vector> tauq, boost::shared_ptr<Vector> taup, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGEBRD(m, n, a->pointer(irrep)[0], lda, d->pointer(irrep), e->pointer(irrep), tauq->pointer(irrep), taup->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGECON, a wrapper to return C_DGECON using objects
int PSI_DGECON(int irrep, char norm, int n, SharedMatrix a, int lda, double anorm, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DGECON(norm, n, a->pointer(irrep)[0], lda, anorm, rcond->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DGEEQU, a wrapper to return C_DGEEQU using objects
int PSI_DGEEQU(int irrep, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> r, boost::shared_ptr<Vector> c, boost::shared_ptr<Vector> rowcnd, boost::shared_ptr<Vector> colcnd, boost::shared_ptr<Vector> amax)
{
    return C_DGEEQU(m, n, a->pointer(irrep)[0], lda, r->pointer(irrep), c->pointer(irrep), rowcnd->pointer(irrep), colcnd->pointer(irrep), amax->pointer(irrep));
}
/// PSI_DGEES, a wrapper to return C_DGEES using objects
int PSI_DGEES(int irrep, char jobvs, char sort, int n, SharedMatrix a, int lda, boost::shared_ptr<IntVector> sdim, boost::shared_ptr<Vector> wr, boost::shared_ptr<Vector> wi, SharedMatrix vs, int ldvs, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGEES(jobvs, sort, n, a->pointer(irrep)[0], lda, sdim->pointer(irrep), wr->pointer(irrep), wi->pointer(irrep), vs->pointer(irrep)[0], ldvs, work->pointer(irrep), lwork);
}
/// PSI_DGEESX, a wrapper to return C_DGEESX using objects
int PSI_DGEESX(int irrep, char jobvs, char sort, char sense, int n, SharedMatrix a, int lda, boost::shared_ptr<IntVector> sdim, boost::shared_ptr<Vector> wr, boost::shared_ptr<Vector> wi, SharedMatrix vs, int ldvs, boost::shared_ptr<Vector> rconde, boost::shared_ptr<Vector> rcondv, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DGEESX(jobvs, sort, sense, n, a->pointer(irrep)[0], lda, sdim->pointer(irrep), wr->pointer(irrep), wi->pointer(irrep), vs->pointer(irrep)[0], ldvs, rconde->pointer(irrep), rcondv->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DGEEV, a wrapper to return C_DGEEV using objects
int PSI_DGEEV(int irrep, char jobvl, char jobvr, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> wr, boost::shared_ptr<Vector> wi, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGEEV(jobvl, jobvr, n, a->pointer(irrep)[0], lda, wr->pointer(irrep), wi->pointer(irrep), vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr, work->pointer(irrep), lwork);
}
/// PSI_DGEEVX, a wrapper to return C_DGEEVX using objects
int PSI_DGEEVX(int irrep, char balanc, char jobvl, char jobvr, char sense, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> wr, boost::shared_ptr<Vector> wi, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, boost::shared_ptr<IntVector> ilo, boost::shared_ptr<IntVector> ihi, boost::shared_ptr<Vector> scale, boost::shared_ptr<Vector> abnrm, boost::shared_ptr<Vector> rconde, boost::shared_ptr<Vector> rcondv, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork)
{
    return C_DGEEVX(balanc, jobvl, jobvr, sense, n, a->pointer(irrep)[0], lda, wr->pointer(irrep), wi->pointer(irrep), vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr, ilo->pointer(irrep), ihi->pointer(irrep), scale->pointer(irrep), abnrm->pointer(irrep), rconde->pointer(irrep), rcondv->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep));
}
/// PSI_DGEGS, a wrapper to return C_DGEGS using objects
int PSI_DGEGS(int irrep, char jobvsl, char jobvsr, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> alphar, boost::shared_ptr<Vector> alphai, boost::shared_ptr<Vector> beta, SharedMatrix vsl, int ldvsl, SharedMatrix vsr, int ldvsr, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGEGS(jobvsl, jobvsr, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, alphar->pointer(irrep), alphai->pointer(irrep), beta->pointer(irrep), vsl->pointer(irrep)[0], ldvsl, vsr->pointer(irrep)[0], ldvsr, work->pointer(irrep), lwork);
}
/// PSI_DGEGV, a wrapper to return C_DGEGV using objects
int PSI_DGEGV(int irrep, char jobvl, char jobvr, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> alphar, boost::shared_ptr<Vector> alphai, boost::shared_ptr<Vector> beta, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGEGV(jobvl, jobvr, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, alphar->pointer(irrep), alphai->pointer(irrep), beta->pointer(irrep), vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr, work->pointer(irrep), lwork);
}
/// PSI_DGEHRD, a wrapper to return C_DGEHRD using objects
int PSI_DGEHRD(int irrep, int n, int ilo, int ihi, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGEHRD(n, ilo, ihi, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGELQF, a wrapper to return C_DGELQF using objects
int PSI_DGELQF(int irrep, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGELQF(m, n, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGELS, a wrapper to return C_DGELS using objects
int PSI_DGELS(int irrep, char trans, int m, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGELS(trans, m, n, nrhs, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, work->pointer(irrep), lwork);
}
/// PSI_DGELSD, a wrapper to return C_DGELSD using objects
int PSI_DGELSD(int irrep, int m, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> s, double rcond, boost::shared_ptr<IntVector> rank, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork)
{
    return C_DGELSD(m, n, nrhs, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, s->pointer(irrep), rcond, rank->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep));
}
/// PSI_DGELSS, a wrapper to return C_DGELSS using objects
int PSI_DGELSS(int irrep, int m, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> s, double rcond, boost::shared_ptr<IntVector> rank, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGELSS(m, n, nrhs, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, s->pointer(irrep), rcond, rank->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGELSX, a wrapper to return C_DGELSX using objects
int PSI_DGELSX(int irrep, int m, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<IntVector> jpvt, double rcond, boost::shared_ptr<IntVector> rank, boost::shared_ptr<Vector> work)
{
    return C_DGELSX(m, n, nrhs, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, jpvt->pointer(irrep), rcond, rank->pointer(irrep), work->pointer(irrep));
}
/// PSI_DGELSY, a wrapper to return C_DGELSY using objects
int PSI_DGELSY(int irrep, int m, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<IntVector> jpvt, double rcond, boost::shared_ptr<IntVector> rank, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGELSY(m, n, nrhs, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, jpvt->pointer(irrep), rcond, rank->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGEQLF, a wrapper to return C_DGEQLF using objects
int PSI_DGEQLF(int irrep, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGEQLF(m, n, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGEQP3, a wrapper to return C_DGEQP3 using objects
int PSI_DGEQP3(int irrep, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<IntVector> jpvt, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGEQP3(m, n, a->pointer(irrep)[0], lda, jpvt->pointer(irrep), tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGEQPF, a wrapper to return C_DGEQPF using objects
int PSI_DGEQPF(int irrep, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<IntVector> jpvt, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work)
{
    return C_DGEQPF(m, n, a->pointer(irrep)[0], lda, jpvt->pointer(irrep), tau->pointer(irrep), work->pointer(irrep));
}
/// PSI_DGERFS, a wrapper to return C_DGERFS using objects
int PSI_DGERFS(int irrep, char trans, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DGERFS(trans, n, nrhs, a->pointer(irrep)[0], lda, af->pointer(irrep)[0], ldaf, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DGERQF, a wrapper to return C_DGERQF using objects
int PSI_DGERQF(int irrep, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGERQF(m, n, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGESDD, a wrapper to return C_DGESDD using objects
int PSI_DGESDD(int irrep, char jobz, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> s, SharedMatrix u, int ldu, SharedMatrix vt, int ldvt, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork)
{
    return C_DGESDD(jobz, m, n, a->pointer(irrep)[0], lda, s->pointer(irrep), u->pointer(irrep)[0], ldu, vt->pointer(irrep)[0], ldvt, work->pointer(irrep), lwork, iwork->pointer(irrep));
}
/// PSI_DGESV, a wrapper to return C_DGESV using objects
int PSI_DGESV(int irrep, int n, int nrhs, SharedMatrix a, int lda, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb)
{
    return C_DGESV(n, nrhs, a->pointer(irrep)[0], lda, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DGESVX, a wrapper to return C_DGESVX using objects
int PSI_DGESVX(int irrep, char fact, char trans, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf, boost::shared_ptr<IntVector> ipiv, char equed, boost::shared_ptr<Vector> r, boost::shared_ptr<Vector> c, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DGESVX(fact, trans, n, nrhs, a->pointer(irrep)[0], lda, af->pointer(irrep)[0], ldaf, ipiv->pointer(irrep), equed, r->pointer(irrep), c->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, rcond->pointer(irrep), ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DGETRF, a wrapper to return C_DGETRF using objects
int PSI_DGETRF(int irrep, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<IntVector> ipiv)
{
    return C_DGETRF(m, n, a->pointer(irrep)[0], lda, ipiv->pointer(irrep));
}
/// PSI_DGETRI, a wrapper to return C_DGETRI using objects
int PSI_DGETRI(int irrep, int n, SharedMatrix a, int lda, boost::shared_ptr<IntVector> ipiv, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGETRI(n, a->pointer(irrep)[0], lda, ipiv->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGETRS, a wrapper to return C_DGETRS using objects
int PSI_DGETRS(int irrep, char trans, int n, int nrhs, SharedMatrix a, int lda, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb)
{
    return C_DGETRS(trans, n, nrhs, a->pointer(irrep)[0], lda, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DGGBAK, a wrapper to return C_DGGBAK using objects
int PSI_DGGBAK(int irrep, char job, char side, int n, int ilo, int ihi, boost::shared_ptr<Vector> lscale, boost::shared_ptr<Vector> rscale, int m, SharedMatrix v, int ldv)
{
    return C_DGGBAK(job, side, n, ilo, ihi, lscale->pointer(irrep), rscale->pointer(irrep), m, v->pointer(irrep)[0], ldv);
}
/// PSI_DGGBAL, a wrapper to return C_DGGBAL using objects
int PSI_DGGBAL(int irrep, char job, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<IntVector> ilo, boost::shared_ptr<IntVector> ihi, boost::shared_ptr<Vector> lscale, boost::shared_ptr<Vector> rscale, boost::shared_ptr<Vector> work)
{
    return C_DGGBAL(job, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, ilo->pointer(irrep), ihi->pointer(irrep), lscale->pointer(irrep), rscale->pointer(irrep), work->pointer(irrep));
}
/// PSI_DGGES, a wrapper to return C_DGGES using objects
int PSI_DGGES(int irrep, char jobvsl, char jobvsr, char sort, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<IntVector> sdim, boost::shared_ptr<Vector> alphar, boost::shared_ptr<Vector> alphai, boost::shared_ptr<Vector> beta, SharedMatrix vsl, int ldvsl, SharedMatrix vsr, int ldvsr, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGGES(jobvsl, jobvsr, sort, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, sdim->pointer(irrep), alphar->pointer(irrep), alphai->pointer(irrep), beta->pointer(irrep), vsl->pointer(irrep)[0], ldvsl, vsr->pointer(irrep)[0], ldvsr, work->pointer(irrep), lwork);
}
/// PSI_DGGESX, a wrapper to return C_DGGESX using objects
int PSI_DGGESX(int irrep, char jobvsl, char jobvsr, char sort, char sense, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<IntVector> sdim, boost::shared_ptr<Vector> alphar, boost::shared_ptr<Vector> alphai, boost::shared_ptr<Vector> beta, SharedMatrix vsl, int ldvsl, SharedMatrix vsr, int ldvsr, boost::shared_ptr<Vector> rconde, boost::shared_ptr<Vector> rcondv, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DGGESX(jobvsl, jobvsr, sort, sense, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, sdim->pointer(irrep), alphar->pointer(irrep), alphai->pointer(irrep), beta->pointer(irrep), vsl->pointer(irrep)[0], ldvsl, vsr->pointer(irrep)[0], ldvsr, rconde->pointer(irrep), rcondv->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DGGEV, a wrapper to return C_DGGEV using objects
int PSI_DGGEV(int irrep, char jobvl, char jobvr, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> alphar, boost::shared_ptr<Vector> alphai, boost::shared_ptr<Vector> beta, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGGEV(jobvl, jobvr, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, alphar->pointer(irrep), alphai->pointer(irrep), beta->pointer(irrep), vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr, work->pointer(irrep), lwork);
}
/// PSI_DGGEVX, a wrapper to return C_DGGEVX using objects
int PSI_DGGEVX(int irrep, char balanc, char jobvl, char jobvr, char sense, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> alphar, boost::shared_ptr<Vector> alphai, boost::shared_ptr<Vector> beta, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, boost::shared_ptr<IntVector> ilo, boost::shared_ptr<IntVector> ihi, boost::shared_ptr<Vector> lscale, boost::shared_ptr<Vector> rscale, boost::shared_ptr<Vector> abnrm, boost::shared_ptr<Vector> bbnrm, boost::shared_ptr<Vector> rconde, boost::shared_ptr<Vector> rcondv, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork)
{
    return C_DGGEVX(balanc, jobvl, jobvr, sense, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, alphar->pointer(irrep), alphai->pointer(irrep), beta->pointer(irrep), vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr, ilo->pointer(irrep), ihi->pointer(irrep), lscale->pointer(irrep), rscale->pointer(irrep), abnrm->pointer(irrep), bbnrm->pointer(irrep), rconde->pointer(irrep), rcondv->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep));
}
/// PSI_DGGGLM, a wrapper to return C_DGGGLM using objects
int PSI_DGGGLM(int irrep, int n, int m, int p, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> x, boost::shared_ptr<Vector> y, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGGGLM(n, m, p, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, d->pointer(irrep), x->pointer(irrep), y->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGGHRD, a wrapper to return C_DGGHRD using objects
int PSI_DGGHRD(int irrep, char compq, char compz, int n, int ilo, int ihi, SharedMatrix a, int lda, SharedMatrix b, int ldb, SharedMatrix q, int ldq, SharedMatrix z, int ldz)
{
    return C_DGGHRD(compq, compz, n, ilo, ihi, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, q->pointer(irrep)[0], ldq, z->pointer(irrep)[0], ldz);
}
/// PSI_DGGLSE, a wrapper to return C_DGGLSE using objects
int PSI_DGGLSE(int irrep, int m, int n, int p, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> c, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> x, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGGLSE(m, n, p, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, c->pointer(irrep), d->pointer(irrep), x->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGGQRF, a wrapper to return C_DGGQRF using objects
int PSI_DGGQRF(int irrep, int n, int m, int p, SharedMatrix a, int lda, boost::shared_ptr<Vector> taua, SharedMatrix b, int ldb, boost::shared_ptr<Vector> taub, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGGQRF(n, m, p, a->pointer(irrep)[0], lda, taua->pointer(irrep), b->pointer(irrep)[0], ldb, taub->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGGRQF, a wrapper to return C_DGGRQF using objects
int PSI_DGGRQF(int irrep, int m, int p, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> taua, SharedMatrix b, int ldb, boost::shared_ptr<Vector> taub, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DGGRQF(m, p, n, a->pointer(irrep)[0], lda, taua->pointer(irrep), b->pointer(irrep)[0], ldb, taub->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DGGSVD, a wrapper to return C_DGGSVD using objects
int PSI_DGGSVD(int irrep, char jobu, char jobv, char jobq, int m, int n, int p, boost::shared_ptr<IntVector> k, boost::shared_ptr<IntVector> l, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> alpha, boost::shared_ptr<Vector> beta, SharedMatrix u, int ldu, SharedMatrix v, int ldv, SharedMatrix q, int ldq, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DGGSVD(jobu, jobv, jobq, m, n, p, k->pointer(irrep), l->pointer(irrep), a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, alpha->pointer(irrep), beta->pointer(irrep), u->pointer(irrep)[0], ldu, v->pointer(irrep)[0], ldv, q->pointer(irrep)[0], ldq, work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DGGSVP, a wrapper to return C_DGGSVP using objects
int PSI_DGGSVP(int irrep, char jobu, char jobv, char jobq, int m, int p, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, double tola, double tolb, boost::shared_ptr<IntVector> k, boost::shared_ptr<IntVector> l, SharedMatrix u, int ldu, SharedMatrix v, int ldv, SharedMatrix q, int ldq, boost::shared_ptr<IntVector> iwork, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work)
{
    return C_DGGSVP(jobu, jobv, jobq, m, p, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, tola, tolb, k->pointer(irrep), l->pointer(irrep), u->pointer(irrep)[0], ldu, v->pointer(irrep)[0], ldv, q->pointer(irrep)[0], ldq, iwork->pointer(irrep), tau->pointer(irrep), work->pointer(irrep));
}
/// PSI_DGTCON, a wrapper to return C_DGTCON using objects
int PSI_DGTCON(int irrep, char norm, int n, boost::shared_ptr<Vector> dl, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> du, boost::shared_ptr<Vector> du2, boost::shared_ptr<IntVector> ipiv, double anorm, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DGTCON(norm, n, dl->pointer(irrep), d->pointer(irrep), du->pointer(irrep), du2->pointer(irrep), ipiv->pointer(irrep), anorm, rcond->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DGTRFS, a wrapper to return C_DGTRFS using objects
int PSI_DGTRFS(int irrep, char trans, int n, int nrhs, boost::shared_ptr<Vector> dl, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> du, boost::shared_ptr<Vector> dlf, boost::shared_ptr<Vector> df, boost::shared_ptr<Vector> duf, boost::shared_ptr<Vector> du2, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DGTRFS(trans, n, nrhs, dl->pointer(irrep), d->pointer(irrep), du->pointer(irrep), dlf->pointer(irrep), df->pointer(irrep), duf->pointer(irrep), du2->pointer(irrep), ipiv->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DGTSV, a wrapper to return C_DGTSV using objects
int PSI_DGTSV(int irrep, int n, int nrhs, boost::shared_ptr<Vector> dl, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> du, SharedMatrix b, int ldb)
{
    return C_DGTSV(n, nrhs, dl->pointer(irrep), d->pointer(irrep), du->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DGTSVX, a wrapper to return C_DGTSVX using objects
int PSI_DGTSVX(int irrep, char fact, char trans, int n, int nrhs, boost::shared_ptr<Vector> dl, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> du, boost::shared_ptr<Vector> dlf, boost::shared_ptr<Vector> df, boost::shared_ptr<Vector> duf, boost::shared_ptr<Vector> du2, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> rcond)
{
    return C_DGTSVX(fact, trans, n, nrhs, dl->pointer(irrep), d->pointer(irrep), du->pointer(irrep), dlf->pointer(irrep), df->pointer(irrep), duf->pointer(irrep), du2->pointer(irrep), ipiv->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, rcond->pointer(irrep));
}
/// PSI_DGTTRF, a wrapper to return C_DGTTRF using objects
int PSI_DGTTRF(int irrep, int n, boost::shared_ptr<Vector> dl, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> du, boost::shared_ptr<Vector> du2, boost::shared_ptr<IntVector> ipiv)
{
    return C_DGTTRF(n, dl->pointer(irrep), d->pointer(irrep), du->pointer(irrep), du2->pointer(irrep), ipiv->pointer(irrep));
}
/// PSI_DGTTRS, a wrapper to return C_DGTTRS using objects
int PSI_DGTTRS(int irrep, char trans, int n, int nrhs, boost::shared_ptr<Vector> dl, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> du, boost::shared_ptr<Vector> du2, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb)
{
    return C_DGTTRS(trans, n, nrhs, dl->pointer(irrep), d->pointer(irrep), du->pointer(irrep), du2->pointer(irrep), ipiv->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DHGEQZ, a wrapper to return C_DHGEQZ using objects
int PSI_DHGEQZ(int irrep, char job, char compq, char compz, int n, int ilo, int ihi, SharedMatrix h, int ldh, SharedMatrix t, int ldt, boost::shared_ptr<Vector> alphar, boost::shared_ptr<Vector> alphai, boost::shared_ptr<Vector> beta, SharedMatrix q, int ldq, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DHGEQZ(job, compq, compz, n, ilo, ihi, h->pointer(irrep)[0], ldh, t->pointer(irrep)[0], ldt, alphar->pointer(irrep), alphai->pointer(irrep), beta->pointer(irrep), q->pointer(irrep)[0], ldq, z->pointer(irrep)[0], ldz, work->pointer(irrep), lwork);
}
/// PSI_DHSEIN, a wrapper to return C_DHSEIN using objects
int PSI_DHSEIN(int irrep, char side, char eigsrc, char initv, int n, SharedMatrix h, int ldh, boost::shared_ptr<Vector> wr, boost::shared_ptr<Vector> wi, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, int mm, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> ifaill, boost::shared_ptr<IntVector> ifailr)
{
    return C_DHSEIN(side, eigsrc, initv, n, h->pointer(irrep)[0], ldh, wr->pointer(irrep), wi->pointer(irrep), vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr, mm, m->pointer(irrep), work->pointer(irrep), ifaill->pointer(irrep), ifailr->pointer(irrep));
}
/// PSI_DHSEQR, a wrapper to return C_DHSEQR using objects
int PSI_DHSEQR(int irrep, char job, char compz, int n, int ilo, int ihi, SharedMatrix h, int ldh, boost::shared_ptr<Vector> wr, boost::shared_ptr<Vector> wi, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DHSEQR(job, compz, n, ilo, ihi, h->pointer(irrep)[0], ldh, wr->pointer(irrep), wi->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), lwork);
}
/// PSI_DORGBR, a wrapper to return C_DORGBR using objects
int PSI_DORGBR(int irrep, char vect, int m, int n, int k, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORGBR(vect, m, n, k, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DORGHR, a wrapper to return C_DORGHR using objects
int PSI_DORGHR(int irrep, int n, int ilo, int ihi, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORGHR(n, ilo, ihi, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DORGLQ, a wrapper to return C_DORGLQ using objects
int PSI_DORGLQ(int irrep, int m, int n, int k, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORGLQ(m, n, k, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DORGQL, a wrapper to return C_DORGQL using objects
int PSI_DORGQL(int irrep, int m, int n, int k, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORGQL(m, n, k, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DORGQR, a wrapper to return C_DORGQR using objects
int PSI_DORGQR(int irrep, int m, int n, int k, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORGQR(m, n, k, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DORGRQ, a wrapper to return C_DORGRQ using objects
int PSI_DORGRQ(int irrep, int m, int n, int k, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORGRQ(m, n, k, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DORGTR, a wrapper to return C_DORGTR using objects
int PSI_DORGTR(int irrep, char uplo, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORGTR(uplo, n, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DORMBR, a wrapper to return C_DORMBR using objects
int PSI_DORMBR(int irrep, char vect, char side, char trans, int m, int n, int k, boost::shared_ptr<Vector> a, int lda, boost::shared_ptr<Vector> tau, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORMBR(vect, side, trans, m, n, k, a->pointer(irrep), lda, tau->pointer(irrep), c->pointer(irrep)[0], ldc, work->pointer(irrep), lwork);
}
/// PSI_DORMHR, a wrapper to return C_DORMHR using objects
int PSI_DORMHR(int irrep, char side, char trans, int m, int n, int ilo, int ihi, boost::shared_ptr<Vector> a, int lda, boost::shared_ptr<Vector> tau, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORMHR(side, trans, m, n, ilo, ihi, a->pointer(irrep), lda, tau->pointer(irrep), c->pointer(irrep)[0], ldc, work->pointer(irrep), lwork);
}
/// PSI_DORMLQ, a wrapper to return C_DORMLQ using objects
int PSI_DORMLQ(int irrep, char side, char trans, int m, int n, int k, boost::shared_ptr<Vector> a, int lda, boost::shared_ptr<Vector> tau, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORMLQ(side, trans, m, n, k, a->pointer(irrep), lda, tau->pointer(irrep), c->pointer(irrep)[0], ldc, work->pointer(irrep), lwork);
}
/// PSI_DORMQL, a wrapper to return C_DORMQL using objects
int PSI_DORMQL(int irrep, char side, char trans, int m, int n, int k, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORMQL(side, trans, m, n, k, a->pointer(irrep)[0], lda, tau->pointer(irrep), c->pointer(irrep)[0], ldc, work->pointer(irrep), lwork);
}
/// PSI_DORMQR, a wrapper to return C_DORMQR using objects
int PSI_DORMQR(int irrep, char side, char trans, int m, int n, int k, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORMQR(side, trans, m, n, k, a->pointer(irrep)[0], lda, tau->pointer(irrep), c->pointer(irrep)[0], ldc, work->pointer(irrep), lwork);
}
/// PSI_DORMR3, a wrapper to return C_DORMR3 using objects
int PSI_DORMR3(int irrep, char side, char trans, int m, int n, int k, int l, boost::shared_ptr<Vector> a, int lda, boost::shared_ptr<Vector> tau, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work)
{
    return C_DORMR3(side, trans, m, n, k, l, a->pointer(irrep), lda, tau->pointer(irrep), c->pointer(irrep)[0], ldc, work->pointer(irrep));
}
/// PSI_DORMRQ, a wrapper to return C_DORMRQ using objects
int PSI_DORMRQ(int irrep, char side, char trans, int m, int n, int k, boost::shared_ptr<Vector> a, int lda, boost::shared_ptr<Vector> tau, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORMRQ(side, trans, m, n, k, a->pointer(irrep), lda, tau->pointer(irrep), c->pointer(irrep)[0], ldc, work->pointer(irrep), lwork);
}
/// PSI_DORMRZ, a wrapper to return C_DORMRZ using objects
int PSI_DORMRZ(int irrep, char side, char trans, int m, int n, int k, int l, boost::shared_ptr<Vector> a, int lda, boost::shared_ptr<Vector> tau, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORMRZ(side, trans, m, n, k, l, a->pointer(irrep), lda, tau->pointer(irrep), c->pointer(irrep)[0], ldc, work->pointer(irrep), lwork);
}
/// PSI_DORMTR, a wrapper to return C_DORMTR using objects
int PSI_DORMTR(int irrep, char side, char uplo, char trans, int m, int n, boost::shared_ptr<Vector> a, int lda, boost::shared_ptr<Vector> tau, SharedMatrix c, int ldc, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DORMTR(side, uplo, trans, m, n, a->pointer(irrep), lda, tau->pointer(irrep), c->pointer(irrep)[0], ldc, work->pointer(irrep), lwork);
}
/// PSI_DPBCON, a wrapper to return C_DPBCON using objects
int PSI_DPBCON(int irrep, char uplo, int n, int kd, SharedMatrix ab, int ldab, double anorm, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DPBCON(uplo, n, kd, ab->pointer(irrep)[0], ldab, anorm, rcond->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DPBEQU, a wrapper to return C_DPBEQU using objects
int PSI_DPBEQU(int irrep, char uplo, int n, int kd, SharedMatrix ab, int ldab, boost::shared_ptr<Vector> s, boost::shared_ptr<Vector> scond, boost::shared_ptr<Vector> amax)
{
    return C_DPBEQU(uplo, n, kd, ab->pointer(irrep)[0], ldab, s->pointer(irrep), scond->pointer(irrep), amax->pointer(irrep));
}
/// PSI_DPBRFS, a wrapper to return C_DPBRFS using objects
int PSI_DPBRFS(int irrep, char uplo, int n, int kd, int nrhs, SharedMatrix ab, int ldab, SharedMatrix afb, int ldafb, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DPBRFS(uplo, n, kd, nrhs, ab->pointer(irrep)[0], ldab, afb->pointer(irrep)[0], ldafb, b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DPBSTF, a wrapper to return C_DPBSTF using objects
int PSI_DPBSTF(int irrep, char uplo, int n, int kd, SharedMatrix ab, int ldab)
{
    return C_DPBSTF(uplo, n, kd, ab->pointer(irrep)[0], ldab);
}
/// PSI_DPBSV, a wrapper to return C_DPBSV using objects
int PSI_DPBSV(int irrep, char uplo, int n, int kd, int nrhs, SharedMatrix ab, int ldab, SharedMatrix b, int ldb)
{
    return C_DPBSV(uplo, n, kd, nrhs, ab->pointer(irrep)[0], ldab, b->pointer(irrep)[0], ldb);
}
/// PSI_DPBSVX, a wrapper to return C_DPBSVX using objects
int PSI_DPBSVX(int irrep, char fact, char uplo, int n, int kd, int nrhs, SharedMatrix ab, int ldab, SharedMatrix afb, int ldafb, char equed, boost::shared_ptr<Vector> s, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DPBSVX(fact, uplo, n, kd, nrhs, ab->pointer(irrep)[0], ldab, afb->pointer(irrep)[0], ldafb, equed, s->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, rcond->pointer(irrep), ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DPBTRF, a wrapper to return C_DPBTRF using objects
int PSI_DPBTRF(int irrep, char uplo, int n, int kd, SharedMatrix ab, int ldab)
{
    return C_DPBTRF(uplo, n, kd, ab->pointer(irrep)[0], ldab);
}
/// PSI_DPBTRS, a wrapper to return C_DPBTRS using objects
int PSI_DPBTRS(int irrep, char uplo, int n, int kd, int nrhs, SharedMatrix ab, int ldab, SharedMatrix b, int ldb)
{
    return C_DPBTRS(uplo, n, kd, nrhs, ab->pointer(irrep)[0], ldab, b->pointer(irrep)[0], ldb);
}
/// PSI_DPOCON, a wrapper to return C_DPOCON using objects
int PSI_DPOCON(int irrep, char uplo, int n, SharedMatrix a, int lda, double anorm, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DPOCON(uplo, n, a->pointer(irrep)[0], lda, anorm, rcond->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DPOEQU, a wrapper to return C_DPOEQU using objects
int PSI_DPOEQU(int irrep, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> s, boost::shared_ptr<Vector> scond, boost::shared_ptr<Vector> amax)
{
    return C_DPOEQU(n, a->pointer(irrep)[0], lda, s->pointer(irrep), scond->pointer(irrep), amax->pointer(irrep));
}
/// PSI_DPORFS, a wrapper to return C_DPORFS using objects
int PSI_DPORFS(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DPORFS(uplo, n, nrhs, a->pointer(irrep)[0], lda, af->pointer(irrep)[0], ldaf, b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DPOSV, a wrapper to return C_DPOSV using objects
int PSI_DPOSV(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb)
{
    return C_DPOSV(uplo, n, nrhs, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb);
}
/// PSI_DPOSVX, a wrapper to return C_DPOSVX using objects
int PSI_DPOSVX(int irrep, char fact, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf, char equed, boost::shared_ptr<Vector> s, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DPOSVX(fact, uplo, n, nrhs, a->pointer(irrep)[0], lda, af->pointer(irrep)[0], ldaf, equed, s->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, rcond->pointer(irrep), ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DPOTRF, a wrapper to return C_DPOTRF using objects
int PSI_DPOTRF(int irrep, char uplo, int n, SharedMatrix a, int lda)
{
    return C_DPOTRF(uplo, n, a->pointer(irrep)[0], lda);
}
/// PSI_DPOTRI, a wrapper to return C_DPOTRI using objects
int PSI_DPOTRI(int irrep, char uplo, int n, SharedMatrix a, int lda)
{
    return C_DPOTRI(uplo, n, a->pointer(irrep)[0], lda);
}
/// PSI_DPOTRS, a wrapper to return C_DPOTRS using objects
int PSI_DPOTRS(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb)
{
    return C_DPOTRS(uplo, n, nrhs, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb);
}
/// PSI_DPTCON, a wrapper to return C_DPTCON using objects
int PSI_DPTCON(int irrep, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, double anorm, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> work)
{
    return C_DPTCON(n, d->pointer(irrep), e->pointer(irrep), anorm, rcond->pointer(irrep), work->pointer(irrep));
}
/// PSI_DPTEQR, a wrapper to return C_DPTEQR using objects
int PSI_DPTEQR(int irrep, char compz, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work)
{
    return C_DPTEQR(compz, n, d->pointer(irrep), e->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep));
}
/// PSI_DPTRFS, a wrapper to return C_DPTRFS using objects
int PSI_DPTRFS(int irrep, int n, int nrhs, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, boost::shared_ptr<Vector> df, boost::shared_ptr<Vector> ef, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work)
{
    return C_DPTRFS(n, nrhs, d->pointer(irrep), e->pointer(irrep), df->pointer(irrep), ef->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep));
}
/// PSI_DPTSV, a wrapper to return C_DPTSV using objects
int PSI_DPTSV(int irrep, int n, int nrhs, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix b, int ldb)
{
    return C_DPTSV(n, nrhs, d->pointer(irrep), e->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DPTSVX, a wrapper to return C_DPTSVX using objects
int PSI_DPTSVX(int irrep, char fact, int n, int nrhs, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, boost::shared_ptr<Vector> df, boost::shared_ptr<Vector> ef, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work)
{
    return C_DPTSVX(fact, n, nrhs, d->pointer(irrep), e->pointer(irrep), df->pointer(irrep), ef->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, rcond->pointer(irrep), ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep));
}
/// PSI_DPTTRF, a wrapper to return C_DPTTRF using objects
int PSI_DPTTRF(int irrep, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e)
{
    return C_DPTTRF(n, d->pointer(irrep), e->pointer(irrep));
}
/// PSI_DPTTRS, a wrapper to return C_DPTTRS using objects
int PSI_DPTTRS(int irrep, int n, int nrhs, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix b, int ldb)
{
    return C_DPTTRS(n, nrhs, d->pointer(irrep), e->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DSBEV, a wrapper to return C_DSBEV using objects
int PSI_DSBEV(int irrep, char jobz, char uplo, int n, int kd, SharedMatrix ab, int ldab, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work)
{
    return C_DSBEV(jobz, uplo, n, kd, ab->pointer(irrep)[0], ldab, w->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep));
}
/// PSI_DSBEVD, a wrapper to return C_DSBEVD using objects
int PSI_DSBEVD(int irrep, char jobz, char uplo, int n, int kd, SharedMatrix ab, int ldab, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DSBEVD(jobz, uplo, n, kd, ab->pointer(irrep)[0], ldab, w->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DSBEVX, a wrapper to return C_DSBEVX using objects
int PSI_DSBEVX(int irrep, char jobz, char range, char uplo, int n, int kd, SharedMatrix ab, int ldab, SharedMatrix q, int ldq, double vl, double vu, int il, int iu, double abstol, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork, boost::shared_ptr<IntVector> ifail)
{
    return C_DSBEVX(jobz, range, uplo, n, kd, ab->pointer(irrep)[0], ldab, q->pointer(irrep)[0], ldq, vl, vu, il, iu, abstol, m->pointer(irrep), w->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), iwork->pointer(irrep), ifail->pointer(irrep));
}
/// PSI_DSBGST, a wrapper to return C_DSBGST using objects
int PSI_DSBGST(int irrep, char vect, char uplo, int n, int ka, int kb, SharedMatrix ab, int ldab, SharedMatrix bb, int ldbb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> work)
{
    return C_DSBGST(vect, uplo, n, ka, kb, ab->pointer(irrep)[0], ldab, bb->pointer(irrep)[0], ldbb, x->pointer(irrep)[0], ldx, work->pointer(irrep));
}
/// PSI_DSBGV, a wrapper to return C_DSBGV using objects
int PSI_DSBGV(int irrep, char jobz, char uplo, int n, int ka, int kb, SharedMatrix ab, int ldab, SharedMatrix bb, int ldbb, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work)
{
    return C_DSBGV(jobz, uplo, n, ka, kb, ab->pointer(irrep)[0], ldab, bb->pointer(irrep)[0], ldbb, w->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep));
}
/// PSI_DSBGVD, a wrapper to return C_DSBGVD using objects
int PSI_DSBGVD(int irrep, char jobz, char uplo, int n, int ka, int kb, SharedMatrix ab, int ldab, SharedMatrix bb, int ldbb, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DSBGVD(jobz, uplo, n, ka, kb, ab->pointer(irrep)[0], ldab, bb->pointer(irrep)[0], ldbb, w->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DSBGVX, a wrapper to return C_DSBGVX using objects
int PSI_DSBGVX(int irrep, char jobz, char range, char uplo, int n, int ka, int kb, SharedMatrix ab, int ldab, SharedMatrix bb, int ldbb, SharedMatrix q, int ldq, double vl, double vu, int il, int iu, double abstol, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork, boost::shared_ptr<IntVector> ifail)
{
    return C_DSBGVX(jobz, range, uplo, n, ka, kb, ab->pointer(irrep)[0], ldab, bb->pointer(irrep)[0], ldbb, q->pointer(irrep)[0], ldq, vl, vu, il, iu, abstol, m->pointer(irrep), w->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), iwork->pointer(irrep), ifail->pointer(irrep));
}
/// PSI_DSBTRD, a wrapper to return C_DSBTRD using objects
int PSI_DSBTRD(int irrep, char vect, char uplo, int n, int kd, SharedMatrix ab, int ldab, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix q, int ldq, boost::shared_ptr<Vector> work)
{
    return C_DSBTRD(vect, uplo, n, kd, ab->pointer(irrep)[0], ldab, d->pointer(irrep), e->pointer(irrep), q->pointer(irrep)[0], ldq, work->pointer(irrep));
}
/// PSI_DSGESV, a wrapper to return C_DSGESV using objects
int PSI_DSGESV(int irrep, int n, int nrhs, boost::shared_ptr<Vector> a, int lda, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx, SharedMatrix work, boost::shared_ptr<IntVector> iter)
{
    return C_DSGESV(n, nrhs, a->pointer(irrep), lda, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, work->pointer(irrep)[0], iter->pointer(irrep));
}
/// PSI_DSTEBZ, a wrapper to return C_DSTEBZ using objects
int PSI_DSTEBZ(int irrep, char range, char order, int n, double vl, double vu, int il, int iu, double abstol, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, boost::shared_ptr<IntVector> m, boost::shared_ptr<IntVector> nsplit, boost::shared_ptr<Vector> w, boost::shared_ptr<IntVector> iblock, boost::shared_ptr<IntVector> isplit, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DSTEBZ(range, order, n, vl, vu, il, iu, abstol, d->pointer(irrep), e->pointer(irrep), m->pointer(irrep), nsplit->pointer(irrep), w->pointer(irrep), iblock->pointer(irrep), isplit->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DSTEDC, a wrapper to return C_DSTEDC using objects
int PSI_DSTEDC(int irrep, char compz, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DSTEDC(compz, n, d->pointer(irrep), e->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DSTEGR, a wrapper to return C_DSTEGR using objects
int PSI_DSTEGR(int irrep, char jobz, char range, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, double vl, double vu, int il, int iu, double abstol, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<IntVector> isuppz, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DSTEGR(jobz, range, n, d->pointer(irrep), e->pointer(irrep), vl, vu, il, iu, abstol, m->pointer(irrep), w->pointer(irrep), z->pointer(irrep)[0], ldz, isuppz->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DSTEIN, a wrapper to return C_DSTEIN using objects
int PSI_DSTEIN(int irrep, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, int m, boost::shared_ptr<Vector> w, boost::shared_ptr<IntVector> iblock, boost::shared_ptr<IntVector> isplit, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork, boost::shared_ptr<IntVector> ifail)
{
    return C_DSTEIN(n, d->pointer(irrep), e->pointer(irrep), m, w->pointer(irrep), iblock->pointer(irrep), isplit->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), iwork->pointer(irrep), ifail->pointer(irrep));
}
/// PSI_DSTEQR, a wrapper to return C_DSTEQR using objects
int PSI_DSTEQR(int irrep, char compz, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work)
{
    return C_DSTEQR(compz, n, d->pointer(irrep), e->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep));
}
/// PSI_DSTERF, a wrapper to return C_DSTERF using objects
int PSI_DSTERF(int irrep, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e)
{
    return C_DSTERF(n, d->pointer(irrep), e->pointer(irrep));
}
/// PSI_DSTEV, a wrapper to return C_DSTEV using objects
int PSI_DSTEV(int irrep, char jobz, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work)
{
    return C_DSTEV(jobz, n, d->pointer(irrep), e->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep));
}
/// PSI_DSTEVD, a wrapper to return C_DSTEVD using objects
int PSI_DSTEVD(int irrep, char jobz, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DSTEVD(jobz, n, d->pointer(irrep), e->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DSTEVR, a wrapper to return C_DSTEVR using objects
int PSI_DSTEVR(int irrep, char jobz, char range, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, double vl, double vu, int il, int iu, double abstol, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<IntVector> isuppz, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DSTEVR(jobz, range, n, d->pointer(irrep), e->pointer(irrep), vl, vu, il, iu, abstol, m->pointer(irrep), w->pointer(irrep), z->pointer(irrep)[0], ldz, isuppz->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DSTEVX, a wrapper to return C_DSTEVX using objects
int PSI_DSTEVX(int irrep, char jobz, char range, int n, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, double vl, double vu, int il, int iu, double abstol, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork, boost::shared_ptr<IntVector> ifail)
{
    return C_DSTEVX(jobz, range, n, d->pointer(irrep), e->pointer(irrep), vl, vu, il, iu, abstol, m->pointer(irrep), w->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), iwork->pointer(irrep), ifail->pointer(irrep));
}
/// PSI_DSYCON, a wrapper to return C_DSYCON using objects
int PSI_DSYCON(int irrep, char uplo, int n, SharedMatrix a, int lda, boost::shared_ptr<IntVector> ipiv, double anorm, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DSYCON(uplo, n, a->pointer(irrep)[0], lda, ipiv->pointer(irrep), anorm, rcond->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DSYEV, a wrapper to return C_DSYEV using objects
int PSI_DSYEV(int irrep, char jobz, char uplo, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> w, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DSYEV(jobz, uplo, n, a->pointer(irrep)[0], lda, w->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DSYEVD, a wrapper to return C_DSYEVD using objects
int PSI_DSYEVD(int irrep, char jobz, char uplo, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> w, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DSYEVD(jobz, uplo, n, a->pointer(irrep)[0], lda, w->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DSYEVR, a wrapper to return C_DSYEVR using objects
int PSI_DSYEVR(int irrep, char jobz, char range, char uplo, int n, SharedMatrix a, int lda, double vl, double vu, int il, int iu, double abstol, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<IntVector> isuppz, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DSYEVR(jobz, range, uplo, n, a->pointer(irrep)[0], lda, vl, vu, il, iu, abstol, m->pointer(irrep), w->pointer(irrep), z->pointer(irrep)[0], ldz, isuppz->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DSYEVX, a wrapper to return C_DSYEVX using objects
int PSI_DSYEVX(int irrep, char jobz, char range, char uplo, int n, SharedMatrix a, int lda, double vl, double vu, int il, int iu, double abstol, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, boost::shared_ptr<IntVector> ifail)
{
    return C_DSYEVX(jobz, range, uplo, n, a->pointer(irrep)[0], lda, vl, vu, il, iu, abstol, m->pointer(irrep), w->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), lwork, iwork->pointer(irrep), ifail->pointer(irrep));
}
/// PSI_DSYGST, a wrapper to return C_DSYGST using objects
int PSI_DSYGST(int irrep, int itype, char uplo, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb)
{
    return C_DSYGST(itype, uplo, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb);
}
/// PSI_DSYGV, a wrapper to return C_DSYGV using objects
int PSI_DSYGV(int irrep, int itype, char jobz, char uplo, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> w, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DSYGV(itype, jobz, uplo, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, w->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DSYGVD, a wrapper to return C_DSYGVD using objects
int PSI_DSYGVD(int irrep, int itype, char jobz, char uplo, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> w, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DSYGVD(itype, jobz, uplo, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, w->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DSYGVX, a wrapper to return C_DSYGVX using objects
int PSI_DSYGVX(int irrep, int itype, char jobz, char range, char uplo, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, double vl, double vu, int il, int iu, double abstol, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> w, SharedMatrix z, int ldz, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, boost::shared_ptr<IntVector> ifail)
{
    return C_DSYGVX(itype, jobz, range, uplo, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, vl, vu, il, iu, abstol, m->pointer(irrep), w->pointer(irrep), z->pointer(irrep)[0], ldz, work->pointer(irrep), lwork, iwork->pointer(irrep), ifail->pointer(irrep));
}
/// PSI_DSYRFS, a wrapper to return C_DSYRFS using objects
int PSI_DSYRFS(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DSYRFS(uplo, n, nrhs, a->pointer(irrep)[0], lda, af->pointer(irrep)[0], ldaf, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DSYSV, a wrapper to return C_DSYSV using objects
int PSI_DSYSV(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DSYSV(uplo, n, nrhs, a->pointer(irrep)[0], lda, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb, work->pointer(irrep), lwork);
}
/// PSI_DSYSVX, a wrapper to return C_DSYSVX using objects
int PSI_DSYSVX(int irrep, char fact, char uplo, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix af, int ldaf, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> rcond)
{
    return C_DSYSVX(fact, uplo, n, nrhs, a->pointer(irrep)[0], lda, af->pointer(irrep)[0], ldaf, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, rcond->pointer(irrep));
}
/// PSI_DSYTRD, a wrapper to return C_DSYTRD using objects
int PSI_DSYTRD(int irrep, char uplo, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> d, boost::shared_ptr<Vector> e, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DSYTRD(uplo, n, a->pointer(irrep)[0], lda, d->pointer(irrep), e->pointer(irrep), tau->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DSYTRF, a wrapper to return C_DSYTRF using objects
int PSI_DSYTRF(int irrep, char uplo, int n, SharedMatrix a, int lda, boost::shared_ptr<IntVector> ipiv, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DSYTRF(uplo, n, a->pointer(irrep)[0], lda, ipiv->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DSYTRI, a wrapper to return C_DSYTRI using objects
int PSI_DSYTRI(int irrep, char uplo, int n, SharedMatrix a, int lda, boost::shared_ptr<IntVector> ipiv, boost::shared_ptr<Vector> work)
{
    return C_DSYTRI(uplo, n, a->pointer(irrep)[0], lda, ipiv->pointer(irrep), work->pointer(irrep));
}
/// PSI_DSYTRS, a wrapper to return C_DSYTRS using objects
int PSI_DSYTRS(int irrep, char uplo, int n, int nrhs, SharedMatrix a, int lda, boost::shared_ptr<IntVector> ipiv, SharedMatrix b, int ldb)
{
    return C_DSYTRS(uplo, n, nrhs, a->pointer(irrep)[0], lda, ipiv->pointer(irrep), b->pointer(irrep)[0], ldb);
}
/// PSI_DTBCON, a wrapper to return C_DTBCON using objects
int PSI_DTBCON(int irrep, char norm, char uplo, char diag, int n, int kd, SharedMatrix ab, int ldab, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DTBCON(norm, uplo, diag, n, kd, ab->pointer(irrep)[0], ldab, rcond->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DTBRFS, a wrapper to return C_DTBRFS using objects
int PSI_DTBRFS(int irrep, char uplo, char trans, char diag, int n, int kd, int nrhs, SharedMatrix ab, int ldab, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DTBRFS(uplo, trans, diag, n, kd, nrhs, ab->pointer(irrep)[0], ldab, b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DTBTRS, a wrapper to return C_DTBTRS using objects
int PSI_DTBTRS(int irrep, char uplo, char trans, char diag, int n, int kd, int nrhs, SharedMatrix ab, int ldab, SharedMatrix b, int ldb)
{
    return C_DTBTRS(uplo, trans, diag, n, kd, nrhs, ab->pointer(irrep)[0], ldab, b->pointer(irrep)[0], ldb);
}
/// PSI_DTGEVC, a wrapper to return C_DTGEVC using objects
int PSI_DTGEVC(int irrep, char side, char howmny, int n, SharedMatrix s, int lds, SharedMatrix p, int ldp, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, int mm, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> work)
{
    return C_DTGEVC(side, howmny, n, s->pointer(irrep)[0], lds, p->pointer(irrep)[0], ldp, vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr, mm, m->pointer(irrep), work->pointer(irrep));
}
/// PSI_DTGEXC, a wrapper to return C_DTGEXC using objects
int PSI_DTGEXC(int irrep, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, SharedMatrix q, int ldq, SharedMatrix z, int ldz, boost::shared_ptr<IntVector> ifst, boost::shared_ptr<IntVector> ilst, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DTGEXC(n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, q->pointer(irrep)[0], ldq, z->pointer(irrep)[0], ldz, ifst->pointer(irrep), ilst->pointer(irrep), work->pointer(irrep), lwork);
}
/// PSI_DTGSEN, a wrapper to return C_DTGSEN using objects
int PSI_DTGSEN(int irrep, int ijob, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, boost::shared_ptr<Vector> alphar, boost::shared_ptr<Vector> alphai, boost::shared_ptr<Vector> beta, SharedMatrix q, int ldq, SharedMatrix z, int ldz, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> pl, boost::shared_ptr<Vector> pr, boost::shared_ptr<Vector> dif, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DTGSEN(ijob, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, alphar->pointer(irrep), alphai->pointer(irrep), beta->pointer(irrep), q->pointer(irrep)[0], ldq, z->pointer(irrep)[0], ldz, m->pointer(irrep), pl->pointer(irrep), pr->pointer(irrep), dif->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DTGSJA, a wrapper to return C_DTGSJA using objects
int PSI_DTGSJA(int irrep, char jobu, char jobv, char jobq, int m, int p, int n, int k, int l, SharedMatrix a, int lda, SharedMatrix b, int ldb, double tola, double tolb, boost::shared_ptr<Vector> alpha, boost::shared_ptr<Vector> beta, SharedMatrix u, int ldu, SharedMatrix v, int ldv, SharedMatrix q, int ldq, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> ncycle)
{
    return C_DTGSJA(jobu, jobv, jobq, m, p, n, k, l, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, tola, tolb, alpha->pointer(irrep), beta->pointer(irrep), u->pointer(irrep)[0], ldu, v->pointer(irrep)[0], ldv, q->pointer(irrep)[0], ldq, work->pointer(irrep), ncycle->pointer(irrep));
}
/// PSI_DTGSNA, a wrapper to return C_DTGSNA using objects
int PSI_DTGSNA(int irrep, char job, char howmny, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, boost::shared_ptr<Vector> s, boost::shared_ptr<Vector> dif, int mm, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork)
{
    return C_DTGSNA(job, howmny, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr, s->pointer(irrep), dif->pointer(irrep), mm, m->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep));
}
/// PSI_DTGSYL, a wrapper to return C_DTGSYL using objects
int PSI_DTGSYL(int irrep, char trans, int ijob, int m, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, SharedMatrix c, int ldc, SharedMatrix d, int ldd, SharedMatrix e, int lde, SharedMatrix f, int ldf, boost::shared_ptr<Vector> dif, boost::shared_ptr<Vector> scale, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork)
{
    return C_DTGSYL(trans, ijob, m, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, c->pointer(irrep)[0], ldc, d->pointer(irrep)[0], ldd, e->pointer(irrep)[0], lde, f->pointer(irrep)[0], ldf, dif->pointer(irrep), scale->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep));
}
/// PSI_DTRCON, a wrapper to return C_DTRCON using objects
int PSI_DTRCON(int irrep, char norm, char uplo, char diag, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> rcond, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DTRCON(norm, uplo, diag, n, a->pointer(irrep)[0], lda, rcond->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DTREVC, a wrapper to return C_DTREVC using objects
int PSI_DTREVC(int irrep, char side, char howmny, int n, SharedMatrix t, int ldt, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, int mm, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> work)
{
    return C_DTREVC(side, howmny, n, t->pointer(irrep)[0], ldt, vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr, mm, m->pointer(irrep), work->pointer(irrep));
}
/// PSI_DTREXC, a wrapper to return C_DTREXC using objects
int PSI_DTREXC(int irrep, char compq, int n, SharedMatrix t, int ldt, SharedMatrix q, int ldq, boost::shared_ptr<IntVector> ifst, boost::shared_ptr<IntVector> ilst, boost::shared_ptr<Vector> work)
{
    return C_DTREXC(compq, n, t->pointer(irrep)[0], ldt, q->pointer(irrep)[0], ldq, ifst->pointer(irrep), ilst->pointer(irrep), work->pointer(irrep));
}
/// PSI_DTRRFS, a wrapper to return C_DTRRFS using objects
int PSI_DTRRFS(int irrep, char uplo, char trans, char diag, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb, SharedMatrix x, int ldx, boost::shared_ptr<Vector> ferr, boost::shared_ptr<Vector> berr, boost::shared_ptr<Vector> work, boost::shared_ptr<IntVector> iwork)
{
    return C_DTRRFS(uplo, trans, diag, n, nrhs, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, x->pointer(irrep)[0], ldx, ferr->pointer(irrep), berr->pointer(irrep), work->pointer(irrep), iwork->pointer(irrep));
}
/// PSI_DTRSEN, a wrapper to return C_DTRSEN using objects
int PSI_DTRSEN(int irrep, char job, char compq, int n, SharedMatrix t, int ldt, SharedMatrix q, int ldq, boost::shared_ptr<Vector> wr, boost::shared_ptr<Vector> wi, boost::shared_ptr<IntVector> m, boost::shared_ptr<Vector> s, boost::shared_ptr<Vector> sep, boost::shared_ptr<Vector> work, int lwork, boost::shared_ptr<IntVector> iwork, int liwork)
{
    return C_DTRSEN(job, compq, n, t->pointer(irrep)[0], ldt, q->pointer(irrep)[0], ldq, wr->pointer(irrep), wi->pointer(irrep), m->pointer(irrep), s->pointer(irrep), sep->pointer(irrep), work->pointer(irrep), lwork, iwork->pointer(irrep), liwork);
}
/// PSI_DTRSNA, a wrapper to return C_DTRSNA using objects
int PSI_DTRSNA(int irrep, char job, char howmny, int n, SharedMatrix t, int ldt, SharedMatrix vl, int ldvl, SharedMatrix vr, int ldvr, boost::shared_ptr<Vector> s, boost::shared_ptr<Vector> sep, int mm, boost::shared_ptr<IntVector> m, SharedMatrix work, int ldwork, boost::shared_ptr<IntVector> iwork)
{
    return C_DTRSNA(job, howmny, n, t->pointer(irrep)[0], ldt, vl->pointer(irrep)[0], ldvl, vr->pointer(irrep)[0], ldvr, s->pointer(irrep), sep->pointer(irrep), mm, m->pointer(irrep), work->pointer(irrep)[0], ldwork, iwork->pointer(irrep));
}
/// PSI_DTRSYL, a wrapper to return C_DTRSYL using objects
int PSI_DTRSYL(int irrep, char trana, char tranb, int isgn, int m, int n, SharedMatrix a, int lda, SharedMatrix b, int ldb, SharedMatrix c, int ldc, boost::shared_ptr<Vector> scale)
{
    return C_DTRSYL(trana, tranb, isgn, m, n, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb, c->pointer(irrep)[0], ldc, scale->pointer(irrep));
}
/// PSI_DTRTRI, a wrapper to return C_DTRTRI using objects
int PSI_DTRTRI(int irrep, char uplo, char diag, int n, SharedMatrix a, int lda)
{
    return C_DTRTRI(uplo, diag, n, a->pointer(irrep)[0], lda);
}
/// PSI_DTRTRS, a wrapper to return C_DTRTRS using objects
int PSI_DTRTRS(int irrep, char uplo, char trans, char diag, int n, int nrhs, SharedMatrix a, int lda, SharedMatrix b, int ldb)
{
    return C_DTRTRS(uplo, trans, diag, n, nrhs, a->pointer(irrep)[0], lda, b->pointer(irrep)[0], ldb);
}
/// PSI_DTZRQF, a wrapper to return C_DTZRQF using objects
int PSI_DTZRQF(int irrep, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau)
{
    return C_DTZRQF(m, n, a->pointer(irrep)[0], lda, tau->pointer(irrep));
}
/// PSI_DTZRZF, a wrapper to return C_DTZRZF using objects
int PSI_DTZRZF(int irrep, int m, int n, SharedMatrix a, int lda, boost::shared_ptr<Vector> tau, boost::shared_ptr<Vector> work, int lwork)
{
    return C_DTZRZF(m, n, a->pointer(irrep)[0], lda, tau->pointer(irrep), work->pointer(irrep), lwork);
}
}

