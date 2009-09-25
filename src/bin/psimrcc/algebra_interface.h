#ifndef _psi_src_bin_psimrcc_algebra_interface_h_
#define _psi_src_bin_psimrcc_algebra_interface_h_

#ifndef FC_SYMBOL
#define FC_SYMBOL 2
#endif

#if FC_SYMBOL==2
#define F_DAXPY daxpy_
#define F_DCOPY dcopy_
#define F_DGEMM dgemm_
#define F_DROT drot_
#define F_DSCAL dscal_
#define F_DGEMV dgemv_
#define F_DSPMV dspmv_
#define F_DDOT  ddot_
#elif FC_SYMBOL==1
#define F_DAXPY daxpy
#define F_DCOPY dcopy
#define F_DGEMM dgemm
#define F_DROT drot
#define F_DSCAL dscal
#define F_DGEMV dgemv
#define F_DSPMV dspmv
#define F_DDOT  ddot
#elif FC_SYMBOL==3
#define F_DAXPY DAXPY
#define F_DCOPY DCOPY
#define F_DGEMM DGEMM
#define F_DROT DROT
#define F_DSCAL DSCAL
#define F_DGEMV DGEMV
#define F_DSPMV DSPMV
#define F_DDOT  DDOT
#elif FC_SYMBOL==4
#define F_DAXPY DAXPY_
#define F_DCOPY DCOPY_
#define F_DGEMM DGEMM_
#define F_DROT DROT_
#define F_DSCAL DSCAL_
#define F_DGEMV DGEMV_
#define F_DSPMV DSPMV_
#define F_DDOT  DDOT_
#endif

namespace psi{ namespace psimrcc{

extern "C" void F_DAXPY(int *length, double *a, double *x, int *inc_x,
                    double *y, int *inc_y);
extern "C" void F_DCOPY(int *length, double *x, int *inc_x,
                    double *y, int *inc_y);
extern "C" void F_DGEMM(const char *transa, const char *transb, int *m, int *n, int *k,
                    double *alpha, double *A, int *lda, double *B, int *ldb,
                    double *beta, double *C, int *ldc);
extern "C" void F_DROT(int *ntot,double *x, int *incx,double *y, int *incy,
                  double *cotheta,double *sintheta);
extern "C" void F_DSCAL(int *n, double *alpha, double *vec, int *inc);
extern "C" void F_DGEMV(char *transa, int *m, int *n, double *alpha, double *A,
                    int *lda, double *X, int *inc_x, double *beta,
                    double *Y, int *inc_y);
extern "C" double F_DDOT(int *n, double *x, int *incx, double *y, int *incy);

void C_DGEMM_12(int m, int n, int k, double alpha,double *A, int nra,
                double *B, int ncb, double beta, double *C, int ncc);
void C_DGEMM_22(int m, int n, int k, double alpha,double *A, int nca,
                double *B, int ncb, double beta, double *C, int ncc);

// void C_DGEMM_11(int m, int n, int k, double alpha,double *A, int nca,
//                 double *B, int ncb, double beta, double *C, int ncc);
// void C_DGEMM_21(int m, int n, int k, double alpha,double *A, int nca,
//                 double *B, int ncb, double beta, double *C, int ncc);


#if FC_SYMBOL==2
#define F_DGEEV dgeev_
#define F_DGESV dgesv_
#define F_DGETRF dgetrf_
#define F_DGETRI dgetri_
#define F_DGESVD dgesvd_
#define F_DSYEV dsyev_
#elif FC_SYMBOL==1
#define F_DGEEV dgeev
#define F_DGESV dgesv
#define F_DGETRF dgetrf
#define F_DGETRI dgetri
#define F_DGESVD dgesvd
#define F_DSYEV dsyev
#elif FC_SYMBOL==3
#define F_DGEEV DGEEV
#define F_DGESV DGESV
#define F_DGETRF DGETRF
#define F_DGETRI DGETRI
#define F_DGESVD DGESVD
#define F_DSYEV DSYEV
#elif FC_SYMBOL==4
#define F_DGEEV DGEEV_
#define F_DGESV DGESV_
#define F_DGETRF DGETRF_
#define F_DGETRI DGETRI_
#define F_DGESVD DGESVD_
#define F_DSYEV DSYEV_
#endif

extern "C" void F_DGEEV(const char *jobvl, const char *jobvr, int *n, double *a, int *lda,
                    double *wr, double *wi, double *vl, int *ldvl, double *vr,
                    int *ldvr, double *work, int *lwork, int *info);
extern "C" void F_DGESV(int *n, int *nrhs, double *A, int *lda, int *ipiv,
                    double *B, int *ldb, int *info);

extern "C" void F_DSYEV(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO );

}}

#endif // _psi_src_bin_psimrcc_algebra_interface_h_
