#ifndef GALINALG_H_
#define GALINALG_H_

#define sgemm_  F77_FUNC(sgemm, SGEMM)
#define dgemm_  F77_FUNC(dgemm, DGEMM)
#define zgemm_  F77_FUNC(zgemm, ZGEMM)
#define cgemm_  F77_FUNC(cgemm, CGEMM)
#define dgetrf_ F77_FUNC(dgetrf,DGETRF)
#define dgetrs_ F77_FUNC(dgetrs,DGETRS)
#define dtrsm_  F77_FUNC(dtrsm, DTRSM)
#define gal_sgemm_  F77_FUNC_(gal_sgemm, GAL_SGEMM)
#define gal_dgemm_  F77_FUNC_(gal_dgemm, GAL_DGEMM)
#define gal_zgemm_  F77_FUNC_(gal_zgemm, GAL_ZGEMM)
#define gal_cgemm_  F77_FUNC_(gal_cgemm, GAL_CGEMM)
#define gal_dgetrf_ F77_FUNC_(gal_dgetrf,GAL_DGETRF)
#define gal_dgetrs_ F77_FUNC_(gal_dgetrs,GAL_DGETRS)
#define gal_dtrsm_  F77_FUNC_(gal_dtrsm, GAL_DTRSM)

#if HAVE_BLAS || ENABLE_F77
#   if BLAS_SIZE == SIZEOF_SHORT
#      define BlasInt short
#   elif BLAS_SIZE == SIZEOF_INT
#      define BlasInt int
#   elif BLAS_SIZE == SIZEOF_LONG
#      define BlasInt long
#   elif BLAS_SIZE == SIZEOF_LONG_LONG
#      define BlasInt long long
#   endif
#else
#   define BlasInt int
#endif

#if HAVE_BLAS
#   if defined(F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS)
extern void sgemm_(char *TRANSA, char *TRANSB, BlasInt *M, BlasInt *N, BlasInt *K, Real *ALPHA, Real *A, BlasInt *LDA, Real *B, BlasInt *LDB, Real *BETA, Real *C, BlasInt *LDC, int alen, int blen);
extern void dgemm_(char *TRANSA, char *TRANSB, BlasInt *M, BlasInt *N, BlasInt *K, DoublePrecision *ALPHA, DoublePrecision *A, BlasInt *LDA, DoublePrecision *B, BlasInt *LDB, DoublePrecision *BETA, DoublePrecision *C, BlasInt *LDC, int alen, int blen);
extern void zgemm_(char *TRANSA, char *TRANSB, BlasInt *M, BlasInt *N, BlasInt *K, DoubleComplex *ALPHA, DoubleComplex *A, BlasInt *LDA, DoubleComplex *B, BlasInt *LDB, DoubleComplex *BETA, DoubleComplex *C, BlasInt *LDC, int alen, int blen);
extern void cgemm_(char *TRANSA, char *TRANSB, BlasInt *M, BlasInt *N, BlasInt *K, SingleComplex *ALPHA, SingleComplex *A, BlasInt *LDA, SingleComplex *B, BlasInt *LDB, SingleComplex *BETA, SingleComplex *C, BlasInt *LDC, int alen, int blen);
#   else
extern void sgemm_(char *TRANSA, int alen, char *TRANSB, int blen, BlasInt *M, BlasInt *N, BlasInt *K, Real *ALPHA, Real *A, BlasInt *LDA, Real *B, BlasInt *LDB, Real *BETA, Real *C, BlasInt *LDC);
extern void dgemm_(char *TRANSA, int alen, char *TRANSB, int blen, BlasInt *M, BlasInt *N, BlasInt *K, DoublePrecision *ALPHA, DoublePrecision *A, BlasInt *LDA, DoublePrecision *B, BlasInt *LDB, DoublePrecision *BETA, DoublePrecision *C, BlasInt *LDC);
extern void zgemm_(char *TRANSA, int alen, char *TRANSB, int blen, BlasInt *M, BlasInt *N, BlasInt *K, DoubleComplex *ALPHA, DoubleComplex *A, BlasInt *LDA, DoubleComplex *B, BlasInt *LDB, DoubleComplex *BETA, DoubleComplex *C, BlasInt *LDC);
extern void cgemm_(char *TRANSA, int alen, char *TRANSB, int blen, BlasInt *M, BlasInt *N, BlasInt *K, SingleComplex *ALPHA, SingleComplex *A, BlasInt *LDA, SingleComplex *B, BlasInt *LDB, SingleComplex *BETA, SingleComplex *C, BlasInt *LDC);
#   endif
#elif ENABLE_F77
#   if defined(F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS)
extern void gal_sgemm_(char *TRANSA, char *TRANSB, BlasInt *M, BlasInt *N, BlasInt *K, Real *ALPHA, Real *A, BlasInt *LDA, Real *B, BlasInt *LDB, Real *BETA, Real *C, BlasInt *LDC, int alen, int blen);
extern void gal_dgemm_(char *TRANSA, char *TRANSB, BlasInt *M, BlasInt *N, BlasInt *K, DoublePrecision *ALPHA, DoublePrecision *A, BlasInt *LDA, DoublePrecision *B, BlasInt *LDB, DoublePrecision *BETA, DoublePrecision *C, BlasInt *LDC, int alen, int blen);
extern void gal_zgemm_(char *TRANSA, char *TRANSB, BlasInt *M, BlasInt *N, BlasInt *K, DoubleComplex *ALPHA, DoubleComplex *A, BlasInt *LDA, DoubleComplex *B, BlasInt *LDB, DoubleComplex *BETA, DoubleComplex *C, BlasInt *LDC, int alen, int blen);
extern void gal_cgemm_(char *TRANSA, char *TRANSB, BlasInt *M, BlasInt *N, BlasInt *K, SingleComplex *ALPHA, SingleComplex *A, BlasInt *LDA, SingleComplex *B, BlasInt *LDB, SingleComplex *BETA, SingleComplex *C, BlasInt *LDC, int alen, int blen);
#   else
extern void gal_sgemm_(char *TRANSA, int alen, char *TRANSB, int blen, BlasInt *M, BlasInt *N, BlasInt *K, Real *ALPHA, Real *A, BlasInt *LDA, Real *B, BlasInt *LDB, Real *BETA, Real *C, BlasInt *LDC);
extern void gal_dgemm_(char *TRANSA, int alen, char *TRANSB, int blen, BlasInt *M, BlasInt *N, BlasInt *K, DoublePrecision *ALPHA, DoublePrecision *A, BlasInt *LDA, DoublePrecision *B, BlasInt *LDB, DoublePrecision *BETA, DoublePrecision *C, BlasInt *LDC);
extern void gal_zgemm_(char *TRANSA, int alen, char *TRANSB, int blen, BlasInt *M, BlasInt *N, BlasInt *K, DoubleComplex *ALPHA, DoubleComplex *A, BlasInt *LDA, DoubleComplex *B, BlasInt *LDB, DoubleComplex *BETA, DoubleComplex *C, BlasInt *LDC);
extern void gal_cgemm_(char *TRANSA, int alen, char *TRANSB, int blen, BlasInt *M, BlasInt *N, BlasInt *K, SingleComplex *ALPHA, SingleComplex *A, BlasInt *LDA, SingleComplex *B, BlasInt *LDB, SingleComplex *BETA, SingleComplex *C, BlasInt *LDC);
#   endif
#else
#   include "xgemm.h"
#endif /* HAVE_BLAS */

#if HAVE_LAPACK
#   if defined(F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS)
extern void dtrsm_(char *side, char *uplo, char *transa, char *diag, BlasInt *m, BlasInt *n, DoublePrecision *alpha, DoublePrecision *a, BlasInt *lda, DoublePrecision *b, BlasInt *ldb, int, int, int, int );
extern void dgetrs_(char *trans, BlasInt *n, BlasInt *nrhs, DoublePrecision *a, Integer *lda, Integer *ipiv, DoublePrecision *b, Integer *ldb, Integer *info, int len );
#   else
extern void dtrsm_(char *side, int, char *uplo, int, char *transa, int, char *diag, int, BlasInt *m, BlasInt *n, DoublePrecision *alpha, DoublePrecision *a, BlasInt *lda, DoublePrecision *b, BlasInt *ldb );
extern void dgetrs_(char *trans, int len, BlasInt *n, BlasInt *nrhs, DoublePrecision *a, Integer *lda, Integer *ipiv, DoublePrecision *b, Integer *ldb, Integer *info );
#   endif
extern void dgetrf_( BlasInt *m, BlasInt *n, DoublePrecision *a, BlasInt *ld, BlasInt *ipiv, BlasInt *info );
#elif ENABLE_F77
#   if defined(F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS)
extern void gal_dtrsm_(char *side, char *uplo, char *transa, char *diag, BlasInt *m, BlasInt *n, DoublePrecision *alpha, DoublePrecision *a, BlasInt *lda, DoublePrecision *b, BlasInt *ldb, int, int, int, int );
extern void gal_dgetrs_(char *trans, BlasInt *n, BlasInt *nrhs, DoublePrecision *a, Integer *lda, Integer *ipiv, DoublePrecision *b, Integer *ldb, Integer *info, int len );
#   else
extern void gal_dtrsm_(char *side, int, char *uplo, int, char *transa, int, char *diag, int, BlasInt *m, BlasInt *n, DoublePrecision *alpha, DoublePrecision *a, BlasInt *lda, DoublePrecision *b, BlasInt *ldb );
extern void gal_dgetrs_(char *trans, int len, BlasInt *n, BlasInt *nrhs, DoublePrecision *a, Integer *lda, Integer *ipiv, DoublePrecision *b, Integer *ldb, Integer *info );
#   endif
extern void gal_dgetrf_( BlasInt *m, BlasInt *n, DoublePrecision *a, BlasInt *ld, BlasInt *ipiv, BlasInt *info );
#else
#endif /* HAVE_LAPACK */

#if HAVE_BLAS
#   if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
#       define BLAS_SGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        sgemm_(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1)
#       define BLAS_DGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        dgemm_(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1)
#       define BLAS_ZGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        zgemm_(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1)
#       define BLAS_CGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        cgemm_(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1)
#   else
#       define BLAS_SGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        sgemm_(ta, 1, tb, 1, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#       define BLAS_DGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        dgemm_(ta, 1, tb, 1, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#       define BLAS_ZGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        zgemm_(ta, 1, tb, 1, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#       define BLAS_CGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        cgemm_(ta, 1, tb, 1, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#   endif
#elif ENABLE_F77
#   if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
#       define BLAS_SGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        gal_sgemm_(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1)
#       define BLAS_DGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        gal_dgemm_(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1)
#       define BLAS_ZGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        gal_zgemm_(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1)
#       define BLAS_CGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        gal_cgemm_(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1)
#   else
#       define BLAS_SGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        gal_sgemm_(ta, 1, tb, 1, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#       define BLAS_DGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        gal_dgemm_(ta, 1, tb, 1, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#       define BLAS_ZGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        gal_zgemm_(ta, 1, tb, 1, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#       define BLAS_CGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
        gal_cgemm_(ta, 1, tb, 1, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#   endif
#else
#   define BLAS_SGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
    xb_sgemm(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#   define BLAS_DGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
    xb_dgemm(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#   define BLAS_ZGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
    xb_zgemm(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#   define BLAS_CGEMM(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc) \
    xb_cgemm(ta, tb, i, j, k, alpha, a, lda, b, ldb, beta, c, ldc)
#endif

#if HAVE_LAPACK
#   if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
#       define LAPACK_DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) \
        dtrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, 1, 1, 1, 1)
#       define LAPACK_DGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info) \
        dgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info, 1)
#   else
#       define LAPACK_DTRSM(side, 1, uplo, 1, transa, 1, diag, 1, m, n, alpha, a, lda, b, ldb) \
        dtrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
#       define LAPACK_DGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info) \
        dgetrs_(trans, 1, n, nrhs, a, lda, ipiv, b, ldb, info)
#   endif
#   define LAPACK_DGETRF dgetrf_
#elif ENABLE_F77
#   if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
#       define LAPACK_DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) \
        gal_dtrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, 1, 1, 1, 1)
#       define LAPACK_DGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info) \
        gal_dgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info, 1)
#   else
#       define LAPACK_DTRSM(side, 1, uplo, 1, transa, 1, diag, 1, m, n, alpha, a, lda, b, ldb) \
        gal_dtrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
#       define LAPACK_DGETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info) \
        gal_dgetrs_(trans, 1, n, nrhs, a, lda, ipiv, b, ldb, info)
#   endif
#   define LAPACK_DGETRF gal_dgetrf_
#elif NOFORT
#endif

#endif /* GALINALG_H_ */
