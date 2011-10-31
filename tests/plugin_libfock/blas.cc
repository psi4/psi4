#include"blas.h"
/**
 * fortran-ordered dgemm
 */
void F_DGEMM(char transa,char transb, integer m, integer n, integer k,
            doublereal alpha,doublereal*A,integer lda,doublereal*B,integer ldb,
            doublereal beta,doublereal*C,integer ldc){
    DGEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
}
/**
 * dcopy
 */
void F_DCOPY(integer n,doublereal*dx,integer incx,doublereal*dy,integer incy){
    DCOPY(n,dx,incx,dy,incy);
}

/**
 * daxpy
 */
void F_DAXPY(integer n,doublereal da,doublereal*dx,integer incx,doublereal*dy,
             integer incy){
    DAXPY(n,da,dx,incx,dy,incy);
}

