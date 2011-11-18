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
 * dnrm2
 */
double F_DNRM2(integer n,doublereal*x,integer incx){
    return DNRM2(n,x,incx);
}
/**
 * ddot
 */
double F_DDOT(integer n,doublereal*dx,integer incx,doublereal*dy,integer incy){
    return DDOT(n,dx,incx,dy,incy);
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

