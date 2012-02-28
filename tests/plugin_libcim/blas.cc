#include"blas.h"
#include<stdlib.h>
#include<stdio.h>
/**
 * fortran-ordered dgemv
 */
void F_DGEMV(char trans,integer m,integer n,doublereal alpha,doublereal*A,integer lda,
            doublereal*X,integer incx,doublereal beta,doublereal*Y,integer incy){
    DGEMV(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy);
}
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

/**
 *  Diagonalize a real symmetric matrix
 */
void Diagonalize(integer N,doublereal*A,doublereal*W){
  char JOBZ = 'V';
  char UPLO = 'U';
  integer LDA = N;
  integer LWORK = 3*N-1;
  doublereal*WORK=(doublereal*)malloc(LWORK*sizeof(doublereal)); 
  integer INFO=0;
  DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO);
}
void Diagonalize2(integer N,doublereal*AP,doublereal*W,doublereal*Z){
  char JOBZ = 'V';
  char UPLO = 'U';
  integer LDZ = N;
  doublereal*WORK=(doublereal*)malloc(3*N*sizeof(doublereal)); 
  integer INFO=0;
  DSPEV(JOBZ,UPLO,N,AP,W,Z,LDZ,WORK,INFO);
}


