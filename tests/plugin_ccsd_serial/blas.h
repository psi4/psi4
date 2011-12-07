#ifndef BLAS_H
#define BLAS_H

/**
 * fortran-ordered blas routines
 */

#define F77NAME(x) x##_
typedef long int integer;
typedef double doublereal;


/**
 * fortran-ordered dgemv
 */
void F_DGEMV(char trans,integer m,integer n,doublereal alpha,doublereal*A,integer lda,
            doublereal*X,integer incx,doublereal beta,doublereal*Y,integer incy);
/**
 * fortran-ordered dgemm
 */
void F_DGEMM(char transa,char transb, integer m, integer n, integer k,
            doublereal alpha,doublereal*A,integer lda,doublereal*B,integer ldb,
            doublereal beta,doublereal*C,integer ldc);
/**
 * ddot
 */
double F_DDOT(integer n,doublereal*dx,integer incx,doublereal*dy,integer incy);

/**
 * dnrm2
 */
double F_DNRM2(integer n,doublereal*x,integer incx);

/**
 * dcopy
 */
void F_DCOPY(integer n,doublereal*dx,integer incx,doublereal*dy,integer incy);

/**
 * daxpy
 */
void F_DAXPY(integer n,doublereal da,doublereal*dx,integer incx,doublereal*dy,
             integer incy);

/**
 * name manging for fortran-ordered dgemv
 */
extern "C" {
    void F77NAME(dgemv)(char&trans,integer&m,integer&n,doublereal&alpha,doublereal*A,integer&lda,
            doublereal*X,integer&incx,doublereal&beta,doublereal*Y,integer&incy);
};
inline void DGEMV(char&trans,integer&m,integer&n,doublereal&alpha,doublereal*A,integer&lda,
            doublereal*X,integer&incx,doublereal&beta,doublereal*Y,integer&incy){
    F77NAME(dgemv)(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy);
}
/**
 * name manging for fortran-ordered dgemm
 */
extern "C" {
    void F77NAME(dgemm)(char&transa,char&transb,integer&m,integer&n,integer&k,
         doublereal&alpha,doublereal*A,integer&lda,doublereal*B,integer&ldb,
         doublereal&beta,doublereal*C,integer&ldc);
};
inline void DGEMM(char&transa,char&transb,integer&m,integer&n,integer&k,
         doublereal&alpha,doublereal*A,integer&lda,doublereal*B,integer&ldb,
         doublereal&beta,doublereal*C,integer&ldc)
{
    F77NAME(dgemm)(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
};
/**
 * name manging dcopy
 */
extern "C" {
    void F77NAME(dcopy)(integer&n,doublereal*dx,integer&incx,doublereal*dy,
         integer&incy);
};
inline void DCOPY(integer&n,doublereal*dx,integer&incx,doublereal*dy,
            integer&incy){
    F77NAME(dcopy)(n,dx,incx,dy,incy);
}
/**
 * name manging daxpy
 */
extern "C" {
   void F77NAME(daxpy)(integer&n,doublereal&da,doublereal*dx,integer&incx,
        doublereal*dy,integer&incy);
};
inline void DAXPY(integer&n,doublereal&da,doublereal*dx,integer&incx,
            doublereal*dy,integer&incy)
{
    F77NAME(daxpy)(n,da,dx,incx,dy,incy);
};
/**
 * name manging dnrm2
 */
extern"C"{
    double F77NAME(dnrm2)(integer&N,doublereal*X,integer&INCX);
};
inline double DNRM2(integer&N,doublereal*X,integer&INCX){
    return F77NAME(dnrm2)(N,X,INCX);
};
/**
 * name manging dgesv
 */
extern"C" {
    void F77NAME(dgesv)(integer &N,integer &NRHS,doublereal*A,integer &LDA,integer*IPIV,doublereal*B,integer &LDB,integer &INFO);
};
inline void DGESV(integer &N,integer &NRHS,doublereal*A,integer &LDA,integer*IPIV,doublereal*B,integer &LDB,integer &INFO){
    F77NAME(dgesv)(N,NRHS,A,LDA,IPIV,B,LDB,INFO);
};
/**
 * name manging ddot
 */
extern "C" {
    double F77NAME(ddot)(integer&n,doublereal*dx,integer&incx,doublereal*dy,integer&incy);
};
inline double DDOT(integer&n,doublereal*dx,integer&incx,doublereal*dy,integer&incy){
    return F77NAME(ddot)(n,dx,incx,dy,incy);
}

#endif
