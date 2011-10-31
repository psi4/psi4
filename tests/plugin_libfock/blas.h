#ifndef BLAS_H
#define BLAS_H

/**
 * fortran-ordered blas routines
 */

#define F77NAME(x) x##_
typedef long int integer;
typedef double doublereal;


/**
 * fortran-ordered dgemm
 */
void F_DGEMM(char transa,char transb, integer m, integer n, integer k,
            doublereal alpha,doublereal*A,integer lda,doublereal*B,integer ldb,
            doublereal beta,doublereal*C,integer ldc);
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

#endif
