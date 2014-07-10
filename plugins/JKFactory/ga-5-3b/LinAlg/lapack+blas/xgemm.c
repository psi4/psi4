#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "xgemm.h"

#define COLMAJOR 'c'
#define ROWMAJOR 'r'

#define MDEBUG 0

void
xb_sgemm (char *transa, char *transb, int *M, int *N, int *K,
	  float *alpha, const float *a, int *p_lda, const float *b, 
	  int *p_ldb, float *beta, float *c, int *p_ldc)
/* 
 * Purpose
 * =======
 *
 * This routine computes the matrix product:
 *
 *      C   <-  alpha * op(A) * op(B)  +  beta * C .
 * 
 * where op(M) represents either M, M transpose, 
 * or M conjugate transpose.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage  of input matrices A, B, and C.
 *
 * transa  (input) enum blas_trans_type
 *         Operation to be done on matrix A before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 *
 * transb  (input) enum blas_trans_type
 *         Operation to be done on matrix B before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 * 
 * m n k   (input) int
 *         The dimensions of matrices A, B, and C.
 *         Matrix C is m-by-n matrix.
 *         Matrix A is m-by-k if A is not transposed, 
 *                     k-by-m otherwise.
 *         Matrix B is k-by-n if B is not transposed, 
 *                     n-by-k otherwise.
 *      
 * alpha   (input) float
 *
 * a       (input) const float*
 *         matrix A.
 * 
 * lda     (input) int
 *         leading dimension of A.
 * 
 * b       (input) const float*
 *         matrix B
 *
 * ldb     (input) int
 *         leading dimension of B.
 *
 * beta    (input) float
 *
 * c       (input/output) float*
 *         matrix C
 *
 * ldc     (input) int
 *         leading dimension of C.
 *
 */
{


    /* Integer Index Variables */
    int i, j, h;

    int ai, bj, ci;
    int aih, bhj, cij;		/* Index into matrices a, b, c during multiply */

    int incai, incaih;		/* Index increments for matrix a */
    int incbj, incbhj;		/* Index increments for matrix b */
    int incci, inccij;		/* Index increments for matrix c */

    /* Input Matrices */
    const float *a_i = a;
    const float *b_i = b;

    /* Output Matrix */
    float *c_i = c;

    /* Input Scalars */
    float alpha_i = *alpha;
    float beta_i = *beta;
    int m=*M, n=*N, k=*K;
    int lda=*p_lda, ldb=*p_ldb, ldc=*p_ldc;

    /* Temporary Floating-Point Variables */
    float a_elem;
    float b_elem;
    float c_elem;
    float prod;
    float sum;
    float tmp1;
    float tmp2;

    char order = COLMAJOR; /* For the time being, it is always COLMAJOR.
			      Eventually, it might change. */

#if MDEBUG
    const float *A = a, *B =b, *C = c, *C1 = c;
    printf("In Sgemm\n");
    printf("m=%d, n=%d, k=%d\n", m, n, k);
    printf("alpha=%f, beta=%f\n", alpha_i, beta_i);
    printf("\n");	  
    /* for(i=0; i<m*k; i++)  printf("%.1f ", *A++); printf("\n\n"); */
    /* for(i=0; i<n*k; i++)  printf("%.1f ", *B++); printf("\n\n"); */
    /* for(i=0; i<n*k; i++)  printf("%.1f ", *C++); printf("\n"); */
#endif
    
    /* Test for error conditions */
    if (m <= 0 || n <= 0 || k <= 0)
    {
	return;
    }

    if (order == COLMAJOR)
    {

	if (ldc < m)
	    return;

	if (*transa == 'n' || *transa == 'N')
	{
	    if (lda < m)
		return;
	}
	else
	{
	    if (lda < k)
		return;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    if (ldb < k)
		return;
	}
	else
	{
	    if (ldb < n)
		return;
	}

    }
    else
    {
	/* row major */
	if (ldc < n)
	    return;

	if (*transa == 'n' || *transa == 'N')
	{
	    if (lda < k)
		return;
	}
	else
	{
	    if (lda < m)
		return;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    if (ldb < n)
		return;
	}
	else
	{
	    if (ldb < k)
		return;
	}
    }

    /* Test for no-op */
    if (alpha_i == 0.0 && beta_i == 1.0)
    {
	return;
    }

    /* Set Index Parameters */
    if (order == COLMAJOR)
    {
	incci = 1;
	inccij = ldc;

	if (*transa == 'n' || *transa == 'N')
	{
	    incai = 1;
	    incaih = lda;
	}
	else
	{
	    incai = lda;
	    incaih = 1;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    incbj = ldb;
	    incbhj = 1;
	}
	else
	{
	    incbj = 1;
	    incbhj = ldb;
	}

    }
    else
    {
	/* row major */
	incci = ldc;
	inccij = 1;

	if (*transa == 'n' || *transa == 'N')
	{
	    incai = lda;
	    incaih = 1;
	}
	else
	{
	    incai = 1;
	    incaih = lda;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    incbj = 1;
	    incbhj = ldb;
	}
	else
	{
	    incbj = ldb;
	    incbhj = 1;
	}

    }



    /* Ajustment to increments */







    /* alpha = 0.  In this case, just return beta * C */
    if (alpha_i == 0.0)
    {

	ci = 0;
	for (i = 0; i < m; i++, ci += incci)
	{
	    cij = ci;
	    for (j = 0; j < n; j++, cij += inccij)
	    {
		c_elem = c_i[cij];
		tmp1 = c_elem * beta_i;
		c_i[cij] = tmp1;
	    }
	}

    }
    else if (alpha_i == 1.0)
    {

	/* Case alpha == 1. */

	if (beta_i == 0.0)
	{
	    /* Case alpha == 1, beta == 0.   We compute  C <--- A * B */

	    ci = 0;
	    ai = 0;
	    for (i = 0; i < m; i++, ci += incci, ai += incai)
	    {

		cij = ci;
		bj = 0;

		for (j = 0; j < n; j++, cij += inccij, bj += incbj)
		{

		    aih = ai;
		    bhj = bj;

		    sum = 0.0;

		    for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		    {
			a_elem = a_i[aih];
			b_elem = b_i[bhj];
			prod = a_elem * b_elem;
			sum = sum + prod;
		    }
		    tmp1 = sum;
		    c_i[cij] = tmp1;
		}
	    }

	}
	else
	{
	    /* Case alpha == 1, but beta != 0.
	       We compute   C <--- A * B + beta * C   */

	    ci = 0;
	    ai = 0;
	    for (i = 0; i < m; i++, ci += incci, ai += incai)
	    {

		cij = ci;
		bj = 0;

		for (j = 0; j < n; j++, cij += inccij, bj += incbj)
		{

		    aih = ai;
		    bhj = bj;

		    sum = 0.0;

		    for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		    {
			a_elem = a_i[aih];
			b_elem = b_i[bhj];
			prod = a_elem * b_elem;
			sum = sum + prod;
		    }

		    c_elem = c_i[cij];
		    tmp2 = c_elem * beta_i;
		    tmp1 = sum;
		    tmp1 = tmp2 + tmp1;
		    c_i[cij] = tmp1;
		}
	    }
	}

    }
    else
    {

	/* The most general form,   C <-- alpha * A * B + beta * C  */
	ci = 0;
	ai = 0;
	for (i = 0; i < m; i++, ci += incci, ai += incai)
	{

	    cij = ci;
	    bj = 0;

	    for (j = 0; j < n; j++, cij += inccij, bj += incbj)
	    {

		aih = ai;
		bhj = bj;

		sum = 0.0;

		for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		{
		    a_elem = a_i[aih];
		    b_elem = b_i[bhj];
		    prod = a_elem * b_elem;
		    sum = sum + prod;
		}

		tmp1 = sum * alpha_i;
		c_elem = c_i[cij];
		tmp2 = c_elem * beta_i;
		tmp1 = tmp1 + tmp2;
		c_i[cij] = tmp1;
	    }
	}

    }

#if MDEBUG
    printf("m=%d, n=%d, k=%d\n", m, n, k);
    printf("lda=%d, ldb=%d, ldc=%d\n", lda, ldb, ldc);
    printf("\n");	  
    for(i=0; i<n*k; i++)  printf("%.1f ", *C1++); printf("\n");
#endif

}


void
xb_dgemm (char *transa, char *transb, int *M, int *N, int *K,
	  double *alpha, const double *a, int *p_lda, const double *b, 
	  int *p_ldb, double *beta, double *c, int *p_ldc)
/* 
 * Purpose
 * =======
 *
 * This routine computes the matrix product:
 *
 *      C   <-  alpha * op(A) * op(B)  +  beta * C .
 * 
 * where op(M) represents either M, M transpose, 
 * or M conjugate transpose.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage  of input matrices A, B, and C.
 *
 * transa  (input) enum blas_trans_type
 *         Operation to be done on matrix A before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 *
 * transb  (input) enum blas_trans_type
 *         Operation to be done on matrix B before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 * 
 * m n k   (input) int
 *         The dimensions of matrices A, B, and C.
 *         Matrix C is m-by-n matrix.
 *         Matrix A is m-by-k if A is not transposed, 
 *                     k-by-m otherwise.
 *         Matrix B is k-by-n if B is not transposed, 
 *                     n-by-k otherwise.
 *      
 * alpha   (input) double
 *
 * a       (input) const double*
 *         matrix A.
 * 
 * lda     (input) int
 *         leading dimension of A.
 * 
 * b       (input) const double*
 *         matrix B
 *
 * ldb     (input) int
 *         leading dimension of B.
 *
 * beta    (input) double
 *
 * c       (input/output) double*
 *         matrix C
 *
 * ldc     (input) int
 *         leading dimension of C.
 *
 */
{


    /* Integer Index Variables */
    int i, j, h;

    int ai, bj, ci;
    int aih, bhj, cij;		/* Index into matrices a, b, c during multiply */

    int incai, incaih;		/* Index increments for matrix a */
    int incbj, incbhj;		/* Index increments for matrix b */
    int incci, inccij;		/* Index increments for matrix c */

    /* Input Matrices */
    const double *a_i = a;
    const double *b_i = b;

    /* Output Matrix */
    double *c_i = c;

    /* Input Scalars */
    double alpha_i = *alpha;
    double beta_i = *beta;
    int m=*M, n=*N, k=*K;
    int lda=*p_lda, ldb=*p_ldb, ldc=*p_ldc;

    /* Temporary Floating-Point Variables */
    double a_elem;
    double b_elem;
    double c_elem;
    double prod;
    double sum;
    double tmp1;
    double tmp2;

    char order = COLMAJOR; /* For the time being, it is always COLMAJOR.
			      Eventually, it might change. */

#if MDEBUG
    const double *A = a, *B =b, *C = c;
    printf("m=%d, n=%d, k=%d\n", m, n, k);
    printf("lda=%d, ldb=%d, ldc=%d\n", lda, ldb, ldc);
    printf("alpha=%lf, beta=%lf\n", alpha_i, beta_i);
    printf("\n");	  
    for(i=0; i<m*k; i++)  printf("%lf ", *A++); printf("\n\n");
    for(i=0; i<n*k; i++)  printf("%lf ", *B++); printf("\n\n");
    for(i=0; i<n*k; i++)  printf("%lf ", *C++); printf("\n");
#endif


    /* Test for error conditions */
    if (m <= 0 || n <= 0 || k <= 0)
    {
	return;
    }

    if (order == COLMAJOR)
    {
	if (ldc < m)
	    return;
	
	if (*transa == 'n' || *transa == 'N')
	{
	    if (lda < m)
		return;
	}
	else
	{
	    if (lda < k)
		return;
	}
	
	if (*transb == 'n' || *transb == 'N')
	{
	    if (ldb < k)
		return;
	}
	else
	{
	    if (ldb < n)
		return;
	}

    }
    else
    {
	/* row major */
	if (ldc < n)
	    return;

	if (*transa == 'n' || *transa == 'N')
	{
	    if (lda < k)
		return;
	}
	else
	{
	    if (lda < m)
		return;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    if (ldb < n)
		return;
	}
	else
	{
	    if (ldb < k)
		return;
	}
    }

    /* Test for no-op */
    if (alpha_i == 0.0 && beta_i == 1.0)
    {
	return;
    }

    /* Set Index Parameters */
    if (order == COLMAJOR)
    {

	incci = 1;
	inccij = ldc;

	if (*transa == 'n' || *transa == 'N')
	{
	    incai = 1;
	    incaih = lda;
	}
	else
	{
	    incai = lda;
	    incaih = 1;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    incbj = ldb;
	    incbhj = 1;
	}
	else
	{
	    incbj = 1;
	    incbhj = ldb;
	}

    }
    else
    {
	/* row major */
	incci = ldc;
	inccij = 1;

	if (*transa == 'n' || *transa == 'N')
	{
	    incai = lda;
	    incaih = 1;
	}
	else
	{
	    incai = 1;
	    incaih = lda;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    incbj = 1;
	    incbhj = ldb;
	}
	else
	{
	    incbj = ldb;
	    incbhj = 1;
	}

    }


    /* Ajustment to increments */


    /* alpha = 0.  In this case, just return beta * C */
    if (alpha_i == 0.0)
    {
#     if MDEBUG
      printf("Case 0 (alpha=0.0, beta==x): C <--- beta * C\n");
#     endif
	ci = 0;
	for (i = 0; i < m; i++, ci += incci)
	{
	    cij = ci;
	    for (j = 0; j < n; j++, cij += inccij)
	    {
		c_elem = c_i[cij];
		tmp1 = c_elem * beta_i;
		c_i[cij] = tmp1;
	    }
	}

    }
    else if (alpha_i == 1.0)
    {

	/* Case alpha == 1. */

	if (beta_i == 0.0)
	{
	    /* Case alpha == 1, beta == 0.   We compute  C <--- A * B */
#         if MDEBUG
	  printf("Case 1 (alpha=1.0, beta==0): C <--- A * B\n");
#         endif	  
	    ci = 0;
	    ai = 0;
	    for (i = 0; i < m; i++, ci += incci, ai += incai)
	    {

		cij = ci;
		bj = 0;

		for (j = 0; j < n; j++, cij += inccij, bj += incbj)
		{

		    aih = ai;
		    bhj = bj;

		    sum = 0.0;

		    for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		    {
			a_elem = a_i[aih];
			b_elem = b_i[bhj];
			prod = a_elem * b_elem;
			sum = sum + prod;
		    }
		    tmp1 = sum;
		    c_i[cij] = tmp1;
		}
	    }

	}
	else
	{
	    /* Case alpha == 1, but beta != 0.
	       We compute   C <--- A * B + beta * C   */
#         if MDEBUG	  
	  printf("Case 2 (alpha=1.0, beta!=0): A * B + beta * C\n");
#         endif
	    ci = 0;
	    ai = 0;
	    for (i = 0; i < m; i++, ci += incci, ai += incai)
	    {

		cij = ci;
		bj = 0;

		for (j = 0; j < n; j++, cij += inccij, bj += incbj)
		{

		    aih = ai;
		    bhj = bj;

		    sum = 0.0;

		    for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		    {
			a_elem = a_i[aih];
			b_elem = b_i[bhj];
			prod = a_elem * b_elem;
			sum = sum + prod;
		    }

		    c_elem = c_i[cij];
		    tmp2 = c_elem * beta_i;
		    tmp1 = sum;
		    tmp1 = tmp2 + tmp1;
		    c_i[cij] = tmp1;
		}
	    }
	}

    }
    else
    {

	/* The most general form,   C <-- alpha * A * B + beta * C  */
#     if MDEBUG      
      printf("Case 3 (alpha=x, beta=x): alpha * A * B + beta * C\n");
#     endif
	ci = 0;
	ai = 0;
	for (i = 0; i < m; i++, ci += incci, ai += incai)
	{

	    cij = ci;
	    bj = 0;

	    for (j = 0; j < n; j++, cij += inccij, bj += incbj)
	    {

		aih = ai;
		bhj = bj;

		sum = 0.0;

		for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		{
		    a_elem = a_i[aih];
		    b_elem = b_i[bhj];
		    prod = a_elem * b_elem;
		    sum = sum + prod;
		}

		tmp1 = sum * alpha_i;
		c_elem = c_i[cij];
		tmp2 = c_elem * beta_i;
		tmp1 = tmp1 + tmp2;
		c_i[cij] = tmp1;
	    }
	}

    }



}


void
xb_zgemm (char * transa, char *transb, int *M, int *N, int *K,
	  const void *alpha, const void *a, int *p_lda, const void *b, 
	  int *p_ldb, const void *beta, void *c, int *p_ldc)
/* 
 * Purpose
 * =======
 *
 * This routine computes the matrix product:
 *
 *      C   <-  alpha * op(A) * op(B)  +  beta * C .
 * 
 * where op(M) represents either M, M transpose, 
 * or M conjugate transpose.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage  of input matrices A, B, and C.
 *
 * transa  (input) enum blas_trans_type
 *         Operation to be done on matrix A before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 *
 * transb  (input) enum blas_trans_type
 *         Operation to be done on matrix B before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 * 
 * m n k   (input) int
 *         The dimensions of matrices A, B, and C.
 *         Matrix C is m-by-n matrix.
 *         Matrix A is m-by-k if A is not transposed, 
 *                     k-by-m otherwise.
 *         Matrix B is k-by-n if B is not transposed, 
 *                     n-by-k otherwise.
 *      
 * alpha   (input) const void*
 *
 * a       (input) const void*
 *         matrix A.
 * 
 * lda     (input) int
 *         leading dimension of A.
 * 
 * b       (input) const void*
 *         matrix B
 *
 * ldb     (input) int
 *         leading dimension of B.
 *
 * beta    (input) const void*
 *
 * c       (input/output) void*
 *         matrix C
 *
 * ldc     (input) int
 *         leading dimension of C.
 *
 */
{


    /* Integer Index Variables */
    int i, j, h;

    int ai, bj, ci;
    int aih, bhj, cij;		/* Index into matrices a, b, c during multiply */

    int incai, incaih;		/* Index increments for matrix a */
    int incbj, incbhj;		/* Index increments for matrix b */
    int incci, inccij;		/* Index increments for matrix c */

    /* Input Matrices */
    const double *a_i = (double *) a;
    const double *b_i = (double *) b;

    /* Output Matrix */
    double *c_i = (double *) c;

    /* Input Scalars */
    double *alpha_i = (double *) alpha;
    double *beta_i = (double *) beta;
    int m=*M, n=*N, k=*K;
    int lda=*p_lda, ldb=*p_ldb, ldc=*p_ldc;

    /* Temporary Floating-Point Variables */
    double a_elem[2];
    double b_elem[2];
    double c_elem[2];
    double prod[2];
    double sum[2];
    double tmp1[2];
    double tmp2[2];

    char order = COLMAJOR; /* For the time being, it is always COLMAJOR.
			      Eventually, it might change. */

#if MDEBUG
    printf("In Zgemm\n");
#endif

    /* Test for error conditions */
    if (m <= 0 || n <= 0 || k <= 0)
    {
	return;
    }

    if (order == COLMAJOR)
    {

	if (ldc < m)
	    return;

	if (*transa == 'n' || *transa == 'N')
	{
	    if (lda < m)
		return;
	}
	else
	{
	    if (lda < k)
		return;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    if (ldb < k)
		return;
	}
	else
	{
	    if (ldb < n)
		return;
	}

    }
    else
    {
	/* row major */
	if (ldc < n)
	    return;

	if (*transa == 'n' || *transa == 'N')
	{
	    if (lda < k)
		return;
	}
	else
	{
	    if (lda < m)
		return;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    if (ldb < n)
		return;
	}
	else
	{
	    if (ldb < k)
		return;
	}
    }

    /* Test for no-op */
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	&& (beta_i[0] == 1.0 && beta_i[1] == 0.0))
    {
	return;
    }

    /* Set Index Parameters */
    if (order == COLMAJOR)
    {
	incci = 1;
	inccij = ldc;

	if (*transa == 'n' || *transa == 'N')
	{
	    incai = 1;
	    incaih = lda;
	}
	else
	{
	    incai = lda;
	    incaih = 1;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    incbj = ldb;
	    incbhj = 1;
	}
	else
	{
	    incbj = 1;
	    incbhj = ldb;
	}

    }
    else
    {
	/* row major */
	incci = ldc;
	inccij = 1;

	if (*transa == 'n' || *transa == 'N')
	{
	    incai = lda;
	    incaih = 1;
	}
	else
	{
	    incai = 1;
	    incaih = lda;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    incbj = 1;
	    incbhj = ldb;
	}
	else
	{
	    incbj = ldb;
	    incbhj = 1;
	}

    }



    /* Ajustment to increments */
    incci *= 2;
    inccij *= 2;
    incai *= 2;
    incaih *= 2;
    incbj *= 2;
    incbhj *= 2;

    /* alpha = 0.  In this case, just return beta * C */
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0)
    {

	ci = 0;
	for (i = 0; i < m; i++, ci += incci)
	{
	    cij = ci;
	    for (j = 0; j < n; j++, cij += inccij)
	    {
		c_elem[0] = c_i[cij];
		c_elem[1] = c_i[cij + 1];
		{
		    tmp1[0] =
			(double) c_elem[0] * beta_i[0] -
			(double) c_elem[1] * beta_i[1];
		    tmp1[1] =
			(double) c_elem[0] * beta_i[1] +
			(double) c_elem[1] * beta_i[0];
		}
		c_i[cij] = tmp1[0];
		c_i[cij + 1] = tmp1[1];
	    }
	}

    }
    else if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0))
    {

	/* Case alpha == 1. */

	if (beta_i[0] == 0.0 && beta_i[1] == 0.0)
	{
	    /* Case alpha == 1, beta == 0.   We compute  C <--- A * B */

	    ci = 0;
	    ai = 0;
	    for (i = 0; i < m; i++, ci += incci, ai += incai)
	    {

		cij = ci;
		bj = 0;

		for (j = 0; j < n; j++, cij += inccij, bj += incbj)
		{

		    aih = ai;
		    bhj = bj;

		    sum[0] = sum[1] = 0.0;

		    for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		    {
			a_elem[0] = a_i[aih];
			a_elem[1] = a_i[aih + 1];
			b_elem[0] = b_i[bhj];
			b_elem[1] = b_i[bhj + 1];
			if (*transa == 'c' || *transa == 'C')
			{
			  a_elem[1] = -a_elem[1];
			}
			if (*transb == 'c' || *transb == 'C')
			{
			    b_elem[1] = -b_elem[1];
			}
			
			{
			    prod[0] =
				(double) a_elem[0] * b_elem[0] -
				(double) a_elem[1] * b_elem[1];
			    prod[1] =
				(double) a_elem[0] * b_elem[1] +
				(double) a_elem[1] * b_elem[0];
			}
			sum[0] = sum[0] + prod[0];
			sum[1] = sum[1] + prod[1];
		    }
		    tmp1[0] = sum[0];
		    tmp1[1] = sum[1];
		    c_i[cij] = tmp1[0];
		    c_i[cij + 1] = tmp1[1];
		}
	    }

	}
	else
	{
	    /* Case alpha == 1, but beta != 0.
	       We compute   C <--- A * B + beta * C   */

	    ci = 0;
	    ai = 0;
	    for (i = 0; i < m; i++, ci += incci, ai += incai)
	    {

		cij = ci;
		bj = 0;

		for (j = 0; j < n; j++, cij += inccij, bj += incbj)
		{

		    aih = ai;
		    bhj = bj;

		    sum[0] = sum[1] = 0.0;

		    for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		    {
			a_elem[0] = a_i[aih];
			a_elem[1] = a_i[aih + 1];
			b_elem[0] = b_i[bhj];
			b_elem[1] = b_i[bhj + 1];
			if (*transa == 'c' || *transa == 'C')
			{
			    a_elem[1] = -a_elem[1];
			}
			if (*transb == 'c' || *transb == 'C')
			{
			    b_elem[1] = -b_elem[1];
			}
			{
			    prod[0] =
				(double) a_elem[0] * b_elem[0] -
				(double) a_elem[1] * b_elem[1];
			    prod[1] =
				(double) a_elem[0] * b_elem[1] +
				(double) a_elem[1] * b_elem[0];
			}
			sum[0] = sum[0] + prod[0];
			sum[1] = sum[1] + prod[1];
		    }

		    c_elem[0] = c_i[cij];
		    c_elem[1] = c_i[cij + 1];
		    {
			tmp2[0] =
			    (double) c_elem[0] * beta_i[0] -
			    (double) c_elem[1] * beta_i[1];
			tmp2[1] =
			    (double) c_elem[0] * beta_i[1] +
			    (double) c_elem[1] * beta_i[0];
		    }
		    tmp1[0] = sum[0];
		    tmp1[1] = sum[1];
		    tmp1[0] = tmp2[0] + tmp1[0];
		    tmp1[1] = tmp2[1] + tmp1[1];
		    c_i[cij] = tmp1[0];
		    c_i[cij + 1] = tmp1[1];
		}
	    }
	}

    }
    else
    {

	/* The most general form,   C <-- alpha * A * B + beta * C  */
	ci = 0;
	ai = 0;
	for (i = 0; i < m; i++, ci += incci, ai += incai)
	{

	    cij = ci;
	    bj = 0;

	    for (j = 0; j < n; j++, cij += inccij, bj += incbj)
	    {

		aih = ai;
		bhj = bj;

		sum[0] = sum[1] = 0.0;

		for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		{
		    a_elem[0] = a_i[aih];
		    a_elem[1] = a_i[aih + 1];
		    b_elem[0] = b_i[bhj];
		    b_elem[1] = b_i[bhj + 1];
		    if (*transa == 'c' || *transa == 'C')
		    {
			a_elem[1] = -a_elem[1];
		    }
		    if (*transb == 'c' || *transb == 'C')
		    {
			b_elem[1] = -b_elem[1];
		    }
		    {
			prod[0] =
			    (double) a_elem[0] * b_elem[0] -
			    (double) a_elem[1] * b_elem[1];
			prod[1] =
			    (double) a_elem[0] * b_elem[1] +
			    (double) a_elem[1] * b_elem[0];
		    }
		    sum[0] = sum[0] + prod[0];
		    sum[1] = sum[1] + prod[1];
		}

		{
		    tmp1[0] =
			(double) sum[0] * alpha_i[0] -
			(double) sum[1] * alpha_i[1];
		    tmp1[1] =
			(double) sum[0] * alpha_i[1] +
			(double) sum[1] * alpha_i[0];
		}
		c_elem[0] = c_i[cij];
		c_elem[1] = c_i[cij + 1];
		{
		    tmp2[0] =
			(double) c_elem[0] * beta_i[0] -
			(double) c_elem[1] * beta_i[1];
		    tmp2[1] =
			(double) c_elem[0] * beta_i[1] +
			(double) c_elem[1] * beta_i[0];
		}
		tmp1[0] = tmp1[0] + tmp2[0];
		tmp1[1] = tmp1[1] + tmp2[1];
		c_i[cij] = tmp1[0];
		c_i[cij + 1] = tmp1[1];
	    }
	}

    }



}


void
xb_cgemm (char * transa, char *transb, int *M, int *N, int *K,
	  const void *alpha, const void *a, int *p_lda, const void *b, 
	  int *p_ldb, const void *beta, void *c, int *p_ldc)
/* 
 * Purpose
 * =======
 *
 * This routine computes the matrix product:
 *
 *      C   <-  alpha * op(A) * op(B)  +  beta * C .
 * 
 * where op(M) represents either M, M transpose, 
 * or M conjugate transpose.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage  of input matrices A, B, and C.
 *
 * transa  (input) enum blas_trans_type
 *         Operation to be done on matrix A before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 *
 * transb  (input) enum blas_trans_type
 *         Operation to be done on matrix B before multiplication.
 *         Can be no operation, transposition, or conjugate transposition.
 * 
 * m n k   (input) int
 *         The dimensions of matrices A, B, and C.
 *         Matrix C is m-by-n matrix.
 *         Matrix A is m-by-k if A is not transposed, 
 *                     k-by-m otherwise.
 *         Matrix B is k-by-n if B is not transposed, 
 *                     n-by-k otherwise.
 *      
 * alpha   (input) const void*
 *
 * a       (input) const void*
 *         matrix A.
 * 
 * lda     (input) int
 *         leading dimension of A.
 * 
 * b       (input) const void*
 *         matrix B
 *
 * ldb     (input) int
 *         leading dimension of B.
 *
 * beta    (input) const void*
 *
 * c       (input/output) void*
 *         matrix C
 *
 * ldc     (input) int
 *         leading dimension of C.
 *
 */
{


    /* Integer Index Variables */
    int i, j, h;

    int ai, bj, ci;
    int aih, bhj, cij;		/* Index into matrices a, b, c during multiply */

    int incai, incaih;		/* Index increments for matrix a */
    int incbj, incbhj;		/* Index increments for matrix b */
    int incci, inccij;		/* Index increments for matrix c */

    /* Input Matrices */
    const float *a_i = (float *) a;
    const float *b_i = (float *) b;

    /* Output Matrix */
    float *c_i = (float *) c;

    /* Input Scalars */
    float *alpha_i = (float *) alpha;
    float *beta_i = (float *) beta;
    int m=*M, n=*N, k=*K;
    int lda=*p_lda, ldb=*p_ldb, ldc=*p_ldc;

    /* Temporary Floating-Point Variables */
    float a_elem[2];
    float b_elem[2];
    float c_elem[2];
    float prod[2];
    float sum[2];
    float tmp1[2];
    float tmp2[2];

    char order = COLMAJOR; /* For the time being, it is always COLMAJOR.
			      Eventually, it might change. */

#if MDEBUG
    printf("In Zgemm\n");
#endif

    /* Test for error conditions */
    if (m <= 0 || n <= 0 || k <= 0)
    {
	return;
    }

    if (order == COLMAJOR)
    {

	if (ldc < m)
	    return;

	if (*transa == 'n' || *transa == 'N')
	{
	    if (lda < m)
		return;
	}
	else
	{
	    if (lda < k)
		return;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    if (ldb < k)
		return;
	}
	else
	{
	    if (ldb < n)
		return;
	}

    }
    else
    {
	/* row major */
	if (ldc < n)
	    return;

	if (*transa == 'n' || *transa == 'N')
	{
	    if (lda < k)
		return;
	}
	else
	{
	    if (lda < m)
		return;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    if (ldb < n)
		return;
	}
	else
	{
	    if (ldb < k)
		return;
	}
    }

    /* Test for no-op */
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	&& (beta_i[0] == 1.0 && beta_i[1] == 0.0))
    {
	return;
    }

    /* Set Index Parameters */
    if (order == COLMAJOR)
    {
	incci = 1;
	inccij = ldc;

	if (*transa == 'n' || *transa == 'N')
	{
	    incai = 1;
	    incaih = lda;
	}
	else
	{
	    incai = lda;
	    incaih = 1;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    incbj = ldb;
	    incbhj = 1;
	}
	else
	{
	    incbj = 1;
	    incbhj = ldb;
	}

    }
    else
    {
	/* row major */
	incci = ldc;
	inccij = 1;

	if (*transa == 'n' || *transa == 'N')
	{
	    incai = lda;
	    incaih = 1;
	}
	else
	{
	    incai = 1;
	    incaih = lda;
	}

	if (*transb == 'n' || *transb == 'N')
	{
	    incbj = 1;
	    incbhj = ldb;
	}
	else
	{
	    incbj = ldb;
	    incbhj = 1;
	}

    }



    /* Ajustment to increments */
    incci *= 2;
    inccij *= 2;
    incai *= 2;
    incaih *= 2;
    incbj *= 2;
    incbhj *= 2;

    /* alpha = 0.  In this case, just return beta * C */
    if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0)
    {

	ci = 0;
	for (i = 0; i < m; i++, ci += incci)
	{
	    cij = ci;
	    for (j = 0; j < n; j++, cij += inccij)
	    {
		c_elem[0] = c_i[cij];
		c_elem[1] = c_i[cij + 1];
		{
		    tmp1[0] =
			(float) c_elem[0] * beta_i[0] -
			(float) c_elem[1] * beta_i[1];
		    tmp1[1] =
			(float) c_elem[0] * beta_i[1] +
			(float) c_elem[1] * beta_i[0];
		}
		c_i[cij] = tmp1[0];
		c_i[cij + 1] = tmp1[1];
	    }
	}

    }
    else if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0))
    {

	/* Case alpha == 1. */

	if (beta_i[0] == 0.0 && beta_i[1] == 0.0)
	{
	    /* Case alpha == 1, beta == 0.   We compute  C <--- A * B */

	    ci = 0;
	    ai = 0;
	    for (i = 0; i < m; i++, ci += incci, ai += incai)
	    {

		cij = ci;
		bj = 0;

		for (j = 0; j < n; j++, cij += inccij, bj += incbj)
		{

		    aih = ai;
		    bhj = bj;

		    sum[0] = sum[1] = 0.0;

		    for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		    {
			a_elem[0] = a_i[aih];
			a_elem[1] = a_i[aih + 1];
			b_elem[0] = b_i[bhj];
			b_elem[1] = b_i[bhj + 1];
			if (*transa == 'c' || *transa == 'C')
			{
			  a_elem[1] = -a_elem[1];
			}
			if (*transb == 'c' || *transb == 'C')
			{
			    b_elem[1] = -b_elem[1];
			}
			
			{
			    prod[0] =
				(float) a_elem[0] * b_elem[0] -
				(float) a_elem[1] * b_elem[1];
			    prod[1] =
				(float) a_elem[0] * b_elem[1] +
				(float) a_elem[1] * b_elem[0];
			}
			sum[0] = sum[0] + prod[0];
			sum[1] = sum[1] + prod[1];
		    }
		    tmp1[0] = sum[0];
		    tmp1[1] = sum[1];
		    c_i[cij] = tmp1[0];
		    c_i[cij + 1] = tmp1[1];
		}
	    }

	}
	else
	{
	    /* Case alpha == 1, but beta != 0.
	       We compute   C <--- A * B + beta * C   */

	    ci = 0;
	    ai = 0;
	    for (i = 0; i < m; i++, ci += incci, ai += incai)
	    {

		cij = ci;
		bj = 0;

		for (j = 0; j < n; j++, cij += inccij, bj += incbj)
		{

		    aih = ai;
		    bhj = bj;

		    sum[0] = sum[1] = 0.0;

		    for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		    {
			a_elem[0] = a_i[aih];
			a_elem[1] = a_i[aih + 1];
			b_elem[0] = b_i[bhj];
			b_elem[1] = b_i[bhj + 1];
			if (*transa == 'c' || *transa == 'C')
			{
			    a_elem[1] = -a_elem[1];
			}
			if (*transb == 'c' || *transb == 'C')
			{
			    b_elem[1] = -b_elem[1];
			}
			{
			    prod[0] =
				(float) a_elem[0] * b_elem[0] -
				(float) a_elem[1] * b_elem[1];
			    prod[1] =
				(float) a_elem[0] * b_elem[1] +
				(float) a_elem[1] * b_elem[0];
			}
			sum[0] = sum[0] + prod[0];
			sum[1] = sum[1] + prod[1];
		    }

		    c_elem[0] = c_i[cij];
		    c_elem[1] = c_i[cij + 1];
		    {
			tmp2[0] =
			    (float) c_elem[0] * beta_i[0] -
			    (float) c_elem[1] * beta_i[1];
			tmp2[1] =
			    (float) c_elem[0] * beta_i[1] +
			    (float) c_elem[1] * beta_i[0];
		    }
		    tmp1[0] = sum[0];
		    tmp1[1] = sum[1];
		    tmp1[0] = tmp2[0] + tmp1[0];
		    tmp1[1] = tmp2[1] + tmp1[1];
		    c_i[cij] = tmp1[0];
		    c_i[cij + 1] = tmp1[1];
		}
	    }
	}

    }
    else
    {

	/* The most general form,   C <-- alpha * A * B + beta * C  */
	ci = 0;
	ai = 0;
	for (i = 0; i < m; i++, ci += incci, ai += incai)
	{

	    cij = ci;
	    bj = 0;

	    for (j = 0; j < n; j++, cij += inccij, bj += incbj)
	    {

		aih = ai;
		bhj = bj;

		sum[0] = sum[1] = 0.0;

		for (h = 0; h < k; h++, aih += incaih, bhj += incbhj)
		{
		    a_elem[0] = a_i[aih];
		    a_elem[1] = a_i[aih + 1];
		    b_elem[0] = b_i[bhj];
		    b_elem[1] = b_i[bhj + 1];
		    if (*transa == 'c' || *transa == 'C')
		    {
			a_elem[1] = -a_elem[1];
		    }
		    if (*transb == 'c' || *transb == 'C')
		    {
			b_elem[1] = -b_elem[1];
		    }
		    {
			prod[0] =
			    (float) a_elem[0] * b_elem[0] -
			    (float) a_elem[1] * b_elem[1];
			prod[1] =
			    (float) a_elem[0] * b_elem[1] +
			    (float) a_elem[1] * b_elem[0];
		    }
		    sum[0] = sum[0] + prod[0];
		    sum[1] = sum[1] + prod[1];
		}

		{
		    tmp1[0] =
			(float) sum[0] * alpha_i[0] -
			(float) sum[1] * alpha_i[1];
		    tmp1[1] =
			(float) sum[0] * alpha_i[1] +
			(float) sum[1] * alpha_i[0];
		}
		c_elem[0] = c_i[cij];
		c_elem[1] = c_i[cij + 1];
		{
		    tmp2[0] =
			(float) c_elem[0] * beta_i[0] -
			(float) c_elem[1] * beta_i[1];
		    tmp2[1] =
			(float) c_elem[0] * beta_i[1] +
			(float) c_elem[1] * beta_i[0];
		}
		tmp1[0] = tmp1[0] + tmp2[0];
		tmp1[1] = tmp1[1] + tmp2[1];
		c_i[cij] = tmp1[0];
		c_i[cij + 1] = tmp1[1];
	    }
	}

    }



}
