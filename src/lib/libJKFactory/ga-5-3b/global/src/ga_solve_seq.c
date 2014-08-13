/**
 * GA_Lu_solve_seq.c: Implemented with CLINPACK routines. Uses LINPACK
 * routines if ENABLE_F77 is not defined,
 * else uses internal or external lapack.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_MATH_H
#   include <math.h>
#endif

#include "globalp.h"
#include "macdecls.h"
#include "ga-papi.h"
#include "ga-wapi.h"
#include "galinalg.h"

#define REAL double
#define ZERO 0.0e0
#define ONE 1.0e0

/*----------------------*/ 
void LP_daxpy(n,da,dx,incx,dy,incy)
/*
     constant times a vector plus a vector.
     jack dongarra, linpack, 3/11/78.
*/
REAL dx[],dy[],da;
int incx,incy,n;
{
	int i,ix,iy;

	if(n <= 0) return;
	if (da == ZERO) return;

	if(incx != 1 || incy != 1) {

		/* code for unequal increments or equal increments
		   not equal to 1 					*/

		ix = 0;
		iy = 0;
		if(incx < 0) ix = (-n+1)*incx;
		if(incy < 0)iy = (-n+1)*incy;
		for (i = 0;i < n; i++) {
			dy[iy] = dy[iy] + da*dx[ix];
			ix = ix + incx;
			iy = iy + incy;
		}
      		return;
	}

	/* code for both increments equal to 1 */

#ifdef ROLL
	for (i = 0;i < n; i++) {
		dy[i] = dy[i] + da*dx[i];
	}
#endif
#ifdef UNROLL

	m = n % 4;
	if ( m != 0) {
		for (i = 0; i < m; i++) 
			dy[i] = dy[i] + da*dx[i];
		if (n < 4) return;
	}
	for (i = m; i < n; i = i + 4) {
		dy[i] = dy[i] + da*dx[i];
		dy[i+1] = dy[i+1] + da*dx[i+1];
		dy[i+2] = dy[i+2] + da*dx[i+2];
		dy[i+3] = dy[i+3] + da*dx[i+3];
	}
#endif
}
   
/*----------------------*/ 

REAL LP_ddot(n,dx,incx,dy,incy)
/*
     forms the dot product of two vectors.
     jack dongarra, linpack, 3/11/78.
*/
REAL dx[],dy[];

int incx,incy,n;
{
	REAL dtemp;
	int i,ix,iy;

	dtemp = ZERO;

	if(n <= 0) return(ZERO);

	if(incx != 1 || incy != 1) {

		/* code for unequal increments or equal increments
		   not equal to 1					*/

		ix = 0;
		iy = 0;
		if (incx < 0) ix = (-n+1)*incx;
		if (incy < 0) iy = (-n+1)*incy;
		for (i = 0;i < n; i++) {
			dtemp = dtemp + dx[ix]*dy[iy];
			ix = ix + incx;
			iy = iy + incy;
		}
		return(dtemp);
	}

	/* code for both increments equal to 1 */

#ifdef ROLL
	for (i=0;i < n; i++)
		dtemp = dtemp + dx[i]*dy[i];
#endif
#ifdef UNROLL

	m = n % 5;
	if (m != 0) {
		for (i = 0; i < m; i++)
			dtemp = dtemp + dx[i]*dy[i];
		if (n < 5) return(dtemp);
	}
	for (i = m; i < n; i = i + 5) {
		dtemp = dtemp + dx[i]*dy[i] +
		dx[i+1]*dy[i+1] + dx[i+2]*dy[i+2] +
		dx[i+3]*dy[i+3] + dx[i+4]*dy[i+4];
	}
#endif
	return(dtemp);
}

/*----------------------*/ 
void LP_dscal(n,da,dx,incx)

/*     scales a vector by a constant.
      jack dongarra, linpack, 3/11/78.
*/
REAL da,dx[];
int n, incx;
{
	int i,nincx;

	if(n <= 0)return;
	if(incx != 1) {

		/* code for increment not equal to 1 */

		nincx = n*incx;
		for (i = 0; i < nincx; i = i + incx)
			dx[i] = da*dx[i];
		return;
	}

	/* code for increment equal to 1 */

#ifdef ROLL
	for (i = 0; i < n; i++)
		dx[i] = da*dx[i];
#endif
#ifdef UNROLL

	m = n % 5;
	if (m != 0) {
		for (i = 0; i < m; i++)
			dx[i] = da*dx[i];
		if (n < 5) return;
	}
	for (i = m; i < n; i = i + 5){
		dx[i] = da*dx[i];
		dx[i+1] = da*dx[i+1];
		dx[i+2] = da*dx[i+2];
		dx[i+3] = da*dx[i+3];
		dx[i+4] = da*dx[i+4];
	}
#endif

}

/*----------------------*/ 
int LP_idamax(n,dx,incx)

/*
     finds the index of element having max. absolute value.
     jack dongarra, linpack, 3/11/78.
*/

REAL dx[];
int incx,n;
{
	REAL dmax;
	int i, ix, itemp = 0;

	if( n < 1 ) return(-1);
	if(n ==1 ) return(0);
	if(incx != 1) {

		/* code for increment not equal to 1 */

		ix = 1;
		dmax = fabs((double)dx[0]);
		ix = ix + incx;
		for (i = 1; i < n; i++) {
			if(fabs((double)dx[ix]) > dmax)  {
				itemp = i;
				dmax = fabs((double)dx[ix]);
			}
			ix = ix + incx;
		}
	}
	else {

		/* code for increment equal to 1 */

		itemp = 0;
		dmax = fabs((double)dx[0]);
		for (i = 1; i < n; i++) {
			if(fabs((double)dx[i]) > dmax) {
				itemp = i;
				dmax = fabs((double)dx[i]);
			}
		}
	}
	return (itemp);
}



/*----------------------*/ 
void LP_dgefa(a,lda,n,ipvt,info)
REAL a[];
int lda,n,ipvt[],*info;

/* We would like to declare a[][lda], but c does not allow it.  In this
function, references to a[i][j] are written a[lda*i+j].  */
/*
     LP_dgefa factors a double precision matrix by gaussian elimination.

     LP_dgefa is usually called by dgeco, but it can be called
     directly with a saving in time if  rcond  is not needed.
     (time for dgeco) = (1 + 9/n)*(time for LP_dgefa) .

     on entry

        a       REAL precision[n][lda]
                the matrix to be factored.

        lda     integer
                the leading dimension of the array  a .

        n       integer
                the order of the matrix  a .

     on return

        a       an upper triangular matrix and the multipliers
                which were used to obtain it.
                the factorization can be written  a = l*u  where
                l  is a product of permutation and unit lower
                triangular matrices and  u  is upper triangular.

        ipvt    integer[n]
                an integer vector of pivot indices.

        info    integer
                = 0  normal value.
                = k  if  u[k][k] .eq. 0.0 .  this is not an error
                     condition for this subroutine, but it does
                     indicate that LP_dgesl or dgedi will divide by zero
                     if called.  use  rcond  in dgeco for a reliable
                     indication of singularity.

     linpack. this version dated 08/14/78 .
     cleve moler, university of new mexico, argonne national lab.

     functions

     blas LP_daxpy,LP_dscal,LP_idamax
*/

{
/*     internal variables	*/

  REAL t;
  int j,k,kp1,l,nm1;


/*     gaussian elimination with partial pivoting	*/
	*info = 0;
	nm1 = n - 1;
	if (nm1 >=  0) {
		for (k = 0; k < nm1; k++) {
			kp1 = k + 1;

          		/* find l = pivot index	*/

			l = LP_idamax(n-k,&a[lda*k+k],1) + k;
			ipvt[k] = l;

			/* zero pivot implies this column already 
			   triangularized */
			if (a[lda*k+l] != ZERO) {

				/* interchange if necessary */

				if (l != k) {
					t = a[lda*k+l];
					a[lda*k+l] = a[lda*k+k];
					a[lda*k+k] = t; 
				}

				/* compute multipliers */

				t = -ONE/a[lda*k+k];
				LP_dscal(n-(k+1),t,&a[lda*k+k+1],1);

				/* row elimination with column indexing */

				for (j = kp1; j < n; j++) {
					t = a[lda*j+l];
					if (l != k) {
						a[lda*j+l] = a[lda*j+k];
						a[lda*j+k] = t;
					}
					LP_daxpy(n-(k+1),t,&a[lda*k+k+1],1,
					      &a[lda*j+k+1],1);
  				} 
  			}
			else { 
            			*info = k;
			}
		} 
	}
	ipvt[n-1] = n-1;
	if (a[lda*(n-1)+(n-1)] == ZERO) *info = n-1;
}

/*----------------------*/ 

void LP_dgesl(a,lda,n,ipvt,b,job)
int lda,n,ipvt[],job;
REAL a[],b[];

/* We would like to declare a[][lda], but c does not allow it.  In this
function, references to a[i][j] are written a[lda*i+j].  */

/*
     LP_dgesl solves the double precision system
     a * x = b  or  trans(a) * x = b
     using the factors computed by dgeco or LP_dgefa.

     on entry

        a       double precision[n][lda]
                the output from dgeco or LP_dgefa.

        lda     integer
                the leading dimension of the array  a .

        n       integer
                the order of the matrix  a .

        ipvt    integer[n]
                the pivot vector from dgeco or LP_dgefa.

        b       double precision[n]
                the right hand side vector.

        job     integer
                = 0         to solve  a*x = b ,
                = nonzero   to solve  trans(a)*x = b  where
                            trans(a)  is the transpose.

    on return

        b       the solution vector  x .

     error condition

        a division by zero will occur if the input factor contains a
        zero on the diagonal.  technically this indicates singularity
        but it is often caused by improper arguments or improper
        setting of lda .  it will not occur if the subroutines are
        called correctly and if dgeco has set rcond .gt. 0.0
        or LP_dgefa has set info .eq. 0 .

     to compute  inverse(a) * c  where  c  is a matrix
     with  p  columns
           dgeco(a,lda,n,ipvt,rcond,z)
           if (!rcond is too small){
           	for (j=0,j<p,j++)
              		LP_dgesl(a,lda,n,ipvt,c[j][0],0);
	   }

     linpack. this version dated 08/14/78 .
     cleve moler, university of new mexico, argonne national lab.

     functions

     blas LP_daxpy,LP_ddot
*/
{
/*     internal variables	*/

	REAL t;
	int k,kb,l,nm1;

	nm1 = n - 1;
	if (job == 0) {

		/* job = 0 , solve  a * x = b
		   first solve  l*y = b    	*/

		if (nm1 >= 1) {
			for (k = 0; k < nm1; k++) {
				l = ipvt[k];
				t = b[l];
				if (l != k){ 
					b[l] = b[k];
					b[k] = t;
				}	
				LP_daxpy(n-(k+1),t,&a[lda*k+k+1],1,&b[k+1],1);
			}
		} 

		/* now solve  u*x = y */

		for (kb = 0; kb < n; kb++) {
		    k = n - (kb + 1);
		    b[k] = b[k]/a[lda*k+k];
		    t = -b[k];
		    LP_daxpy(k,t,&a[lda*k+0],1,&b[0],1);
		}
	}
	else { 

		/* job = nonzero, solve  trans(a) * x = b
		   first solve  trans(u)*y = b 			*/

		for (k = 0; k < n; k++) {
			t = LP_ddot(k,&a[lda*k+0],1,&b[0],1);
			b[k] = (b[k] - t)/a[lda*k+k];
		}

		/* now solve trans(l)*x = y	*/

		if (nm1 >= 1) {
			for (kb = 1; kb < nm1; kb++) {
				k = n - (kb+1);
				b[k] = b[k] + LP_ddot(n-(k+1),&a[lda*k+k+1],1,&b[k+1],1);
				l = ipvt[k];
				if (l != k) {
					t = b[l];
					b[l] = b[k];
					b[k] = t;
				}
			}
		}
	}
}




/**
 * solve the set of linear equations 
 *
 *     AX = B
 *
 * with possibly multiple rhs stored as columns of matrix B
 * the matrix A is not destroyed
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_lu_solve_seq = pnga_lu_solve_seq
#endif
void pnga_lu_solve_seq(char *trans, Integer g_a, Integer g_b) {

  logical oactive;  /* true iff this process participates */
  Integer dimA1, dimA2, typeA;
  Integer dimB1, dimB2, typeB;
  Integer me;
#if HAVE_LAPACK || ENABLE_F77
  BlasInt blas_dimA1, blas_dimA2, blas_dimB1, blas_dimB2, info=0;
#else
  Integer info=0;
#endif
  Integer dims[2], ndim;
  Integer lo[2], hi[2];

  /** check environment */
  me     = pnga_nodeid();
  
  /** check GA info for input arrays */
  pnga_check_handle(g_a, "ga_lu_solve: a");
  pnga_check_handle(g_b, "ga_lu_solve: b");
  pnga_inquire(g_a, &typeA, &ndim, dims);
  dimA1 = dims[0];
  dimA2 = dims[1];
  pnga_inquire(g_b, &typeB, &ndim, dims);
  dimB1 = dims[0];
  dimB2 = dims[1];
  
  GA_PUSH_NAME("ga_lu_solve_seq");

  if (dimA1 != dimA2) 
    pnga_error("ga_lu_solve: g_a must be square matrix ", 1);
  else if(dimA1 != dimB1) 
    pnga_error("ga_lu_solve: dims of A and B do not match ", 1);
  else if(typeA != C_DBL || typeB != C_DBL) 
    pnga_error("ga_lu_solve: wrong type(s) of A and/or B ", 1);
  
  pnga_sync();
  oactive = (me == 0);

  if (oactive) {
    DoublePrecision *adra, *adrb;
    Integer *adri;
    Integer one=1; 

    /** allocate a,b, and work and ipiv arrays */
    adra = (DoublePrecision*) ga_malloc(dimA1*dimA2, F_DBL, "a");
    adrb = (DoublePrecision*) ga_malloc(dimB1*dimB2, F_DBL, "b");
    adri = (Integer*) ga_malloc(GA_MIN(dimA1,dimA2), F_INT, "ipiv");

    /** Fill local arrays from global arrays */   
    lo[0] = one;
    hi[0] = dimA1;
    lo[1] = one;
    hi[1] = dimA2;
    pnga_get(g_a, lo, hi, adra, &dimA1);
    lo[0] = one;
    hi[0] = dimB1;
    lo[1] = one;
    hi[1] = dimB2;
    pnga_get(g_b, lo, hi, adrb, &dimB1);

    /** LU factorization */
#if HAVE_LAPACK || ENABLE_F77
    blas_dimA1 = dimA1;
    blas_dimA2 = dimA2;
    blas_dimB1 = dimB1;
    blas_dimB2 = dimB2;
    LAPACK_DGETRF(&blas_dimA1, &blas_dimA2, adra, &blas_dimA1, adri, &info);
#else
    {  int info_t;
      LP_dgefa(adra, (int)dimA1, (int)dimA2, (int*)adri, &info_t);
      info = info_t;
    }
#endif

    /** SOLVE */
    if(info == 0) {
#if HAVE_LAPACK || ENABLE_F77
      LAPACK_DGETRS(trans, &blas_dimA1, &blas_dimB2, adra, &blas_dimA1, 
          adri, adrb, &blas_dimB1, &info);
#else
      DoublePrecision *p_b;
      Integer i;
      int job=0;
      if(*trans == 't' || *trans == 'T') job = 1; 
      for(i=0; i<dimB2; i++) {
        p_b = adrb + i*dimB1;
        LP_dgesl(adra, (int)dimA1, (int)dimA2, (int*)adri, p_b, job);
      }
#endif

      if(info == 0) {
        lo[0] = one;
        hi[0] = dimB1;
        lo[1] = one;
        hi[1] = dimB2;
        pnga_put(g_b, lo, hi, adrb, &dimB1);
      } else
        pnga_error(" ga_lu_solve: LP_dgesl failed ", -info);

    }
    else
      pnga_error(" ga_lu_solve: LP_dgefa failed ", -info);

    /** deallocate work arrays */
    ga_free(adri);
    ga_free(adrb);
    ga_free(adra);
  }

  pnga_sync();
  
  GA_POP_NAME;
}
