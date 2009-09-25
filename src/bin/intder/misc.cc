/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "3dmatrix.h"


#define EXTERN
#include "globals.h"
#include "params.h"

using namespace psi::intder;

namespace psi { namespace intder {
extern Params gParams;
}}

#include <libciomr/libciomr.h>
#include <physconst.h>

namespace psi { namespace intder {

//Calculates the dot product of two three-dimensional vectors
//SCAPRO(U, V, D) from intder2000.f
double dot_prod(double *e2,double *e1)
{
  int i;
  double dotp = 0;

  for(i = 0; i < 3; i++)
    dotp += (e1[i] * e2[i]);

  return dotp;
}

//Returns vector product of two 3 dimensional vectors
double *vect_prod(double *m1, double *m2)
{
  double *cprod = (double *)malloc(3 * sizeof(double));

  cprod[0] = m1[1]*m2[2] - m1[2]*m2[1];
  cprod[1] = -(m1[0]*m2[2] - m1[2]*m2[0]);
  cprod[2] = m1[0]*m2[1] - m1[1]*m2[0];

  return cprod;
}

//MAT1 from intder2000.f, see README
void mat1(double **EM, double *V)
{
	int i;

	for(i=0;i<3;i++)
		EM[i][i] = 0.0;

	EM[1][0] = -V[2];
	EM[2][0] = V[1];
	EM[2][1] = -V[0];
	EM[0][1] = -EM[1][0];
	EM[0][2] = -EM[2][0];
	EM[1][2] = -EM[2][1];
}

//MAT2 from intder2000.f
void mat2(double **AA, double *VK)
{
	int i;

	for(i=0;i<3;i++)
		AA[i][i] = 0.0;

	AA[0][1] = VK[0];
	AA[0][2] = VK[1];
	AA[1][2] = VK[2];
	AA[1][0] = -AA[0][1];
	AA[2][0] = -AA[0][2];
	AA[2][1] = -AA[1][2];

}

//EXPMAT from intder2000.f, takes the exponent of a 3x3 matrix by expanding over a power series
double **exp_matrix(double *AA)
{
  int i, j, l;
  int iter = 0;
  int MAXITER = 50;  //Maximum orders of power series
  double tolerance = 1.0E-50;
  double norm = 0.0;

  double **kmat;
  double **tempmat1;
  double **tempmat2;
  double **tempmat3;
  
  kmat = init_matrix(3,3);
  tempmat1 = init_matrix(3,3);
  tempmat2 = init_matrix(3,3);
  tempmat3 = init_matrix(3,3);

  for(i = 0; i < 3; i++) {
    tempmat1[i][i]=1.0;
    tempmat2[i][i]=1.0;
  }
  
  kmat[0][1] = AA[0];
  kmat[0][2] = AA[1];
  kmat[1][2] = AA[2];
  kmat[1][0] = -(AA[0]);
  kmat[2][0] = -(AA[1]);
  kmat[2][1] = -(AA[2]);

  do {
    iter++;
    tempmat3 = init_matrix(3,3);
    
    for(i = 0; i < 3; i++) 
      {
	for(j = 0; j < 3; j++) {
	  for(l = 0; l < 3; l++) {
	    tempmat3[i][j] += tempmat1[i][l] * kmat[l][j];
	    tempmat3[i][j] /= (double)iter;
	    if(fabs(tempmat3[i][j]) > norm) 
	      norm = fabs(tempmat3[i][j]);
	  }
	}
      }
    
    for(i = 0; i < 3; i++) 
      {
	for(j = 0; j < 3; j++) {   
	  tempmat2[i][j] += tempmat3[i][j];
	  tempmat1[i][j] = tempmat3[i][j];
	}
      }
	  
  } while(iter < MAXITER && norm > tolerance);
  
  if(norm > tolerance)
    fprintf(outfile, "\n!!! Convergence not reached in exp_matrix after %i iterations. CNORM = %lf  !!!\n", iter, norm);
  
  free_matrix(tempmat1,3);
  free_matrix(tempmat3,3);

  return tempmat2;
}


//TRIPRO from intder2000.f, returns triple products of unit vectors (a pointer to a 3x3x3 3dmatrix)
C3DMatrix *tripro(void)
{
	int i, j, k;
	C3DMatrix *product = new C3DMatrix(3, 3, 3);
	double *vect;
	double **rmat;

	vect = new double[3];
	rmat = init_matrix(3,3);

	for(k=0;k<3;k++)
	{
		for(i=0;i<3;i++)
			vect[i] = 0.0;

		vect[k] = 1.0;
		mat1(rmat, vect);
		for(j=0;j<3;j++)
			for(i=0;i<3;i++)
				product->Set(i, j, k, rmat[i][j]);
	}
		
	delete[] vect;
	free_matrix(rmat,3);
	return product;
}	

void fill3a(int nx, int ny, C3DMatrix *f3)
{
	int m, n, p;

	for(p=0; p<ny; p++)	{
		for(n=0; n<p-1; n++)	{
			for(m=0; m<n-1; m++)	{
				f3->Set(n, m, p, f3->Get(m, n, p));
				f3->Set(n, p, m, f3->Get(m, n, p));
				f3->Set(m, p, n, f3->Get(m, n, p));
				f3->Set(p, m, n, f3->Get(m, n, p));
				f3->Set(p, n, m, f3->Get(m, n, p));
			}
			f3->Set(n, p, n, f3->Get(m, n, p));
			f3->Set(p, n, n, f3->Get(m, n, p));
		}
		for(m=0; m<p-1; m++)	{
			f3->Set(p, m, p, f3->Get(m, p, p));
			f3->Set(p, p, m, f3->Get(m, p, p));
		}
	}
}

void fill3b(int nx, int ny, C3DMatrix *f3)
{
	int m, n, p;

	for(m=0; m<ny; m++) {
		for(n=0; n<m-1; n++) {
			for(p=0; p<n-1; p++) {
				f3->Set(n, m, p, f3->Get(m, n, p));
				f3->Set(n, p, m, f3->Get(m, n, p));
				f3->Set(m, p, n, f3->Get(m, n, p));
				f3->Set(p, m, n, f3->Get(m, n, p));
				f3->Set(p, n, m, f3->Get(m, n, p));
			}
			f3->Set(n, m, n, f3->Get(m, n, n));
			f3->Set(n, n, m, f3->Get(m, n, n));
		}

		for(p=0; p<m-1; p++) {
			f3->Set(m, p, m, f3->Get(m, m, p));
			f3->Set(p, m, m, f3->Get(m, m, p));
		}
	}
}

void fill4a(int nx, int ny, C4DMatrix *f4)
{
  int m, n, p, q;

  for (q=0; q<ny; q++) {
    for (p=0; p<=q; p++) {
      for (n=0; n<=p; n++) {
        for (m=0; m<=n; m++) {
          f4->Set(n,m,p,q, f4->Get(m,n,p,q));
          f4->Set(n,p,m,q, f4->Get(m,n,p,q));
          f4->Set(n,p,q,m, f4->Get(m,n,p,q));
          f4->Set(m,p,n,q, f4->Get(m,n,p,q));
          f4->Set(p,m,n,q, f4->Get(m,n,p,q));
          f4->Set(p,n,m,q, f4->Get(m,n,p,q));
          f4->Set(p,n,q,m, f4->Get(m,n,p,q));
          f4->Set(m,p,q,n, f4->Get(m,n,p,q));
          f4->Set(p,m,q,n, f4->Get(m,n,p,q));
          f4->Set(p,q,m,n, f4->Get(m,n,p,q));
          f4->Set(p,q,n,m, f4->Get(m,n,p,q));
          f4->Set(m,n,q,p, f4->Get(m,n,p,q));
          f4->Set(n,m,q,p, f4->Get(m,n,p,q));
          f4->Set(n,q,m,p, f4->Get(m,n,p,q));
          f4->Set(n,q,p,m, f4->Get(m,n,p,q));
          f4->Set(m,q,n,p, f4->Get(m,n,p,q));
          f4->Set(q,m,n,p, f4->Get(m,n,p,q));
          f4->Set(q,n,m,p, f4->Get(m,n,p,q));
          f4->Set(q,n,p,m, f4->Get(m,n,p,q));
          f4->Set(m,q,p,n, f4->Get(m,n,p,q));
          f4->Set(q,m,p,n, f4->Get(m,n,p,q));
          f4->Set(q,p,m,n, f4->Get(m,n,p,q));
          f4->Set(q,p,n,m, f4->Get(m,n,p,q));
        }
      }
    }
  }
}

//'Fill4A' for C5DMatrix
void fill4a(int nx, int ny, C5DMatrix *f4)
{
  int m, n, p, q, r;

	for(r=0; r<nx; r++) {
	  for (q=0; q<ny; q++) {
  	  for (p=0; p<=q; p++) {
    	  for (n=0; n<=p; n++) {
      	  for (m=0; m<=n; m++) {
       	    f4->Set(n,m,p,q,r, f4->Get(m,n,p,q,r));
         		f4->Set(n,p,m,q,r, f4->Get(m,n,p,q,r));
         		f4->Set(n,p,q,m,r, f4->Get(m,n,p,q,r));
         		f4->Set(m,p,n,q,r, f4->Get(m,n,p,q,r));
         		f4->Set(p,m,n,q,r, f4->Get(m,n,p,q,r));
           	f4->Set(p,n,m,q,r, f4->Get(m,n,p,q,r));
           	f4->Set(p,n,q,m,r, f4->Get(m,n,p,q,r));
            f4->Set(m,p,q,n,r, f4->Get(m,n,p,q,r));
            f4->Set(p,m,q,n,r, f4->Get(m,n,p,q,r));
            f4->Set(p,q,m,n,r, f4->Get(m,n,p,q,r));
            f4->Set(p,q,n,m,r, f4->Get(m,n,p,q,r));
            f4->Set(m,n,q,p,r, f4->Get(m,n,p,q,r));
            f4->Set(n,m,q,p,r, f4->Get(m,n,p,q,r));
            f4->Set(n,q,m,p,r, f4->Get(m,n,p,q,r));
            f4->Set(n,q,p,m,r, f4->Get(m,n,p,q,r));
            f4->Set(m,q,n,p,r, f4->Get(m,n,p,q,r));
            f4->Set(q,m,n,p,r, f4->Get(m,n,p,q,r));
            f4->Set(q,n,m,p,r, f4->Get(m,n,p,q,r));
            f4->Set(q,n,p,m,r, f4->Get(m,n,p,q,r));
            f4->Set(m,q,p,n,r, f4->Get(m,n,p,q,r));
            f4->Set(q,m,p,n,r, f4->Get(m,n,p,q,r));
            f4->Set(q,p,m,n,r, f4->Get(m,n,p,q,r));
            f4->Set(q,p,n,m,r, f4->Get(m,n,p,q,r));
					}
        }
      }
    }
  }
}


/*** SYM_MATRIX_INVERT inverts a matrix by diagonalization
 *  **A = matrix to be inverted
 *  dim = dimension of A
 *  print_det = print determinant if 1, nothing if 0
 *  redundant = zero eigenvalues allowed if 1
 *  returns: inverse of A          written by Rollin A King */

double **symm_matrix_invert(double **A, int dim, int print_det, int redundant) 
{
  int i;
  double **A_inv, **A_vects, *A_vals, **A_temp, det=1.0;
  
  A_inv   = block_matrix(dim,dim);
  A_temp  = block_matrix(dim,dim);
  A_vects = block_matrix(dim,dim);
  A_vals  = init_array(dim);
  
  sq_rsp(dim,dim,A,A_vals,1,A_vects,EVAL_TOL);
  
  if (redundant == 0) {
    for (i=0;i<dim;++i) {
      det *= A_vals[i];
      A_inv[i][i] = 1.0/A_vals[i];
    }
    if (print_det)
      fprintf(outfile,"Determinant: %10.6e\n",det);
    if (fabs(det) < 1E-10) {
      fprintf(outfile,"Determinant: %10.6e\n",det);
      fprintf(outfile,"Determinant is too small...aborting.\n");
      fclose(outfile);
      exit(2);
    }
  }
  else {
    for (i=0;i<dim;++i) {
      det *= A_vals[i];
       if (fabs(A_vals[i]) > REDUNDANT_EVAL_TOL)
        A_inv[i][i] = 1.0/A_vals[i];
      else
        A_inv[i][i] = 0.0;
    }
    if (print_det)
      fprintf(outfile,"Determinant: %10.6e\n",det);
  }

  mmult(A_inv,0,A_vects,1,A_temp,0,dim,dim,dim,0);
  mmult(A_vects,0,A_temp,0,A_inv,0,dim,dim,dim,0);

  free(A_vals);
  free_block(A_vects);
  free_block(A_temp);
  return A_inv;
}

//What on Earth is the point of this function? 
double dot_x(double val1, int idim1, double val2, int idim2, int n)
{
  int i = 0;
  double d = 1;

  for(i = 0; i < n; i++) {
    d += (val1 * val2);
    val1 += idim1;
    val2 += idim2;
  }
      
  return d;
}

//this is now replaced with mat1, mat2, tripro
/*
double ***Christoffel1(void)
  //Function for Christoffel symbol of first kind, used to differentiate vectors, known as
  //TRIPRO / MAT1 / MAT2 functions of intder2000. 
  //What this does is return a 3 x 3 x 3 matrix (or an array of matrices) whose values
  //are the triple products of unit vectors (meaning they have values of 0 or 1)

  /* {
  int i,j,k;
  static double ***Chris = NULL;
  double *evect = NULL;

  if (Chris == NULL)
    {
      Chris = new double[3][3];
      evect = new double[3];

      for(i=0; i < 3; i++)
	for (j=0; j < 3; j++)
	  Chris[i][j] = (double*)new double[3];

      for(i = 0; i < 3; i++) 
	for(j = 0; j < 3; j++) {
	  evect[j] = ((i == j) ? 1.0 : 0.0);
	    for(k = 0; k < 3; k++) {
	      Chris[1][0][i] = -evect[2];
	      Chris[0][1][i] = -Chris[0][1][i];
	      Chris[2][0][i] = evect[1];
	      Chris[0][2][i] = -Chris[2][0][i];
	      Chris[1][2][i] = -evect[0];
	      Chris[2][1][i] = -Chris[2][1][i];
	    }
	}
      delete evect;
    }
    return Chris;
  }
*/

}} // namespace psi::intder
