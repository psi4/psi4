/*! \file
    \ingroup OPTKING
    \brief miscellaneous little matrix and print functions
*/

#define EXTERN
#include "globals.h"
#undef EXTERN

#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>

namespace psi { namespace optking {

double nuclear_repulsion(double *fatomic_num, double *coord) {
  int i, j, dim;
  double dist, tval = 0.0;

  dim = optinfo.natom;
  for (i=0; i<dim; ++i)
    for (j=0; j<i; ++j) {
      dist = sqrt(
          SQR(coord[3*i+0]-coord[3*j+0])
          + SQR(coord[3*i+1]-coord[3*j+1])
          + SQR(coord[3*i+2]-coord[3*j+2]) );

      tval += fatomic_num[i]*fatomic_num[j] / dist;
    }
//  fprintf(outfile,"returning repulsion: %15.10lf \n", tval);
  return tval;
}

/*** PRINT_MAT2   prints a matrix to output file ***/
void print_mat2(double **matrix, int rows, int cols, FILE *of) {
  int i,j,col;
  for (i=0;i<rows;++i) {
    col = 0;
    for (j=0;j<cols;++j) {
      if (col == 9) {
        fprintf(outfile,"\n");
        col = 0;
      }
      fprintf(of,"%15.10lf",matrix[i][j]);
      ++col;
    }
    fprintf(outfile,"\n");
  } 
  return;
} 

/*** PRINT_MAT5   prints a matrix to output file ***/
void print_mat5(double **matrix, int rows, int cols, FILE *of) {
  int i,j,col;
  for (i=0;i<rows;++i) {
    col = 0;
    for (j=0;j<cols;++j) {
      if (col == 14) {
        fprintf(outfile,"\n");
        col = 0;
      }
      fprintf(of," %9.5lf",matrix[i][j]);
      ++col;
    }
    fprintf(outfile,"\n");
  } 
  return;
} 

/*** CROSS_PRODUCT   computes cross product of two vectors ***/
void cross_product(double *u,double *v,double *out)
{
  out[0] = u[1]*v[2]-u[2]*v[1];
  out[1] = -1.0*(u[0]*v[2]-u[2]*v[0]);
  out[2] = u[0]*v[1]-u[1]*v[0];
  return;
} 

/*** SCALAR_MULT   performs scalar multiplication of a vector***/ 
void scalar_mult(double a, double *vect, int dim) {
  int i;
  if (dim == 0) return;
  //if ( (fabs(a) < 1.0e-10) || (fabs(a) > 1.0e10) )
    //fprintf(outfile,"Warning: scalar_mult() scaling by %9.5e\n",a);
  for (i=0;i<dim;++i)
    vect[i] *= a;
  return;
}

/*** SCALAR_DIV   performs scalar division of a vector ***/
void scalar_div(double a, double *vect) {
  int i;
  for (i=0;i<3;++i)
    vect[i] /= a;
  return;
}

void punt(const char *message) {
  fprintf(outfile,"\nerror: %s\n", message);
  fprintf(outfile,"         *** stopping execution ***\n");
  fprintf(stderr,"\n OPTKING error: %s\n", message);
  fprintf(stderr,"                 *** stopping execution ***\n");
  exit_io();
}

void open_PSIF(void) {
  psio_open(PSIF_OPTKING, PSIO_OPEN_OLD);
  return;
}

void close_PSIF(void) {
  psio_close(PSIF_OPTKING, 1);
  return;
}

/*** SWAP_TORS -- canonical torsion order is atom a < atom d ***/
void swap_tors(int *a, int *b, int *c, int *d) {
  int p;
  if (*a > *d) {
    p = *a;
    *a = *d;
    *d = p;
    p = *b;
    *b = *c;
    *c = p;
  }
}

/*** SWAP -- canonical order of bonds is atom a < atom b ***/
void swap(int *a, int *b) {
  int c;
  if (*a > *b) {
    c = *a;
    *a = *b;
    *b = c;
  }
  return;
}

/*** DIV_INT Does little go into big? ***/
int div_int(int big, int little) {
  if (little > big) return 0;
  for(;;) {
    big -= little;
    if (big == 0) return 1;
    if (big < 0) return 0;
  }
}

// dot arrays of doubles
void dot_array(double *a, double *b, long int n, double *value) {
  long int i;
  double tval = 0.0;
  for (i=0; i < n; i++) tval += a[i]*b[i];
  *value = tval;
}

/*** SYM_MATRIX_INVERT inverts a matrix by diagonalization
 *  **A = matrix to be inverted
 *  dim = dimension of A
 *  print_det = print determinant if 1, nothing if 0
 *  redundant = zero eigenvalues allowed if 1
 *  returns: inverse of A ***/
double **symm_matrix_invert(double **A, int dim, int print_det, int redundant) {
  int i;
  double **A_inv, **A_vects, *A_vals, **A_temp, det=1.0;

  A_inv   = init_matrix(dim,dim);
  A_temp  = init_matrix(dim,dim);
  A_vects = init_matrix(dim,dim);
  A_vals  = init_array(dim);

  opt_sq_rsp(dim,dim,A,A_vals,1,A_vects,EVAL_TOL);

  if (redundant == 0) {
    for (i=0;i<dim;++i) {
      det *= A_vals[i];
      A_inv[i][i] = 1.0/A_vals[i];
    }
    if (print_det)
      fprintf(outfile,"Determinant: %10.6e\n",det);
    if (fabs(det) < 1E-10)
      fprintf(outfile,"Warning, determinant is small: %10.6e\n",det);
  }
  else {
    for (i=0;i<dim;++i) {
      det *= A_vals[i];
      if (fabs(A_vals[i]) > REDUNDANT_EVAL_TOL) {
        A_inv[i][i] = 1.0/A_vals[i];
      }
      else {
fprintf(outfile,"\ndetected redundant eval - setting inverse element to 0 \n");
        A_inv[i][i] = 0.0;
      }
    }
    if (print_det)
      fprintf(outfile,"Determinant: %10.6e\n",det);
  }

  opt_mmult(A_inv,0,A_vects,1,A_temp,0,dim,dim,dim,0);
  opt_mmult(A_vects,0,A_temp,0,A_inv,0,dim,dim,dim,0);

  free_array(A_vals);
  free_matrix(A_vects);
  free_matrix(A_temp);
  return A_inv;
} 

double energy_chkpt(void) {
  double energy;

  chkpt_init(PSIO_OPEN_OLD);
  /* energy = chkpt_rd_escf(); */
  energy = chkpt_rd_etot();
  chkpt_close();
  return energy;
}
     

/* This routine transposes matrices and calls lapack dgeev() in libqt *
 * to diagonalize a square nonsymmetric matrix.  The eigenvalues      *
 * are returned in random order.                                      */
void dgeev_optking(int L, double **G, double *lambda, double **alpha) {

  int i, j, lwork, info;
  double *evals_i, *work, **left_evects, tval, temp;

  evals_i = init_array(L); 
  left_evects = init_matrix(L,L);

  work = init_array(20*L);
  lwork = 20*L;          

  for (i=0; i<L; ++i)
    for (j=0; j<i; ++j) {
      temp = G[i][j];
      G[i][j] = G[j][i];
      G[j][i] = temp;
    }

  i = C_DGEEV(L, G, L, lambda, evals_i, left_evects,
    L, alpha, L, work, lwork, info);

  for (i=0; i<L; ++i)
    for (j=0; j<i; ++j) {
      temp = alpha[i][j];
      alpha[i][j] = alpha[j][i];
      alpha[j][i] = temp;
    }

  free_array(work);
  free_array(evals_i);
  free_matrix(left_evects);
  return;
}

double **mass_mat(double *masses) {
    int i, dim;
    double **u;

    dim = 3*optinfo.natom;
    u = init_matrix(dim,dim);

    for (i=0; i<dim; ++i) {
      if (masses[i] == 0.0)
        u[i][i] = 0.0;
      else
        u[i][i] = 1.0/masses[i];
    }
    return u;
}

void print_evects(double **evects, double *evals, int nrow, int ncol, FILE *out)
{
      int ii,jj,kk,nn;
      int i,j;

      ii=0;jj=0;
L200:
      ii++;
      jj++;
      kk=10*jj;
      nn=ncol;
      if (nn > kk) nn=kk;
      fprintf (out,"\n");
      for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
      fprintf (out,"\n");
      for (i=0; i < nrow; i++) {
         fprintf (out,"\n%5d",i+1);
         for (j=ii-1; j < nn; j++) {
            fprintf (out,"%12.7f",evects[i][j]);
            }
         }
      fprintf (out,"\n");
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",evals[j]);
         }
      fprintf (out,"\n");
      if (ncol <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
}

void opt_report(FILE *of) {
  int i;
  double tval, tval2, tval3, tval4;
  char keyword[30];
  fprintf(of,"\n\t            ****  Optimization Summary  ****\n");
  fprintf(of,"\t----------------------------------------------------------------------\n");
  fprintf(of,"\t Step         Energy             Delta(E)      RMS force    MAX force \n");
  fprintf(of,"\t----------------------------------------------------------------------\n");
  open_PSIF();
  tval2 = 0;
  for (i=0; i<optinfo.iteration; ++i) {
    sprintf(keyword,"Energy %d", i);
    psio_read_entry(PSIF_OPTKING, keyword, (char *) &tval, sizeof(double));
    sprintf(keyword,"RMS force %d", i);
    psio_read_entry(PSIF_OPTKING, keyword, (char *) &tval3, sizeof(double));
    sprintf(keyword,"MAX force %d", i);
    psio_read_entry(PSIF_OPTKING, keyword, (char *) &tval4, sizeof(double));
    fprintf(of,"\t %3d  %18.12lf  %18.12lf  %10.2e   %10.2e\n", i+1, tval, tval - tval2, tval3 ,tval4);
    tval2 = tval;
  }
  fprintf(of,"\t----------------------------------------------------------------------\n");
  fprintf(of,"\n");
  close_PSIF();
}

void print_mat(double **a, int m, int n, FILE *out) {
  int ii=0,jj=0,kk=0,nn,ll;
  int i,j,k;

  while (n > kk)  {
    ii++;
    jj++;
    kk=10*jj;
    nn=n;
    if (nn > kk) nn=kk;
    ll = 2*(nn-ii+1)+1;
    fprintf (out,"\n");
    for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
    fprintf (out,"\n");
    for (i=0; i < m; i++) {
      fprintf (out,"\n%5d",i+1);
      for (j=ii-1; j < nn; j++) {
        fprintf (out,"%12.7f",a[i][j]);
      }
    }
    fprintf (out,"\n");
    ii=kk;
  }
  fflush(out);
  return;
}

void eivout(double **a, double *b, int m, int n, FILE *out) {
  int ii=0,jj=0,kk=0,nn;
  int i,j;

  while (n > kk) {
    ii++;
    jj++;
    kk=10*jj;
    nn=n;
    if (nn > kk) nn=kk;
    fprintf (out,"\n");
    for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
    fprintf (out,"\n");
    for (i=0; i < m; i++) {
      fprintf (out,"\n%5d",i+1);
      for (j=ii-1; j < nn; j++) {
        fprintf (out,"%12.7f",a[i][j]);
      }
    }
    fprintf (out,"\n");
    fprintf (out,"\n     ");
    for (j=ii-1; j < nn; j++) {
      fprintf(out,"%12.7f",b[j]);
    }
    fprintf (out,"\n");
    ii=kk;
  }
  fflush(out);
  return;
}








}} /* namespace psi::optking */

