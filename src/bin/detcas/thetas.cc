/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "globaldefs.h"
#include "globals.h"

namespace psi { namespace detcas { 

void rotate_test(int dim, int npairs, int *p_arr, int *q_arr, 
                 double *theta_arr);


/*
** postmult_by_U()
** 
** Postmultiply a block matrix (must be contiguous memory!) by a
**  unitary matrix U parameterized as a series of Givens rotations
**
*/
void postmult_by_U(int irrep, int dim, double **mat,
                   int npairs, int *p_arr, int *q_arr, double *theta_arr)
{

  int p, q, i, j, pair;
  double theta, sintheta, costheta;

  /* apply the transformation C = C^0 U, assuming square matrices */
  /* U is a series of Givens rotation matrices */
  for (pair=0; pair<npairs; pair++) {
    p = *p_arr++; 
    q = *q_arr++;
    theta = *theta_arr++;
    theta = -theta; /* the DROT routine rotates around the other way    */
                    /* compared to our def of G(theta) in the VBD paper */
    costheta = cos(theta);
    sintheta = sin(theta);
    C_DROT(dim,&(mat[0][q]),dim,&(mat[0][p]),dim,costheta,sintheta);
  }
}


/*
** postmult_by_exp_R()
** 
** Postmultiply a block matrix (must be contiguous memory!) by a
**  unitary matrix U parameterized as e^R with R an antisymmetric
**  matrix (based on formalism of notes by Yukio Yamaguchi)
**
*/
void postmult_by_exp_R(int irrep, int dim, double **mat,
  int npairs, int *p_arr, int *q_arr, double *theta_arr)
{

  int p, q, pair;
  double **R, **U, **C0, tval;

  U = block_matrix(dim,dim);
  R = block_matrix(dim,dim);
  C0= block_matrix(dim,dim);

  /* hm, warning, can't handle nbfso != nmo */
  for (p=0; p<dim; p++) {
    for (q=0; q<dim; q++) {
      C0[p][q] = mat[p][q];
    }
  }

  for (pair=0; pair<npairs; pair++) {
    p = *p_arr++; 
    q = *q_arr++;
    tval = *theta_arr++;
    R[p][q] = tval;
    U[p][q] = tval;
    R[q][p] = -tval;
    U[q][p] = -tval;
  }

  /* U = U + 1/2 R*R */
  /* can you pass the same matrix twice to DGEMM?? */
  C_DGEMM('n','n',dim,dim,dim,0.5,&(R[0][0]),dim,&(R[0][0]),dim,1.0,
          &(U[0][0]),dim);

  /* U += I */
  for (p=0; p<dim; p++) {
    U[p][p] += 1.0;
  }

  schmidt(U,dim,dim,outfile);

  /* C = C0 * U */
  C_DGEMM('n','n',dim,dim,dim,1.0,&(C0[0][0]),dim,&(U[0][0]),dim,0.0,
          &(mat[0][0]),dim);  
  
  free_block(C0);
  free_block(U);
  free_block(R);
}



/*
** premult_by_U()
** 
** Premultiply a block matrix (must be contiguous memory!) by a
**  unitary matrix U parameterized as a series of Givens rotations
**
*/
void premult_by_U(int irrep, int dim, double **mat,
                  int npairs, int *ppair, int *qpair, double *theta_arr)
{

  int p, q, i, j, pair;
  double theta, sintheta, costheta;

  /* apply the transformation X = U X, assuming square matrix */
  /* U is a series of Givens rotation matrices */
  /* n.b. we need to go backwards through the list of Givens rotations */
  for (pair=npairs-1; pair>=0; pair--) {
    p = ppair[pair]; 
    q = qpair[pair];
    theta = theta_arr[pair];
    costheta = cos(theta);
    sintheta = sin(theta);
    C_DROT(dim,&(mat[q][0]),1,&(mat[p][0]),1,costheta,sintheta);
  }

}




/*
** rotate_test
**
** This does about the same thing as rotate_orbs_irrep except that
** it does it on a unit matrix so the results can be checked easily
*/
void rotate_test(int dim, int npairs, int *p_arr, int *q_arr, 
                 double *theta_arr)
{

  int p, q, i, j, pair;
  double **tmpmat, theta, sintheta, costheta;

  /* set up a unit matrix */
  tmpmat = block_matrix(dim,dim);
  for (i=0; i<dim; i++) {
    tmpmat[i][i] = 1.0;
  }
 
  /* print new coefficients */
  fprintf(outfile, "\n\tOld molecular orbitals\n");
  print_mat(tmpmat, dim, dim, outfile);


  /* now apply U, as a series of Givens rotation matrices */
  for (pair=0; pair<npairs; pair++) {
    p = *p_arr++; 
    q = *q_arr++;
    theta = *theta_arr++;
    theta = -theta;  /* DROT rotates around the other way */
    costheta = cos(theta);
    sintheta = sin(theta);
    fprintf(outfile, "\nApplying rotation (%2d,%2d) = %12.6lf\n", p, q, theta);
    fprintf(outfile, "Cos(theta)=%12.6lf, Sin(theta)=%12.6lf\n", 
            costheta, sintheta);
    C_DROT(dim,&(tmpmat[0][q]),dim,&(tmpmat[0][p]),dim,costheta,sintheta);
    if (Params.print_lvl > 3) {
      fprintf(outfile, "\n\tMatrix after transformation:\n");
      print_mat(tmpmat, dim, dim, outfile);
    }
  }

  /* print new coefficients */
  fprintf(outfile, "\n\tNew molecular orbitals\n");
  print_mat(tmpmat, dim, dim, outfile);

  free_block(tmpmat);

}


/*
** read_thetas()
**
** Read in the theta array from disk.  If there is none, assume they're
**  all set to 0.
*/
void read_thetas(int npairs)
{

  FILE *fp;

  CalcInfo.theta_cur = init_array(npairs);

  /* look for the thetas on disk...if they're around, read them in */
  ffileb_noexit(&fp,"thetas.dat",2);
  if (fp != NULL) {
    if (Params.print_lvl > 2)
      fprintf(outfile, "\nReading orbital rotation angles\n");
    if (fread(CalcInfo.theta_cur, sizeof(double), npairs, fp) != npairs) {
      fprintf(outfile, "Error reading angles.\n");
      zero_arr(CalcInfo.theta_cur, npairs);
    }
    fclose(fp);
  }

    
}



/*
** write_thetas()
**
** Write the theta array to disk. 
*/
void write_thetas(int npairs)
{

  FILE *fp;

  ffileb_noexit(&fp,"thetas.dat",0);
  if (fp != NULL) {
    if (Params.print_lvl > 2)
      fprintf(outfile, "\nWriting orbital rotation angles\n");
    if (fwrite(CalcInfo.theta_cur, sizeof(double), npairs, fp) != npairs) {
      fprintf(outfile, "Error writing angles.\n");
    }
    fclose(fp);
  }

  else {
    fprintf(outfile, "Error opening thetas.dat for writing\n");
  }
    
}


/*
** calc_de_dtheta()
**
** This function calculates dE / dTheta = dE/dU * dU/dTheta
**
** Feed in a vector of thetas (theta) and a matrix dE/dU and write
**   to a vector of dE/d(theta)  (dET).  Actually, this is called
**   for an irrep at a time.
**
** Based on similar code from an old version of the VBD code.
** Can use blas calls to drot later.
**
** C. David Sherrill
** May 1998
**
*/
void calc_dE_dT(int n, double **dEU, int npairs, int *ppair, int *qpair,
                double *theta, double *dET)
{
  int i,a,m,l,pair;
  double temp, costheta, sintheta;
  double **Uleft, **Uright, **Scratch;

  zero_arr(dET, npairs);


  /* form temporary U matrices */

  Uleft = block_matrix(n ,n);
  Uright = block_matrix(n, n);
  Scratch = block_matrix(n, n);

  /* init Uright and Uleft to unit matrix */
  for(m=0; m < n; m++)  {
    Uright[m][m] = 1.0;
    Uleft[m][m] = 1.0;
  }

  /* init Uleft to the U matrix */

  for (pair=0; pair<npairs; pair++) {
    a = ppair[pair];
    i = qpair[pair]; 

    costheta = cos(theta[pair]);
    sintheta = sin(theta[pair]);

    /* in-place postmultication of Uleft by G_{ia} */
    for(m=0; m < n; m++)  {
      temp = Uleft[m][i];
      Uleft[m][i] = temp*costheta - Uleft[m][a]*sintheta;
      Uleft[m][a] = Uleft[m][a]*costheta + temp*sintheta;
      }
  }

  if (Params.print_lvl > 3) {
    fprintf(outfile, "dE/dU after backtransform: \n");
    print_mat(dEU, n, n, outfile);
  }

  /* Loop over i,a to form dE/d(theta): we are working right to left in this
     algorithm, hence we must go backwards through the theta list */
  
  for (pair=npairs-1; pair>=0; pair--) {
    a = ppair[pair];
    i = qpair[pair];

    costheta = cos(theta[pair]);
    sintheta = sin(theta[pair]);

    /*
    fprintf(outfile, "Derivative (i=%d, a=%d)\n", i, a);
    fprintf(outfile, "Cos = %lf, Sin=%lf\n", costheta, sintheta);
    */

    /* post-multiply Uleft by G(+) */
    for(m=0; m < n; m++)  {
      temp = Uleft[m][i];
      Uleft[m][i] = temp*costheta + Uleft[m][a]*sintheta;
      Uleft[m][a] = Uleft[m][a]*costheta - temp*sintheta;
    }

    /*
    fprintf(outfile, "Uleft after postmultiplication by G(+)(%d,%d)\n",
      i, a);
    print_mat(Uleft, n, n, outfile);
    */

    /* Now do Uleft*dG/d(theta)*Uright series of multi */
    for(l=0; l < n; l++)  {
      for(m=0; m < n; m++)  {
        Scratch[l][m] = (-sintheta*Uleft[l][i] - costheta*Uleft[l][a])
          *Uright[i][m] + (costheta*Uleft[l][i] - sintheta*Uleft[l][a])
            *Uright[a][m];
      }
    }

    /*
    fprintf(outfile, "Uleft * dG/dTheta(%d,%d) * Uright\n", i, a);
    print_mat(Scratch, n, n, outfile);
    */

    for(l=0; l < n; l++)  {
      for(m=0; m < n; m++)  {
        dET[pair] += dEU[l][m]*Scratch[l][m];
        /*
        if ((fabs(Scratch[l][m]) > 0.0001) && ((l<nocc && m<nocc) ||
          (l >= nocc && m >= nocc))) {
          fprintf(outfile, "nonzero element for theta(%d, %d) element", i, a);
          fprintf(outfile, "%d %d\n", l, m);
        }
        */
      }
    }

    /*
    fprintf(outfile, "dE/dTheta(%d,%d) = %12.6lf\n", i, a, dET[pair]);
    */

    /* pre-multiply Uright by G */
    for(m=0; m < n; m++)  {
      temp = Uright[i][m];
      Uright[i][m] = temp*costheta + Uright[a][m]*sintheta;
      Uright[a][m] = Uright[a][m]*costheta - temp*sintheta;
    }

    /*
    fprintf(outfile, "Uright after premultiplication by G \n");
    print_mat(Uright, n, n, outfile);
    */

  }

  /* free memory */
  free_block(Uleft);
  free_block(Uright);
  free_block(Scratch);
}

}} // end namespace psi::detcas

