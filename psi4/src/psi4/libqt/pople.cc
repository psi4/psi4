/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
  \file
  \brief Pople's method for solving linear equations
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "qt.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

#define ZERO 1e-13

/*!
** POPLE(): Uses Pople's method to iteratively solve linear equations
**          Ax = b
**
** Matt Leininger, April 1998
**
**  \param A         = matrix
**  \param x         = initially has vector b, but returns vector x.
**  \param dimen     = dimension of vector x.
**  \param num_vecs  = number of vectors x to obtain.
**  \param tolerance = cutoff threshold for norm of expansion vector.
**
** Returns: 0 for success, 1 for failure
** \ingroup QT
*/
int pople(double **A, double *x, int dimen, int /*num_vecs*/, double tolerance,
          std::string out, int print_lvl)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   double det, tval;
   double **Bmat; /* Matrix of expansion vectors */
   double **Ab;   /* Matrix of A x expansion vectors */
   double **M;    /* Pople subspace matrix */
   double *n;     /* overlap of b or transformed b in Ax = b
                     with expansion vector zero */
   double *r;     /* residual vector */
   double *b;     /* b vector in Ax = b, or transformed b vector */
   double *sign;  /* sign array  to insure diagonal element of A are positive */
   double **Mtmp; /* tmp M matrix passed to flin */
   int i, j, L=0, I;
   double norm, rnorm=1.0, *dotprod, *alpha;
   double *dvec;
   int llast=0, last=0, maxdimen;
   int quit=0;

   maxdimen = 200;
   /* initialize working array */
   Bmat = block_matrix(maxdimen, dimen);
   alpha = init_array(maxdimen);
   M = block_matrix(maxdimen, maxdimen);
   Mtmp = block_matrix(maxdimen, maxdimen);
   dvec = init_array(dimen);
   Ab = block_matrix(maxdimen, dimen);
   r = init_array(dimen);
   b = init_array(dimen);
   n = init_array(dimen);
   dotprod = init_array(maxdimen);
   sign = init_array(dimen);

   if (print_lvl > 6) {
   printer->Printf("\n\n Using Pople's Method for solving linear equations.\n");
   printer->Printf("     --------------------------------------------------\n");
   printer->Printf("         Iter             Norm of Residual Vector      \n");
   printer->Printf("        ------           -------------------------     \n");
   }

   norm = 0.0;
   for (i=0; i<dimen; i++) norm += x[i]*x[i];
   if (norm < ZERO) quit = 1;

   if (!quit) {
       for (i=0; i<dimen; i++) {
           if (A[i][i] > ZERO) sign[i] = 1.0;
           else sign[i] = -1.0;
           x[i] *= sign[i];
         }

       for (i=0; i<dimen; i++) {
           Bmat[0][i] = x[i];
           b[i] = x[i];
           dvec[i] = sqrt(fabs(A[i][i]));
           b[i] /= dvec[i];
           /*   outfile->Printf("A[%d][%d] = %lf\n",i,i, A[i][i]);
                outfile->Printf("dvec[%d] = %lf\n",i, dvec[i]);
                outfile->Printf("x[%d] = %lf\n",i, x[i]);
           */
         }

       if (print_lvl > 8) {
           printer->Printf(" A matrix in POPLE(LIBQT):\n");
           print_mat(A, dimen, dimen, out);
         }

       /* Constructing P matrix */
       for (i=0; i<dimen; i++) {
           for (j=0; j<dimen; j++) {
               if (i==j) A[i][i] = 0.0;
               else A[i][j] = -A[i][j]*sign[i];
             }
         }
       if (print_lvl > 8) {
           printer->Printf(" P matrix in POPLE(LIBQT):\n");
           print_mat(A, dimen, dimen, out);
         }

       /* Precondition P matrix with diagonal elements of A */
       for (i=0; i<dimen; i++)
           for (j=0; j<dimen; j++) {
               tval = 1.0/(dvec[j]*dvec[i]);
               A[i][j] *= tval;
             }

       if (print_lvl > 8) {
           printer->Printf(" Preconditioned P matrix in POPLE(LIBQT):\n");
           print_mat(A, dimen, dimen, out);
         }

       /*
          for (i=0; i<dimen; i++) {
          for (j=0; j<dimen; j++) {
          if (i==j) A[i][i] = 1.0 - A[i][i];
          else A[i][j] = -A[i][j];
          }
          }
       */

       for (i=0; i<dimen; i++) Bmat[0][i] /= dvec[i];


       dot_arr(Bmat[0],Bmat[0],dimen,&norm);
       norm = sqrt(norm);
       for (i=0; i<dimen; i++) {
           x[i] = Bmat[0][i];
           Bmat[0][i] /= norm;
         }

       dot_arr(Bmat[0],x,dimen,&(n[0]));

       while (!last) {
           /* form A*b_i */
           for (i=0; i<dimen; i++)
               dot_arr(A[i], Bmat[L], dimen, &(Ab[L][i]));


           /* Construct M matrix */
           /* M = delta_ij - <b_new|Ab_j> */
           zero_mat(M, maxdimen, maxdimen);
           for (i=0; i<=L; i++) {
               for (j=0; j<=L; j++) {
                   dot_arr(Bmat[i], Ab[j], dimen, &(dotprod[i]));
                   if (i==j) M[i][j] = 1.0 - dotprod[i];
                   else M[i][j] = -dotprod[i];
                 }
             }

           if (llast) last = llast;

           for (i=0; i<=L; i++) {
               alpha[i] = n[i];
               for (j=0; j<=L; j++)
                   Mtmp[i][j] = M[i][j];
             }

           flin(Mtmp, alpha, L+1, 1, &det);


           /* Need to form and backtransform x to orig basis x = D^(-1/2) x */
           zero_arr(x, dimen);
           for (i=0; i<=L; i++)
               for (I=0; I<dimen; I++)
                   x[I] += alpha[i]*Bmat[i][I]/dvec[I];

           /* Form residual vector Ax - b = r */
           zero_arr(r, dimen);
           for (I=0; I<dimen; I++)
               for (i=0; i<=L; i++)
                   r[I] += alpha[i]*(Bmat[i][I] - Ab[i][I]);

           for (I=0; I<dimen; I++) r[I] -= b[I];

           dot_arr(r, r, dimen, &rnorm);
           rnorm = sqrt(rnorm);
           if (print_lvl > 6) {
               printer->Printf(
                 "        %3d                     %10.3E\n",L+1,rnorm);

             }

           if (L+1>dimen) {
               printer->Printf("POPLE: Too many vectors in expansion space.\n");
               return 1;
             }

           /* place residual in b vector space */
           if (L+1>= maxdimen) {
               printer->Printf(
                 "POPLE (LIBQT): Number of expansion vectors exceeds"
                       " maxdimen (%d)\n", L+1);
               return 1;
             }
           for (i=0; i<dimen; i++) Bmat[L+1][i] = Ab[L][i];
           /*
              for (i=0; i<dimen; i++) Bmat[L+1][i] = r[i];
              for (i=0; i<dimen; i++) Bmat[L+1][i] = r[i]/(dvec[i]*dvec[i]);
           */


           /* Schmidt orthonormalize new b vec to expansion space */
           for (j=0; j<=L; j++) {
               dot_arr(Bmat[j], Bmat[L+1], dimen, &(dotprod[j]));
               for (I=0; I<dimen; I++)
                   Bmat[L+1][I] -= dotprod[j] * Bmat[j][I];
             }

           /* Normalize new expansion vector */
           dot_arr(Bmat[L+1], Bmat[L+1], dimen, &norm);
           norm = sqrt(norm);
           for (I=0; I<dimen; I++) Bmat[L+1][I] /= norm;

           /* check orthogonality of subspace */
           if (0) {
               for (i=0; i<=L+1; i++) {
                   for (j=0; j<=i; j++) {
                       dot_arr(Bmat[i], Bmat[j], dimen, &tval);
                       printer->Printf( "Bvec[%d] * Bvec[%d] = %f\n",i,j,tval);
                     }
                 }
             }


           if (rnorm<tolerance) llast = last = 1;
           L++;
         }

       zero_arr(x, dimen);
       for (i=0; i<=L; i++)
           for (I=0; I<dimen; I++)
               x[I] += alpha[i]*Bmat[i][I];

       /* Need to backtransform x to orginal basis x = D^(-1/2) x */
       for (I=0; I<dimen; I++)
           x[I] /= dvec[I];

     }
   free_block(Bmat);
   free(alpha);
   free_block(M);
   free_block(Mtmp);
   free(n);
   free(dvec);
   free_block(Ab);
   free(r);
   free(b);
   free(dotprod);

   return 0;
}

}
