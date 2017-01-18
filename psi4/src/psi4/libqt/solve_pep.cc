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
  \brief Solve a 2x2 pseudo-eigenvalue problem
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/psi4-dec.h"
namespace psi {

#define A_MIN 1.0E-10

/*!
** solve_2x2_pep(): Solve a 2x2 pseudo-eigenvalue problem of the form
**    [ H11 - E    H12 - E*S ]  [c1]
**    [ H12 - E*S  H22 - E   ]  [c2]  = 0
**
**  \param H     =  matrix to get eigenvalues of
**  \param S     =  overlap between states 1 & 2
**  \param evals =  pointer to array to hold 2 eigenvalues
**  \param evecs =  matrix to hold 2 eigenvectors
**
** Returns: none
** \ingroup QT
*/
void solve_2x2_pep(double **H, double S, double *evals, double **evecs)
{
   int i ;
   double a, b, c, p, q, x ;
   double norm, tval, tval2;

   /* put in quadratic form */
   a = 1.0 - S * S ;
   b = 2.0 * S * H[0][1] - H[0][0] - H[1][1] ;
   c = H[0][0] * H[1][1] - H[0][1] * H[0][1] ;

   /* solve the quadratic equation for E0 and E1 */

   tval = b*b ;
   tval -= 4.0 * a * c ;
   tval2 = sqrt(tval) ;
   if (tval < 0.0) {
      outfile->Printf( "(solve_2x2_pep): radical less than 0\n") ;
      return ;
      }
   else if (fabs(a) < A_MIN) {
      printf("min a reached\n") ;
      evals[0] = evals[1] = H[1][1] ;
      }
   else {
      evals[0] = -b / (2.0 * a) ;
      evals[0] -= tval2 / (2.0 * a) ;
      evals[1] = -b / (2.0 * a) ;
      evals[1] += tval2 / (2.0 * a) ;
      }

   /* Make sure evals[0] < evals[1] */
   if (evals[1] < evals[0]) {
      tval = evals[0] ;
      evals[0] = evals[1] ;
      evals[1] = tval ;
      }

   /* Make sure evals[0] < H[1][1] */
   if (evals[0] > H[1][1]) {
      printf("Warning: using H11 as eigenvalues\n") ;
      evals[0] = evals[1] = H[1][1] ;
      printf("Got evals[0] = H[1][1] = %12.7f\n", evals[0]) ;
      }

   /* get the eigenvectors */
   for (i=0; i<2; i++) {
      p = H[0][0] - evals[i] ;
      q = H[0][1] - S * evals[i] ;
      x = -p/q ;
      norm = 1.0 + x * x ;
      norm = sqrt(norm) ;
      evecs[i][0] = 1.0 / norm ;
      evecs[i][1] = x / norm ;
      }

   /* test
   for (i=0; i<2; i++) {
      p = H[i][0] * evecs[0][0] + H[i][1] * evecs[0][1] ;
      if (i==0) q = evecs[0][0] + S * evecs[0][1] ;
      else q = S * evecs[0][0] + evecs[0][1] ;
      q *= evals[0] ;
      printf("2x2 check %d: LHS = %12.6f  RHS = %12.6f\n", i, p, q) ;
      }
    */

}



/***
main()
{
double **H, **evecs, *evals ;
double S ;
void solve_2x2_pep() ;

   H = (double **) malloc (2 * sizeof(double *)) ;
   evecs = (double **) malloc (2 * sizeof(double *)) ;
   H[0] = (double *) malloc (2 * sizeof(double)) ;
   H[1] = (double *) malloc (2 * sizeof(double)) ;
   evecs[0] = (double *) malloc (2 * sizeof(double)) ;
   evecs[1] = (double *) malloc (2 * sizeof(double)) ;
   evals = (double *) malloc (2 * sizeof(double)) ;
   H[0][0] = -83.92663122885 ; H[0][1] = -83.92636151402 ;
   H[1][0] = -83.92636151402 ; H[1][1] = -83.92663613012 ;
   S = 0.999996726987 ;

   solve_2x2_pep(H, S, evals, evecs) ;
   printf("eval0 = %16.8f\n", evals[0]) ;
   printf("evec0 = (%12.6f %12.6f)\n", evecs[0][0], evecs[0][1]) ;
   printf("\neval1 = %16.8f\n", evals[1]) ;
   printf("evec1 = (%12.6f %12.6f)\n", evecs[1][0], evecs[1][1]) ;

}
***/

}
