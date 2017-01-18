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
** \file
** \brief Multiply two matrices (superceded by C_DGEMM)
** \ingroup CIOMR
*/

#include <cmath>
#include "psi4/libciomr/libciomr.h"

namespace psi {

static int keep_nr=0;
static int keep_nl=0;
static int keep_nc=0;
static double **aa,**bb;

/*!
** mmult():
** a reasonably fast matrix multiply (at least on the DEC3100)
** written by ETS
**
** \param AF = first matrix to multiply
** \param ta = if 1, transpose AF before multiplying; otherwise, 0
** \param BF = second matrix to multiply
** \param tb = if 1, transpose BF before multiplying; otherwise 0
** \param CF = matrix to hold result of AF*BF
** \param tc = if 1, transpose CF after the multiplication; otherwise 0
** \param nr = number of rows of AF
** \param nl = number of cols of AF and rows of BF
** \param nc = number of cols of BF
** \param add = if 1, add AF*BF to the matrix passed in as CF; else 0
**
** Returns: none
**
** nr,nl,nc are the number of rows,links,and columns in the
**          final matrices to be multiplied together
**          if ta=0 AF should have the dimensions nr x nl
**          if ta=1 AF should have the dimensions nl x nr
**          if tb=0 BF should have the dimensions nl x nc
**          if tb=1 BF should have the dimensions nc x nl
**          if tc=0 CF should have the dimensions nr x nc
**          if tc=1 CF should have the dimensions nc x nr
**
** \ingroup CIOMR
*/
void mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
  int nr, int nl, int nc, int add)
{
   int odd_nr,odd_nc,odd_nl;
   int i,j,k;
   double t00,t01,t10,t11;
   double **a,**b;
   double *att,*bt;
   double *at1,*bt1;

   if(!aa) {
      aa = (double **) init_matrix(nr,nl);
      bb = (double **) init_matrix(nc,nl);
      keep_nr = nr;
      keep_nl = nl;
      keep_nc = nc;
      }

   if(nl > keep_nl) {
      free_matrix(aa,keep_nr);
      free_matrix(bb,keep_nc);
      keep_nl = nl;
      keep_nr = (nr > keep_nr) ? nr : keep_nr;
      keep_nc = (nc > keep_nc) ? nc : keep_nc;
      aa = (double **) init_matrix(keep_nr,keep_nl);
      bb = (double **) init_matrix(keep_nc,keep_nl);
      }
   if(nr > keep_nr) {
      free_matrix(aa,keep_nr);
      keep_nr = nr;
      aa = (double **) init_matrix(keep_nr,keep_nl);
      }
   if(nc > keep_nc) {
      free_matrix(bb,keep_nc);
      keep_nc = nc;
      bb = (double **) init_matrix(keep_nc,keep_nl);
      }

   odd_nr = (nr)%2;
   odd_nc = (nc)%2;
   odd_nl = (nl)%2;

   a=aa;
   if(ta)
      for(i=0; i < nr ; i++)
         for(j=0; j < nl ; j++)
            a[i][j] = AF[j][i];
   else
      a=AF;

   b=bb;
   if(tb)
      b=BF;
   else
      for(i=0; i < nc ; i++)
         for(j=0; j < nl ; j++)
            b[i][j] = BF[j][i];

   for(j=0; j < nc-1 ; j+=2) {
      for(i=0; i < nr-1 ; i+=2) {
         att=a[i]; bt=b[j];
         at1=a[i+1]; bt1=b[j+1];
         if(add) {
            if(tc) {
               t00 = CF[j][i];
               t01 = CF[j+1][i];
               t10 = CF[j][i+1];
               t11 = CF[j+1][i+1];
               }
            else {
               t00 = CF[i][j];
               t01 = CF[i][j+1];
               t10 = CF[i+1][j];
               t11 = CF[i+1][j+1];
               }
            }
         else t00=t01=t10=t11=0.0;
         for(k=nl; k ; k--,att++,bt++,at1++,bt1++) {
            t00 += *att * *bt;
            t01 += *att * *bt1;
            t10 += *at1 * *bt;
            t11 += *at1 * *bt1;
            }
         if(tc) {
            CF[j][i]=t00;
            CF[j+1][i]=t01;
            CF[j][i+1]=t10;
            CF[j+1][i+1]=t11;
            }
         else {
            CF[i][j]=t00;
            CF[i][j+1]=t01;
            CF[i+1][j]=t10;
            CF[i+1][j+1]=t11;
            }
         }
      if(odd_nr) {
         att=a[i]; bt=b[j];
         bt1=b[j+1];
         if(add) {
            if(tc) {
               t00 = CF[j][i];
               t01 = CF[j+1][i];
               }
            else {
               t00 = CF[i][j];
               t01 = CF[i][j+1];
               }
            }
         else t00=t01=0.0;
         for(k= nl; k ; k--,att++,bt++,bt1++) {
            t00 += *att * *bt;
            t01 += *att * *bt1;
            }
         if(tc) {
            CF[j][i]=t00;
            CF[j+1][i]=t01;
            }
         else {
            CF[i][j]=t00;
            CF[i][j+1]=t01;
            }
         }
      }
   if(odd_nc) {
      for(i=0; i < nr-1 ; i+=2) {
         att=a[i]; bt=b[j];
         at1=a[i+1];
         if(add) {
            if(tc) {
               t00 = CF[j][i];
               t10 = CF[j][i+1];
               }
            else {
               t00 = CF[i][j];
               t10 = CF[i+1][j];
               }
            }
         else t00=t10=0.0;
         for(k= nl; k ; k--,att++,bt++,at1++) {
            t00 += *att * *bt;
            t10 += *at1 * *bt;
            }
         if(tc) {
            CF[j][i]=t00;
            CF[j][i+1]=t10;
            }
         else {
            CF[i][j]=t00;
            CF[i+1][j]=t10;
            }
         }
      if(odd_nr) {
         att=a[i]; bt=b[j];
         if(add)
            t00 = (tc) ? CF[j][i] : CF[i][j];
         else t00=0.0;
         for(k=nl; k ; k--,att++,bt++)
            t00 += *att * *bt;
         if(tc) CF[j][i]=t00;
         else CF[i][j]=t00;
         }
      }
   }

}
