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
  \brief Sort eigenvalues and eigenvectors in ascending or descending order
  \ingroup CIOMR
*/


#include <cstdlib>

namespace psi {

/*!
** eigsort(): Sort the eigenvalues in d and eigenvectors in v in ascending
** (n>0) or descending (n<0) order.  abs(n) is the number of eigenvalues. 
** 
** \param d   = array of eigenvalues
** \param v   = matrix of eigenvectors (each column is an eigenvector)
**              Note: seems to assume v is a square matrix, could be a
**              problem if, e.g., nmo != nso.
** \param n   = abs(n) is the number of eigenvalues/cols of v.  Use
**              n>0 to sort in ascending order, n<0 to sort in descending
**              order
**
** Returns: none
**
** \ingroup CIOMR
*/
void eigsort(double *d, double **v, int n)
{
  int i,j,k;
  double p;

  /* Modified by Ed Valeev - if n is negative,
     sort eigenvalues in descending order */

  if (n >= 0) {
    for (i=0; i < n-1 ; i++) {
      k=i;
      p=d[i];
      for (j=i+1; j < n; j++) {
	if (d[j] < p) {
	  k=j;
	  p=d[j];
	}
      }
      if (k != i) {
	d[k]=d[i];
	d[i]=p;
	for (j=0; j < n; j++) {
	  p=v[j][i];
	  v[j][i]=v[j][k];
	  v[j][k]=p;
	}
      }
    }
  }
  else {
    n = abs(n);
    for (i=0; i < n-1 ; i++) {
      k=i;
      p=d[i];
      for (j=i+1; j < n; j++) {
	if (d[j] > p) {
	  k=j;
	  p=d[j];
	}
      }
      if (k != i) {
	d[k]=d[i];
	d[i]=p;
	for (j=0; j < n; j++) {
	  p=v[j][i];
	  v[j][i]=v[j][k];
	  v[j][k]=p;
	}
      }
    }
  }
}


/*!
** mosort(): Minor modification of eigsort() to also sort a series of
** irrep labels.
**
** \param d   = array of eigenvalues
** \param v   = matrix of eigenvectors (each column is an eigenvector)
** \param sym = array of symmetry ID's (irreps)
** \param nso = number of rows in v
** \param nmo = abs(nmo) is the number of eigenvalues/cols of v.  Use
**              nmo>0 to sort in ascending order, nmo<0 to sort in descending
**              order
** 
** Returns:none
**
** TDC, 6/03
** \ingroup CIOMR
*/
void mosort(double *d, double **v, int *sym, int nso, int nmo)
{
  int i, j, k, l;
  double p;

  if(nmo > 0) {
    for (i=0; i < nmo-1 ; i++) {
      k=i;
      p=d[i];
      for (j=i+1; j < nmo; j++) {
	if (d[j] < p) {
	  k=j;
	  p=d[j];
	}
      }
      if (k != i) {
	d[k]=d[i];
	d[i]=p;

	l = sym[i];
	sym[i] = sym[k];
	sym[k] = l;

	for (j=0; j < nso; j++) {
	  p=v[j][i];
	  v[j][i]=v[j][k];
	  v[j][k]=p;
	}
      }
    }
  }
  else if(nmo < 0) {
    nmo = abs(nmo);
    for (i=0; i < nmo-1 ; i++) {
      k=i;
      p=d[i];
      for (j=i+1; j < nmo; j++) {
	if (d[j] > p) {
	  k=j;
	  p=d[j];
	}
      }
      if (k != i) {
	d[k]=d[i];
	d[i]=p;

	l = sym[i];
	sym[i] = sym[k];
	sym[k] = l;

	for (j=0; j < nso; j++) {
	  p=v[j][i];
	  v[j][i]=v[j][k];
	  v[j][k]=p;
	}
      }
    }
  }
}

}
