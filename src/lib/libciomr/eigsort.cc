/*!
  \file eigsort.cc
  \ingroup (CIOMR)
*/

#include "includes.h"

extern "C" {

/*!
** eigsort: Sort the eigenvalues in d and eigenvectors in v in ascending
** (n>0) or descending (n<0) order.  abs(n) is the number of eigenvalues. 
** \ingroup (CIOMR)
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

/*
** MOSORT: Minor modification of eigsort() above to also sort a series of
** irrep labels.
**
** TDC, 6/03
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

} /* extern "C" */
