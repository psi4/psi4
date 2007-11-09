/*!
  \file eigout.cc
  \ingroup (CIOMR)
*/

#include "includes.h"

extern "C" {

/*!
** eigout(): Print out eigenvectors and eigenvalues.  Prints an n x m
** matrix of eigenvectors.  Under each eigenvector, the corresponding
** elements of two arrays, b and c, will also be printed.  This is
** useful for printing, for example, the SCF eigenvectors with their
** associated eigenvalues (orbital energies) and also the population.
**
** \ingroup (CIOMR)
*/
void eigout(double **a, double *b, double *c, int m, int n, FILE *out)
   {
      int ii,jj,kk,nn;
      int i,j;

      ii=0;jj=0;
L200:
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
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",c[j]);
         }
      fprintf (out,"\n");
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
      }

} /* extern "C" */
