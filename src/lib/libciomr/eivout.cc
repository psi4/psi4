/*!
  \file eivout.cc
  \ingroup (CIOMR)
*/
 
#include "includes.h"

extern "C" {

/*!
** eivout: Print out eigenvectors and eigenvalues to the output file
**
** \param a = eigenvectors
** \param b = eigenvalues
** \param m = rows of a
** \param n = columns of a
** \param out = output file pointer
**
** \ingroup (CIOMR)
*/
void eivout(double **a, double *b, int m, int n, FILE *out)
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
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
      }

} /* extern "C" */
