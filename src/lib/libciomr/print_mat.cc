/*!
** \file
** \brief Print a matrix of doubles
** \ingroup CIOMR
*/

#include <cstdio>

namespace psi {

/*!
** print_mat: Print a matrix a of dimensions mxn to file pointer out.
**
** \param a   = matrix to print
** \param m   = number of rows in matrix
** \param n   = number of columns in matrix
** \param out = file pointer for output
**
** Returns: none
**
** \ingroup CIOMR
*/
void print_mat(double **a, int m, int n, FILE *out)
   {
      int ii,jj,kk,nn,ll;
      int i,j,k;

      ii=0;jj=0;
L200:
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
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
      }

}

