/*!
** \file
** \brief Print a lower-triangle array of doubles
** \ingroup CIOMR
*/

#include <cstdio>

namespace psi {

/*!
** print_array(): Prints a lower-triangle of a symmetric matrix packed as
**  an array of doubles.
**
** \param a   = array (packed lower triangle of matrix) to print
** \param m   = dimension of matrix (mxm)
** \param out = file pointer for output
**
** Returns: none
**
** \ingroup CIOMR
*/
void print_array(double *a, int m, FILE *out)
   {
      int ii,jj,kk,mm,nn,ll;
      int i,j,k,i1,i2;

      ii=0;jj=0;
L200:
      ii++;
      jj++;
      kk=10*jj;
      nn = kk + kk*(kk-1)/2;
      mm=m;
      if (m > kk) mm=kk;
      ll = 2*(mm-ii+1)+1;
      fprintf (out,"\n");
      for (i=ii; i <= mm; i++) fprintf(out,"       %5d",i);
      fprintf (out,"\n");
      for (i=ii; i <= m; i++) {
         i1=i*(i-1)/2+ii;
         i2=i+i*(i-1)/2;
         if (i2 > nn) i2 = i1+9;
         fprintf (out,"\n%5d",i);
         for (j=i1; j <= i2; j++) {
            fprintf (out,"%12.7f",a[j-1]);
            }
         }
      if (m <= kk) {
         fprintf(out,"\n");
         fflush(out);
         return;
         }
      ii=kk; goto L200;
      }

}

