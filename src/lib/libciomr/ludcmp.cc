
#include <cstdlib>
#include <cmath>
#include "libciomr.h"

namespace psi {

void ludcmp(double** a,int n,int* indx,double* d)
   {
      int i,imax,j,k;
      double big,dum,sum,temp;
      double *vv;

      vv = (double *) init_array(n);

      *d = 1.0;

      for (i=0; i < n ; i++) {
         big=0.0;
         for (j=0; j < n; j++) {
            if ((temp=fabs(a[i][j])) > big) big=temp;
            }
         if (big == 0.0) {
            *d = 0.0;
            return;
            }
         vv[i] = 1.0/big;
         }
      for (j=0; j < n ; j++) {
         for (i=0; i < j ; i++) {
            sum = a[i][j];
            for (k=0; k < i ; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            }
         big = 0.0;
         for (i=j ; i < n ; i++) {
            sum=a[i][j];
            for (k=0; k < j ; k++) sum -= a[i][k]*a[k][j];
            a[i][j] = sum;
            if ((dum=vv[i]*fabs(sum)) >= big) {
               big = dum;
               imax = i;
               }
            }
         if (j != imax) {
            for (k=0; k < n; k++) {
               dum=a[imax][k];
               a[imax][k]=a[j][k];
               a[j][k]=dum;
               }
            *d = -(*d);
            vv[imax]=vv[j];
            }
         indx[j]=imax;
         if (a[j][j] == 0.0) a[j][j] = 1.0e-20;
         if (j != n-1) {
            dum = 1.0/a[j][j];
            for (i=j+1; i < n ; i++) a[i][j] *= dum;
            }
         }
      free(vv);
      }

}

