/*!
** \file
** \brief Converts a symmetric matrix to tridiagonal form for use in tqli
** \ingroup CIOMR
*/

#include <cmath>

#define DSIGN(a,b) ((b) >= 0.0) ? (fabs(a)) : (-fabs(a))

namespace psi {
  
/*!
** tred2(): converts symmetric matrix to a tridagonal form for use in tqli 
**
** if matz = 0, only find eigenvalues, else find both eigenvalues and
** eigenvectors 
**
** Returns: none
** \ingroup CIOMR
*/
void tred2(int n, double** a, double* d, double* e, int matz) {
    int i, j, k, l, il, ik, jk, kj;
    double f, g, h, hh, scale, scale_inv, h_inv;
    double temp;
    
    if (n == 1)
      return;
    
    for (i=n-1; i > 0; i--) {
      l = i-1;
      h = 0.0;
      scale = 0.0;
      if (l) {
        for (k=0; k <= l; k++) {
          scale += fabs(a[i][k]);
        }
        if (scale == 0.0) {
          e[i] = a[i][l];
        } else {
          scale_inv=1.0/scale;
          for (k=0; k <= l; k++) {
            a[i][k] *= scale_inv;
            h += a[i][k]*a[i][k];
          }
          f=a[i][l];
          g= -(DSIGN(sqrt(h),f));
          e[i] = scale*g;
          h -= f*g;
          a[i][l] = f-g;
          f = 0.0;
          h_inv=1.0/h;
          for (j=0; j <= l; j++) {
            if (matz)
              a[j][i] = a[i][j]*h_inv;
            g = 0.0;
            for (k=0; k <= j; k++) {
              g += a[j][k]*a[i][k];
            }
            if (l > j) {
              for (k=j+1; k <= l; k++) {
                g += a[k][j]*a[i][k];
              }
            }
            e[j] = g*h_inv;
            f += e[j]*a[i][j];
          }
          hh = f/(h+h);
          for (j=0; j <= l; j++) {
            f = a[i][j];
            g = e[j] - hh*f;
            e[j] = g;
            for (k=0; k <= j; k++) {
              a[j][k] -= (f*e[k]+ g*a[i][k]);
            }
          }
        }
      } else {
        e[i] = a[i][l];
      }
      d[i] = h;
    }
    if (matz)
      d[0] = 0.0;
    e[0] = 0.0;
    
    for (i=0; i < n; i++) {
      l = i-1;
      if (matz) {
        if (d[i]) {
          for (j=0; j <= l; j++) {
            g = 0.0;
            for (k=0; k <= l; k++) {
              g += a[i][k]*a[k][j];
            }
            for (k=0; k <= l; k++) {
              a[k][j] -= g*a[k][i];
            }
          }
        }
      }
      d[i] = a[i][i];
      if (matz) {
        a[i][i] = 1.0;
        if (l >= 0) {
          for (j=0; j<= l; j++) {
            a[i][j] = 0.0;
            a[j][i] = 0.0;
          }
        }
      }
    }
}

}

