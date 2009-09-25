/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** test_bfgs.c
**
** Try the Numerical Recipies code 
**
** C. David Sherrill
** Georgia Institute of Technology
** March 2004
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>

namespace psi { namespace detcas {

void test_bfgs(void)
{
  int ndim, i, j;
  int iter, maxiter = 20;
  double E, E_last;
  double *x_cur, *x_last, *g_cur, *g_last;
  double *dx, *dg, *hdg;
  double **hessin;
  double fac, fad, fae;
  double dfunc(double *x, double *g);

  ndim = 2;
  E_last = 0.0;

  x_cur = init_array(ndim);
  x_last = init_array(ndim);
  g_cur = init_array(ndim);
  g_last = init_array(ndim);
  dx = init_array(ndim);
  dg = init_array(ndim);
  hdg = init_array(ndim);
  hessin = block_matrix(ndim,ndim);

  /* guess */
  x_cur[0] = 1.0;
  x_cur[1] = 5.0;

  for (i=0; i<ndim; i++) {
    hessin[i][i] = 1.0;
  }

  /* get value and gradient */
  E = dfunc(x_cur,g_cur);

  for (iter=0; iter<maxiter; iter++) {

    printf("Iteration %d.  Energy = %12.6lf\n", iter, E);
    printf("x and gradient\n");
    for (i=0; i<ndim; i++)
      printf("%12.6lf %12.6lf\n", x_cur[i], g_cur[i]);
    printf("\n");

    if (fabs(E - E_last) < 0.00001) {
      printf("Converged\n");
      exit(1);
    }

    for (i=0; i<ndim; i++) {
      dg[i] = g_cur[i] - g_last[i];
      dx[i] = x_cur[i] - x_last[i];
    }

    for (i=0; i<ndim; i++) {
      hdg[i]=0.0;
      for (j=0; j<ndim; j++) hdg[i] += hessin[i][j]*dg[j];
    }

    fac = fae = 0.0;
    for (i=0; i<ndim; i++) {
      fac += dg[i]*dx[i];
      fae += dg[i]*hdg[i];
    }
    fac = 1.0/fac;
    fad = 1.0/fae;
    for (i=0; i<ndim; i++) dg[i] = fac*dx[i] - fad*hdg[i];
    for (i=0; i<ndim; i++) {
      for (j=0; j<ndim; j++) {
        hessin[i][j] += fac*dx[i]*dx[j] - fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j];
      }
    }

    printf("Hessian inverse:\n");
    for (i=0; i<ndim; i++) {
      for (j=0; j<ndim; j++) {
        printf("%12.6lf ", hessin[i][j]);
      }
      printf("\n");
    }
    printf("\n");
 
    for (i=0; i<ndim; i++) {
      x_last[i] = x_cur[i];
      g_last[i] = g_cur[i];
    }

    for (i=0; i<ndim; i++) {
      for (j=0; j<ndim; j++) {
        x_cur[i] -= hessin[i][j] * g_cur[j];
      }
    }

    E_last = E;
    E = dfunc(x_cur,g_cur);


  } /* end iteration */

  free(x_cur);  free(x_last);  free(g_cur);  free(g_last);
  free(dx);  free(dg);  free(hdg);
  free_block(hessin);
  exit(0); 
}

double dfunc(double *x, double *g)
{
  double E;
  /* right now try (x-1)^2 + 10y^2, min is 0 at x=1,y=0 */
  E = (x[0] - 1.0) * (x[0] - 1.0) + 10.0 * x[1] * x[1];
  g[0] = 2.0 * (x[0] - 1.0);
  g[1] = 20.0 * x[1];
  return(E);
}

}} // end namespace psi::detcas

