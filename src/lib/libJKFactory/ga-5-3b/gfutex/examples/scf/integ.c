#include <math.h>
#include "cscc.h"
#include "ga.h"

double fm[2001][5];
double rdelta, delta, delo2;

double exprjh(double x)
{
//$Id: integ.F,v 1.1 2005-03-08 23:58:03 d3g293 Exp $
  double ret;
     
//     dumb solution to underflow problems on sun

  if (x < -37.0)  {
    ret = 0.0;
  } 
  else {
    ret = exp(x);
  }

  return ret;
}

void setfm(void)
{
  int i, ii;
  double t[2001];
  double et[2001], rr, tt;
  int MAXm;

  delta = 0.014;
  delo2 = delta * 0.5;
  rdelta = 1.0 / delta;
  //MAXm = 4;
  MAXm = 3; //Fortran Array Index 1==> C Array Index 0==>

  for (i = 0; i < 2001; i++) {
    tt = delta * (double)i; //i-1 ---> i
    et[i] = exprjh(-tt);
    t[i] = tt * 2.0;
    fm[i][MAXm + 1] = 0.0;
  }

  for (i = 199; i > MAXm; i--) {
    rr = 1.0 / (double)(2 * i + 1); //+1 ---> +1
    for (ii = 0; ii < 2001; ii++) {
      fm[ii][MAXm + 1] = (et[ii] + t[ii] * fm[ii][MAXm + 1]) * rr;
    }
  }

  for (i = MAXm; i >= 0; i--) {
    rr = 1.0 / (double) (2 * i + 1); //-1 ---> +1
    for (ii = 0; ii < 2001; ii++) {
      fm[ii][i] = (et[ii] + t[ii] * fm[ii][i+1]) * rr;
    }
  }

  return;
}

void f0(double *value, double t)
{
  const double fac0 = 0.88622692545276;
  const double rhalf = 0.5;
  const double rthird = (1.0 / 3.0);
  const double rquart = 0.25;

  double t0 = 28.0; //fortran: data --> C: long double or double?

  int n;
  double x;

  //     computes f0 to a relative accuracy of better than 4.e-13 for all t.
  //     uses 4th order taylor expansion on grid out to t=28.0
  //     asymptotic expansion accurate for t greater than 28

  if (t >= t0)
    *value = fac0 / sqrt(t);
  else {
    n = (int) ((t + delo2) * rdelta);
    x = delta * (double) n - t;
    //n=n+1; //c index 0, fortran index 1
    *value = fm[n][0] + x * (fm[n][1] + x * rhalf *
			     (fm[n][2] + x * rthird * (fm[n][3] + x * rquart * fm[n][4])));
  }

  return;
}

void addin(double *g, int *i, int *j, int *k, int *l, double *fock, double *dens, int *iky) //seems unused
{
    static double g2, g4, gg;
    static int ij, ik, il, jk, jl, kl;
    static double aij, aik, ajk, ail, gil;

    gg = *g;
    g2 = gg + gg;
    g4 = g2 + g2;
    ik = iky[*i] + *k;
    il = iky[*i] + *l;
    ij = iky[*i] + *j;
    jk = iky[MAX(*j,*k)] + MIN(*j,*k);
    jl = iky[MAX(*j,*l)] + MIN(*j,*l);
    kl = iky[*k] + *l;
    aij = g4 * dens[kl] + fock[ij];
    fock[kl] = g4 * dens[ij] + fock[kl];
    fock[ij] = aij;
    gil = gg;
    if (*i == *k || *j == *l) 
   {
	gg = g2;
    }
    if (*j == *k) 
    {
	gil = g2;
    }
    ajk = fock[jk] - gil * dens[il];
    ail = fock[il] - gil * dens[jk];
    aik = fock[ik] - gg * dens[jl];
    fock[jl] -= gg * dens[ik];
    fock[jk] = ajk;
    fock[il] = ail;
    fock[ik] = aik;

    return;
}

void dfill(int *n, double *val, double *a, int *ia) //seems unused
{
    int i;

    if (*ia == 1) 
    {
	for (i = 0; i < *n; i++) 
        {
	    a[i] = *val;
	}
    } 
    else 
    {
	for (i = 0; i < (*n - 1) * (*ia) + 1; i += *ia) 
        {
	    a[i] = *val;
	}
    }

    return;
}
