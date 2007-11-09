/*!
  \file eri.c
  By Edward Valeev
  \ingroup (QT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include "qt.h"

extern "C" {
	
#define MAXFAC 100
#define EPS 1.0E-17     /* Absolute precision in computing Fm(t)
                           (see recursion:calc_fij() ) */
#define MIN(a,b) ((a)>(b) ? (b) : (a))

static double *df;
static double *fac;
static double **bc;
static void calc_f(double *, int, double); 


/*!
  eri()

  This is a very inefficient function for computing ERIs
  over primitive Gaussian functions. The argument 
  list is self-explanatory, except for norm_flag:

  \param norm_flag:  tells what kind of normalization to use,
         0 - no normalization, >0 - normalized ERI
  \ingroup (QT)
*/

double eri(unsigned int l1, unsigned int m1, unsigned int n1, 
           double alpha1, double A[3],
	   unsigned int l2, unsigned int m2, unsigned int n2, 
           double alpha2, double B[3],
	   unsigned int l3, unsigned int m3, unsigned int n3, 
           double alpha3, double C[3],
	   unsigned int l4, unsigned int m4, unsigned int n4, 
           double alpha4, double D[3], int norm_flag)
{

  const double gammap = alpha1 + alpha2;
  const double Px = (alpha1*A[0] + alpha2*B[0])/gammap;
  const double Py = (alpha1*A[1] + alpha2*B[1])/gammap;
  const double Pz = (alpha1*A[2] + alpha2*B[2])/gammap;
  const double PAx = Px - A[0];
  const double PAy = Py - A[1];
  const double PAz = Pz - A[2];
  const double PBx = Px - B[0];
  const double PBy = Py - B[1];
  const double PBz = Pz - B[2];
  const double AB2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]); 
  
  const double gammaq = alpha3 + alpha4;
  const double gammapq = gammap*gammaq/(gammap+gammaq);
  const double Qx = (alpha3*C[0] + alpha4*D[0])/gammaq;
  const double Qy = (alpha3*C[1] + alpha4*D[1])/gammaq;
  const double Qz = (alpha3*C[2] + alpha4*D[2])/gammaq;
  const double QCx = Qx - C[0];
  const double QCy = Qy - C[1];
  const double QCz = Qz - C[2];
  const double QDx = Qx - D[0];
  const double QDy = Qy - D[1];
  const double QDz = Qz - D[2];
  const double CD2 = (C[0]-D[0])*(C[0]-D[0]) + (C[1]-D[1])*(C[1]-D[1]) + (C[2]-D[2])*(C[2]-D[2]);

  const double PQx = Px - Qx;
  const double PQy = Py - Qy;
  const double PQz = Pz - Qz;
  const double PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;

  int u1,u2,v1,v2,w1,w2,tx,ty,tz,txmax,tymax,tzmax;
  int i,j,k;
  int lp,lq,mp,mq,np,nq;
  int zeta;
  double *flp, *flq, *fmp, *fmq, *fnp, *fnq;
  double *F;
  double K1, K2;
  double Gx,Gy,Gz;
  double pfac;
  double result = 0.0;
  double tmp;
  int u1max,u2max,v1max,v2max,w1max,w2max;

  K1 = exp(-alpha1*alpha2*AB2/gammap);
  K2 = exp(-alpha3*alpha4*CD2/gammaq);
  pfac = 2*pow(M_PI,2.5)*K1*K2/(gammap*gammaq*sqrt(gammap+gammaq));


  if (fac == NULL) {
    fac = init_array(MAXFAC);
    fac[0] = 1.0;
    for(i=1;i<MAXFAC;i++)
      fac[i] = fac[i-1]*i;
    bc = block_matrix(MAXFAC,MAXFAC);
    for(i=0;i<MAXFAC;i++)
      for(j=0;j<=i;j++)
	bc[i][j] = fac[i]/(fac[i-j]*fac[j]);
  }

  if (norm_flag > 0) {
    pfac *= norm_const(l1,m1,n1,alpha1,A);
    pfac *= norm_const(l2,m2,n2,alpha2,B);
    pfac *= norm_const(l3,m3,n3,alpha3,C);
    pfac *= norm_const(l4,m4,n4,alpha4,D);
  }
  

  F = init_array(l1+l2+l3+l4+m1+m2+m3+m4+n1+n2+n3+n4+1);
  calc_f(F,l1+l2+l3+l4+m1+m2+m3+m4+n1+n2+n3+n4,PQ2*gammapq);

  
  flp = init_array(l1+l2+1);
  for(k=0;k<=l1+l2;k++)
    for(i=0;i<=MIN(k,l1);i++) {
      j = k-i;
      if (j > l2) continue;
      tmp = bc[l1][i]*bc[l2][j];
      if (l1-i > 0)
	tmp *= pow(PAx,l1-i);
      if (l2-j > 0)
	tmp *= pow(PBx,l2-j);
      flp[k] += tmp;
    }
  fmp = init_array(m1+m2+1);
  for(k=0;k<=m1+m2;k++)
    for(i=0;i<=MIN(k,m1);i++) {
      j = k-i;
      if (j > m2) continue;
      tmp = bc[m1][i]*bc[m2][j];
      if (m1-i > 0)
	tmp *= pow(PAy,m1-i);
      if (m2-j > 0)
	tmp *= pow(PBy,m2-j);
      fmp[k] += tmp;
    }
  fnp = init_array(n1+n2+1);
  for(k=0;k<=n1+n2;k++)
    for(i=0;i<=MIN(k,n1);i++) {
      j = k-i;
      if (j > n2) continue;
      tmp = bc[n1][i]*bc[n2][j];
      if (n1-i > 0)
	tmp *= pow(PAz,n1-i);
      if (n2-j > 0)
	tmp *= pow(PBz,n2-j);
      fnp[k] += tmp;
    }
  flq = init_array(l3+l4+1);
  for(k=0;k<=l3+l4;k++)
    for(i=0;i<=MIN(k,l3);i++) {
      j = k-i;
      if (j > l4) continue;
      tmp = bc[l3][i]*bc[l4][j];
      if (l3-i > 0)
	tmp *= pow(QCx,l3-i);
      if (l4-j > 0)
	tmp *= pow(QDx,l4-j);
      flq[k] += tmp;
    }
  fmq = init_array(m3+m4+1);
  for(k=0;k<=m3+m4;k++)
    for(i=0;i<=MIN(k,m3);i++) {
      j = k-i;
      if (j > m4) continue;
      tmp = bc[m3][i]*bc[m4][j];
      if (m3-i > 0)
	tmp *= pow(QCy,m3-i);
      if (m4-j > 0)
	tmp *= pow(QDy,m4-j);
      fmq[k] += tmp;
    }
  fnq = init_array(n3+n4+1);
  for(k=0;k<=n3+n4;k++)
    for(i=0;i<=MIN(k,n3);i++) {
      j = k-i;
      if (j > n4) continue;
      tmp = bc[n3][i]*bc[n4][j];
      if (n3-i > 0)
	tmp *= pow(QCz,n3-i);
      if (n4-j > 0)
	tmp *= pow(QDz,n4-j);
      fnq[k] += tmp;
    }


  for(lp=0;lp<=l1+l2;lp++)
    for(lq=0;lq<=l3+l4;lq++) {
      u1max = lp/2;
      u2max = lq/2;
      for(u1=0;u1<=u1max;u1++)
	for(u2=0;u2<=u2max;u2++) {
	  Gx = pow(-1,lp)*flp[lp]*flq[lq]*fac[lp]*fac[lq]*pow(gammap,u1-lp)*pow(gammaq,u2-lq)*fac[lp+lq-2*u1-2*u2]*
	       pow(gammapq,lp+lq-2*u1-2*u2)/(fac[u1]*fac[u2]*fac[lp-2*u1]*fac[lq-2*u2]);
	  for(mp=0;mp<=m1+m2;mp++)
	    for(mq=0;mq<=m3+m4;mq++) {
	      v1max = mp/2;
	      v2max = mq/2;
	      for(v1=0;v1<=v1max;v1++)
		for(v2=0;v2<=v2max;v2++) {
		  Gy = pow(-1,mp)*fmp[mp]*fmq[mq]*fac[mp]*fac[mq]*pow(gammap,v1-mp)*pow(gammaq,v2-mq)*fac[mp+mq-2*v1-2*v2]*
		       pow(gammapq,mp+mq-2*v1-2*v2)/(fac[v1]*fac[v2]*fac[mp-2*v1]*fac[mq-2*v2]);
		  for(np=0;np<=n1+n2;np++)
		    for(nq=0;nq<=n3+n4;nq++) {
		      w1max = np/2;
		      w2max = nq/2;
		      for(w1=0;w1<=w1max;w1++)
			for(w2=0;w2<=w2max;w2++) {
			  Gz = pow(-1,np)*fnp[np]*fnq[nq]*fac[np]*fac[nq]*pow(gammap,w1-np)*pow(gammaq,w2-nq)*fac[np+nq-2*w1-2*w2]*
			       pow(gammapq,np+nq-2*w1-2*w2)/(fac[w1]*fac[w2]*fac[np-2*w1]*fac[nq-2*w2]);
			  txmax = (lp+lq-2*u1-2*u2)/2;
			  tymax = (mp+mq-2*v1-2*v2)/2;
			  tzmax = (np+nq-2*w1-2*w2)/2;
			  for(tx=0;tx<=txmax;tx++)
			  for(ty=0;ty<=tymax;ty++)
			  for(tz=0;tz<=tzmax;tz++) {
			    zeta = lp+lq+mp+mq+np+nq-2*u1-2*u2-2*v1-2*v2-2*w1-2*w2-tx-ty-tz;
			    result += Gx*Gy*Gz*F[zeta]*
				      pow(-1,tx+ty+tz)*
				      pow(PQx,lp+lq-2*u1-2*u2-2*tx)*
				      pow(PQy,mp+mq-2*v1-2*v2-2*ty)*
				      pow(PQz,np+nq-2*w1-2*w2-2*tz)/
				      (pow(4,u1+u2+tx+v1+v2+ty+w1+w2+tz)*
				      pow(gammapq,tx)*
				      pow(gammapq,ty)*
				      pow(gammapq,tz)*
				      fac[lp+lq-2*u1-2*u2-2*tx]*fac[tx]*
				      fac[mp+mq-2*v1-2*v2-2*ty]*fac[ty]*
				      fac[np+nq-2*w1-2*w2-2*tz]*fac[tz]);
			  }
			}
		    }
		}
	    }
	}
    }

  free(F);
  free(flp);
  free(fmp);
  free(fnp);
  free(flq);
  free(fmq);
  free(fnq);

  return result*pfac;
}


/*!
  calc_f()

  This function computes infamous integral Fn(t). For its definition
  see Obara and Saika paper, or Shavitt's chapter in the
  Methods in Computational Physics book (see reference below).
  This piece of code is from Dr. Justin Fermann's program CINTS 

 \ingroup (QT)
*/
void calc_f(double *F, int n, double t)
{
  int i, m, k;
  int m2;
  double t2;
  double num;
  double sum;
  double term1, term2;
  static double K = 1.0/M_2_SQRTPI;
  double et;


  if (t>20.0){   /* For big t's do upward recursion */
    t2 = 2*t;
    et = exp(-t);
    t = sqrt(t);
    F[0] = K*erf(t)/t;
    for(m=0; m<=n-1; m++){
      F[m+1] = ((2*m + 1)*F[m] - et)/(t2);
    }
  }
  else {        /* For smaller t's compute F with highest n using
                   asymptotic series (see I. Shavitt in
                   Methods in Computational Physics, ed. B. Alder eta l,
                   vol 2, 1963, page 8) */
    et = exp(-t);
    t2 = 2*t;
    m2 = 2*n;
    num = df[m2];
    i=0;
    sum = 1.0/(m2+1);
    do{
      i++;
      num = num*t2;
      term1 = num/df[m2+2*i+2];
      sum += term1;
    } while (fabs(term1) > EPS && i < MAXFAC);
    F[n] = sum*et;
    for(m=n-1;m>=0;m--){        /* And then do downward recursion */
      F[m] = (t2*F[m+1] + et)/(2*m+1);
    }
  }
}


/*!
  norm_const()

  \ingroup (QT)
*/  
double norm_const(unsigned int l1, unsigned int m1, unsigned int n1, 
                  double alpha1, double A[3])
{
  int i;
  
  if (df == NULL) {
    df = init_array(2*MAXFAC);
    df[0] = 1.0;
    df[1] = 1.0;
    df[2] = 1.0;
    for(i=3; i<MAXFAC*2; i++) {
      df[i] = (i-1)*df[i-2];
    }
  }
  
  return pow(2*alpha1/M_PI,0.75)*pow(4*alpha1,0.5*(l1+m1+n1))/sqrt(df[2*l1]*df[2*m1]*df[2*n1]);
}

} /* extern "C" */