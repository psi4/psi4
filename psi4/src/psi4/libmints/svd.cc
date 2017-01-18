/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "psi4/psi4-dec.h"
typedef double real;
typedef int integer;

#define r_sign(a,b) ((*b<0.0)?-fabs(*a):fabs(*a))
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define dmax(a,b) (((a)>(b))?(a):(b))

int
sing_(double *q, int *lq, int *iq, double *s, double *p,
     int *lp, int *ip, double *a, int *la, int *m, int *n, double *w);

static int
bidag2_(double *d, double *b, double *q, int *lq, int *iq,
        double *p, int *lp, int *ip, double *a, int *la, int *m, int *n);
static int
hsr1_(double *a, int *la, int *n);
static int
hsr2_(double *a, int *la, int *n);
static int
hsr3_(double *a, int *la, int *m, int *n);
static int
hsr4_(double *a, int *la, int *m, int *n);
static int
hsr5_(double *a, int *la, int *m, int *n);
static int
singb_(double *d, int *n, double *u, int *iu,
       double *q, int *lq, int *mq, int *iq,
       double *p, int *lp, int *mp, int *ip, double *e, double *f);
static int
sng0_(double *q, int *lq, int *m, double *p, int *lp, int *n,
      int *l, int *j, int *k, double *x, double *y);
static int
sft_(double *s, double *a, double *b, double *c,
     double *d, double *e2, double *e1, double *e0, double *f2, double *f1);
static int
sng1_(double *q, int *lq, int *m, double *p,
      int *lp, int *n, int *l, int *j, int *k, double *x, double *y);
static int
scl_(double *d, double *u, int *n, double *q,
     int *lq, int *mq, double *p, int *lp, int *mp,
     double *e, double *f, double *b,
     int *j, int *k, int *jl,
     int *jr);
static int
eig3_(double *ea, double *eb, double *a, double *b, double *y, double *z);
static int
sort2_(double *x, double *y, double *w, int *n);
static int
fgv_(double *x, double *y, double *s, double *p, double *q,
     double *a, double *b);

#ifdef TEST

main(int argc, char**argv)
{
  int i,j;
  double *tmp;

  int m = atoi(argv[1]);
  int n = atoi(argv[2]);

  int l = ((m<n)?m:n);

  double *q = (double*)malloc(sizeof(double)*m*m);
  double *s = (double*)malloc(sizeof(double)*n);
  double *p = (double*)malloc(sizeof(double)*n*n);
  double *a = (double*)malloc(sizeof(double)*n*m);
  double *w = (double*)malloc(sizeof(double)*3*m);

  int lq = m;
  int lp = n;
  int iq = 3;
  int ip = 3;
  int la = lq;

  tmp = a;
  for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {
          *tmp++ = drand48() * ((drand48()<0.5)?1.0:-1.0);
          /* *tmp++ = 1.0; */
        }
    }

  printf("A:\n");
  for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {
          printf(" % 6.4f", a[j*m + i]);
        }
      printf("\n");
    }

  sing_(q, &lq, &iq, s, p,
        &lp, &ip, a, &la, &m, &n, w);

  printf("Q:\n");
  for (i=0; i<m; i++) {
      for (j=0; j<m; j++) {
          printf(" % 6.4f",q[j*m+i]);
        }
      printf("\n");
    }

  printf("P:\n");
  for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
          printf(" % 6.4f",p[j*n+i]);
        }
      printf("\n");
    }

  printf("S:\n");
  tmp = s;
  for (j=0; j<l; j++) {
      printf(" % 6.4f",*tmp++);
    }
  printf("\n");

  printf("QQt:\n");
  for (i=0; i<m; i++) {
      for (j=0; j<m; j++) {
          int k;
          double tmp = 0.0;
          for (k=0; k<m; k++) {
              tmp += q[k*m+i]*q[k*m+j];
            }
          printf(" % 6.4f",tmp);
        }
      printf("\n");
    }

  printf("QtQ:\n");
  for (i=0; i<m; i++) {
      for (j=0; j<m; j++) {
          int k;
          double tmp = 0.0;
          for (k=0; k<m; k++) {
              tmp += q[i*m+k]*q[j*m+k];
            }
          printf(" % 6.4f",tmp);
        }
      printf("\n");
    }

  printf("PPt:\n");
  for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
          int k;
          double tmp = 0.0;
          for (k=0; k<l; k++) {
              tmp += p[k*n+i]*p[k*n+j];
            }
          printf(" % 6.4f",tmp);
        }
      printf("\n");
    }

  printf("PtP:\n");
  for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
          int k;
          double tmp = 0.0;
          for (k=0; k<l; k++) {
              tmp += p[i*n+k]*p[j*n+k];
            }
          printf(" % 6.4f",tmp);
        }
      printf("\n");
    }

  printf("QSPt:\n");
  for (i=0; i<m; i++) {
      for (j=0; j<n; j++) {
          int k;
          double tmp = 0.0;
          for (k=0; k<l; k++) {
              tmp += q[k*m+i]*s[k]*p[k*n+j];
            }
          printf(" % 6.4f",tmp);
        }
      printf("\n");
    }


  return 0;
}

#endif /* TEST */

/*      ________________________________________________________ */
/*     |                                                        | */
/*     |COMPUTE SINGULAR VALUE DECOMPOSITION OF A GENERAL MATRIX| */
/*     |      A = Q TIMES DIAGONAL MATRIX TIMES P TRANSPOSE     | */
/*     |                                                        | */
/*     |    INPUT:                                              | */
/*     |                                                        | */
/*     |         S     --ARRAY WITH AT LEAST N ELEMENTS         | */
/*     |                                                        | */
/*     |         LQ    --LEADING (ROW) DIMENSION OF ARRAY Q     | */
/*     |                                                        | */
/*     |         IQ    --AN INTEGER WHICH INDICATES WHICH COL-  | */
/*     |                 UMNS OF Q TO COMPUTE (= 0 MEANS NONE,  | */
/*     |                 = 1 MEANS FIRST L, = 2 MEANS LAST M-L, | */
/*     |                 = 3 MEANS ALL M WHERE L = MIN(M,N))    | */
/*     |                                                        | */
/*     |         LP    --LEADING (ROW) DIMENSION OF ARRAY P     | */
/*     |                                                        | */
/*     |         IP    --AN INTEGER (LIKE IQ) WHICH INDICATES   | */
/*     |                 WHICH COLUMNS OF P TO COMPUTE          | */
/*     |                                                        | */
/*     |         A     --ARRAY CONTAINING COEFFICIENT MATRIX    | */
/*     |                                                        | */
/*     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     | */
/*     |                                                        | */
/*     |         M     --ROW DIMENSION OF MATRIX STORED IN A    | */
/*     |                                                        | */
/*     |         N     --COLUMN DIMENSION OF MATRIX STORED IN A | */
/*     |                                                        | */
/*     |         W     --WORK ARRAY(LENGTH AT LEAST MAX(M,3L-1))| */
/*     |                                                        | */
/*     |    OUTPUT:                                             | */
/*     |                                                        | */
/*     |         Q     --Q FACTOR IN THE SINGULAR VALUE DECOMP. | */
/*     |                                                        | */
/*     |         S     --SINGULAR VALUES IN DECREASING ORDER    | */
/*     |                                                        | */
/*     |         P     --P FACTOR IN THE SINGULAR VALUE DECOMP. | */
/*     |                                                        | */
/*     |         A     --THE HOUSEHOLDER VECTORS USED IN THE    | */
/*     |                 REDUCTION PROCESS                      | */
/*     |                                                        | */
/*     |    NOTE:                                               | */
/*     |                                                        | */
/*     |         EITHER P OR Q CAN BE IDENTIFIED WITH A BUT NOT | */
/*     |         BOTH. WHEN EITHER P OR Q ARE IDENTIFIED WITH A,| */
/*     |         THE HOUSEHOLDER VECTORS IN A ARE DESTROYED. IF | */
/*     |         IQ = 2, Q MUST HAVE M COLUMNS EVEN THOUGH THE  | */
/*     |         OUTPUT ARRAY HAS JUST M-L COLUMNS. SIMILARLY IF| */
/*     |         IP = 2, P MUST HAVE N COLUMNS EVEN THOUGH THE  | */
/*     |         OUTPUT ARRAY HAS JUST N-L COLUMNS.             | */
/*     |                                                        | */
/*     |    BUILTIN FUNCTIONS: MIN0                             | */
/*     |    PACKAGE SUBROUTINES: BIDAG2,EIG3,FGV,HSR1-HSR5,SCL, | */
/*     |                         SFT,SINGB,SNG0,SNG1,SORT2      | */
/*     |________________________________________________________| */

int
sing_(real *q, int *lq, int *iq, double *s, double *p,
      int *lp, int *ip, double *a, int *la, int *m, int *n, double *w)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, p_dim1, p_offset, i__1;

    /* Local variables */
    static integer i, l;
    static integer jl, iu, jr;

    /* Parameter adjustments */
    --w;
    a_dim1 = *la;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    p_dim1 = *lp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    --s;
    q_dim1 = *lq;
    q_offset = q_dim1 + 1;
    q -= q_offset;

    /* Function Body */
    if (*iq >= 0) {
	goto L20;
    }
L10:
    psi::outfile->Printf("ERROR: INPUT PARAMETER IQ FOR SUBROUTINE SING\n");
    psi::outfile->Printf("EITHER LESS THAN 0 OR GREATER THAN 3\n");
    abort();
L20:
    if (*iq > 3) {
	goto L10;
    }
    jl = 0;
    if (*iq == 0) {
	goto L30;
    }
    if (*iq == 2) {
	goto L30;
    }
    jl = 1;
L30:
    if (*ip >= 0) {
	goto L50;
    }
L40:
    psi::outfile->Printf("ERROR: INPUT PARAMETER IP FOR SUBROUTINE SING\n");
    psi::outfile->Printf("EITHER LESS THAN 0 OR GREATER THAN 3\n");
    abort();
L50:
    if (*ip > 3) {
	goto L40;
    }
    jr = 0;
    if (*ip == 0) {
	goto L60;
    }
    if (*ip == 2) {
	goto L60;
    }
    jr = 1;
L60:
    bidag2_(&s[1], &w[1], &q[q_offset], lq, iq, &p[p_offset], lp, ip, &a[
	    a_offset], la, m, n);
    l = min(*m,*n);
    if (l > 1) {
	goto L80;
    }
    if (s[1] >= (double)0.) {
	return 0;
    }
    s[1] = -(double)s[1];
    if ((real) jl == (double)0.) {
	return 0;
    }
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
/* L70: */
	q[i + q_dim1] = -(double)q[i + q_dim1];
    }
    return 0;
L80:
    iu = 0;
    if (*m >= *n) {
	goto L90;
    }
    iu = 1;
L90:
    singb_(&s[1], &l, &w[1], &iu, &q[q_offset], lq, m, &jl, &p[p_offset], lp,
	    n, &jr, &w[l], &w[l + l]);
    return 0;
} /* sing_ */


/*      ________________________________________________________ */
/*     |                                                        | */
/*     |      REDUCE A GENERAL MATRIX A TO BIDIAGONAL FORM      | */
/*     |     A = Q TIMES BIDIAGONAL MATRIX TIMES P TRANSPOSE    | */
/*     |                                                        | */
/*     |    INPUT:                                              | */
/*     |                                                        | */
/*     |         D     --ARRAY WITH AT LEAST N ELEMENTS         | */
/*     |                                                        | */
/*     |         B     --ARRAY WITH AT LEAST M ELEMENTS         | */
/*     |                                                        | */
/*     |         LQ    --LEADING (ROW) DIMENSION OF ARRAY Q     | */
/*     |                                                        | */
/*     |         IQ    --AN INTEGER WHICH INDICATES WHICH COL-  | */
/*     |                 UMNS OF Q TO COMPUTE (= 0 MEANS NONE,  | */
/*     |                 = 1 MEANS FIRST L, = 2 MEANS LAST M-L, | */
/*     |                 = 3 MEANS ALL M WHERE L = MIN(M,N))    | */
/*     |                                                        | */
/*     |         LP    --LEADING (ROW) DIMENSION OF ARRAY P     | */
/*     |                                                        | */
/*     |         IP    --AN INTEGER (LIKE IQ) WHICH INDICATES   | */
/*     |                 WHICH COLUMNS OF P TO COMPUTE          | */
/*     |                                                        | */
/*     |         A     --ARRAY CONTAINING COEFFICIENT MATRIX    | */
/*     |                                                        | */
/*     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     | */
/*     |                                                        | */
/*     |         M     --ROW DIMENSION OF MATRIX STORED IN A    | */
/*     |                                                        | */
/*     |         N     --COLUMN DIMENSION OF MATRIX STORED IN A | */
/*     |                                                        | */
/*     |    OUTPUT:                                             | */
/*     |                                                        | */
/*     |         D     --DIAGONAL OF BIDIAGONAL FORM            | */
/*     |                                                        | */
/*     |         B     --SUPERDIAGONAL (IF M .GE. N) OR SUBDIAG-| */
/*     |                 ONAL (IF M .LT. N) OF BIDIAGONAL FORM  | */
/*     |                                                        | */
/*     |         A     --THE HOUSEHOLDER VECTORS USED IN THE    | */
/*     |                 REDUCTION PROCESS                      | */
/*     |                                                        | */
/*     |         Q     --THE Q FACTOR IN THE BIDIAGONALIZATION  | */
/*     |                                                        | */
/*     |         P     --THE P FACTOR IN THE BIDIAGONALIZATION  | */
/*     |                                                        | */
/*     |    NOTE:                                               | */
/*     |                                                        | */
/*     |         EITHER P OR Q CAN BE IDENTIFIED WITH A BUT NOT | */
/*     |         BOTH. WHEN EITHER P OR Q ARE IDENTIFIED WITH A,| */
/*     |         THEN THE HOUSEHOLDER VECTORS IN A ARE DESTROYED| */
/*     |                                                        | */
/*     |    BUILTIN FUNCTIONS: ABS,MIN0,SQRT                    | */
/*     |    PACKAGE SUBROUTINES: HSR1-HSR5                      | */
/*     |________________________________________________________| */

static int
bidag2_(double *d, double *b, double *q, int *lq, int *iq,
        double *p, int *lp, int *ip, double *a, int *la, int *m, int *n)
{
    /* System generated locals */
    integer a_dim1, a_offset, q_dim1, q_offset, p_dim1, p_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    static integer h, i, j, k, l;
    static real r, s, t, u;
    static integer jp, jq;

    /* Parameter adjustments */
    a_dim1 = *la;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    p_dim1 = *lp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    q_dim1 = *lq;
    q_offset = q_dim1 + 1;
    q -= q_offset;
    --b;
    --d;

    /* Function Body */
    l = min(*m,*n);
    if (*iq >= 0) {
	goto L20;
    }
L10:
    psi::outfile->Printf("ERROR: INPUT PARAMETER IQ FOR SUBROUTINE BIDAG2\n");
    psi::outfile->Printf("EITHER LESS THAN 0 OR GREATER THAN 3\n");
    abort();
L20:
    if (*iq > 3) {
	goto L10;
    }
    jq = *iq;
    if (*iq <= 1) {
	goto L30;
    }
    if (*iq == 3) {
	goto L30;
    }
    if (*m == l) {
	jq = 0;
    }
L30:
    if (*ip >= 0) {
	goto L50;
    }
L40:
    psi::outfile->Printf("ERROR: INPUT PARAMETER IP FOR SUBROUTINE BIDAG2\n");
    psi::outfile->Printf("EITHER LESS THAN 0 OR GREATER THAN 3\n");
    abort();
L50:
    if (*ip > 3) {
	goto L40;
    }
    jp = *ip;
    if (*ip <= 1) {
	goto L60;
    }
    if (*ip == 3) {
	goto L60;
    }
    if (*n == l) {
	jp = 0;
    }
L60:
    k = 1;
    h = 2;
    if (*m < *n) {
	goto L330;
    }
    if (*m > 1) {
	goto L70;
    }
    d[1] = a[a_dim1 + 1];
    if (*iq > 0) {
	q[q_dim1 + 1] = (double)1.;
    }
    if (*ip > 0) {
	p[p_dim1 + 1] = (double)1.;
    }
    return 0;
L70:
    j = k;
    k = h;
    i__1 = *m;
    for (i = k; i <= i__1; ++i) {
/* L80: */
	if (a[i + j * a_dim1] != (double)0.) {
	    goto L110;
	}
    }
    d[j] = a[j + j * a_dim1];
    a[j + j * a_dim1] = (double)0.;
    i__1 = *n;
    for (i = k; i <= i__1; ++i) {
/* L90: */
	d[i] = a[j + i * a_dim1];
    }
    if (jq == 0) {
	goto L200;
    }
    i__1 = *m;
    for (i = j; i <= i__1; ++i) {
/* L100: */
	q[i + j * q_dim1] = (double)0.;
    }
    goto L200;
L110:
    t = (r__1 = a[j + j * a_dim1], fabs(r__1));
    if (t != (double)0.) {
	u = (double)1. / t;
    }
    r = (double)1.;
    i__1 = *m;
    for (l = i; l <= i__1; ++l) {
	s = (r__1 = a[l + j * a_dim1], fabs(r__1));
	if (s <= t) {
	    goto L120;
	}
	u = (double)1. / s;
/* Computing 2nd power */
	r__1 = t * u;
	r = r * (r__1 * r__1) + (double)1.;
	t = s;
	goto L130;
L120:
/* Computing 2nd power */
	r__1 = s * u;
	r += r__1 * r__1;
L130:
	;
    }
    s = t * sqrt(r);
    r = a[j + j * a_dim1];
    u = (double)1. / sqrt(s * (s + fabs(r)));
    if (r < (double)0.) {
	s = -(double)s;
    }
    d[j] = -(double)s;
    a[j + j * a_dim1] = u * (r + s);
    i__1 = *m;
    for (i = k; i <= i__1; ++i) {
/* L140: */
	a[i + j * a_dim1] *= u;
    }
    if (jq == 0) {
	goto L160;
    }
    i__1 = *m;
    for (i = j; i <= i__1; ++i) {
/* L150: */
	q[i + j * q_dim1] = a[i + j * a_dim1];
    }
L160:
    if (k > *n) {
	goto L620;
    }
    i__1 = *n;
    for (l = k; l <= i__1; ++l) {
	t = (double)0.;
	i__2 = *m;
	for (i = j; i <= i__2; ++i) {
/* L170: */
	    t += a[i + j * a_dim1] * a[i + l * a_dim1];
	}
	a[j + l * a_dim1] -= t * a[j + j * a_dim1];
	d[l] = a[j + l * a_dim1];
	i__2 = *m;
	for (i = k; i <= i__2; ++i) {
/* L180: */
	    a[i + l * a_dim1] -= t * a[i + j * a_dim1];
	}
/* L190: */
    }
L200:
    h = k + 1;
    if (k < *n) {
	goto L210;
    }
    if (k > *n) {
	goto L620;
    }
    if (*m == *n) {
	goto L610;
    }
    b[j] = a[j + *n * a_dim1];
    goto L70;
L210:
    i__1 = *n;
    for (i = h; i <= i__1; ++i) {
/* L220: */
	if (d[i] != (double)0.) {
	    goto L240;
	}
    }
    b[j] = d[k];
    a[j + k * a_dim1] = (double)0.;
    if (*ip == 0) {
	goto L70;
    }
    i__1 = *n;
    for (i = k; i <= i__1; ++i) {
/* L230: */
	p[i + j * p_dim1] = (double)0.;
    }
    goto L70;
L240:
    t = (r__1 = d[k], fabs(r__1));
    if (t != (double)0.) {
	u = (double)1. / t;
    }
    r = (double)1.;
    i__1 = *n;
    for (l = i; l <= i__1; ++l) {
	s = (r__1 = d[l], fabs(r__1));
	if (s <= t) {
	    goto L250;
	}
	u = (double)1. / s;
/* Computing 2nd power */
	r__1 = t * u;
	r = r * (r__1 * r__1) + (double)1.;
	t = s;
	goto L260;
L250:
/* Computing 2nd power */
	r__1 = s * u;
	r += r__1 * r__1;
L260:
	;
    }
    s = t * sqrt(r);
    r = d[k];
    u = (double)1. / sqrt(s * (s + fabs(r)));
    if (r < (double)0.) {
	s = -(double)s;
    }
    d[k] = u * (r + s);
    i__1 = *n;
    for (i = h; i <= i__1; ++i) {
/* L270: */
	d[i] *= u;
    }
    if (*ip == 0) {
	goto L290;
    }
    i__1 = *n;
    for (i = k; i <= i__1; ++i) {
/* L280: */
	p[i + j * p_dim1] = d[i];
    }
L290:
    b[j] = -(double)s;
    i__1 = *m;
    for (i = k; i <= i__1; ++i) {
/* L300: */
	b[i] = (double)0.;
    }
    i__1 = *n;
    for (l = k; l <= i__1; ++l) {
	t = d[l];
	a[j + l * a_dim1] = t;
	i__2 = *m;
	for (i = k; i <= i__2; ++i) {
/* L310: */
	    b[i] += t * a[i + l * a_dim1];
	}
    }
    i__2 = *n;
    for (l = k; l <= i__2; ++l) {
	t = d[l];
	i__1 = *m;
	for (i = k; i <= i__1; ++i) {
/* L320: */
	    a[i + l * a_dim1] -= t * b[i];
	}
    }
    goto L70;
L330:
    i__1 = *n;
    for (i = k; i <= i__1; ++i) {
/* L340: */
	d[i] = a[k + i * a_dim1];
    }
L350:
    j = k;
    k = h;
    i__1 = *n;
    for (i = k; i <= i__1; ++i) {
/* L360: */
	if (d[i] != (double)0.) {
	    goto L370;
	}
    }
    u = d[j];
    d[j] = (double)0.;
    goto L440;
L370:
    t = (r__1 = d[j], fabs(r__1));
    if (t != (double)0.) {
	u = (double)1. / t;
    }
    r = (double)1.;
    i__1 = *n;
    for (l = i; l <= i__1; ++l) {
	s = (r__1 = d[l], fabs(r__1));
	if (s <= t) {
	    goto L380;
	}
	u = (double)1. / s;
/* Computing 2nd power */
	r__1 = t * u;
	r = r * (r__1 * r__1) + (double)1.;
	t = s;
	goto L390;
L380:
/* Computing 2nd power */
	r__1 = s * u;
	r += r__1 * r__1;
L390:
	;
    }
    s = t * sqrt(r);
    r = d[j];
    u = (double)1. / sqrt(s * (s + fabs(r)));
    if (r < (double)0.) {
	s = -(double)s;
    }
    d[j] = u * (r + s);
    i__1 = *n;
    for (i = k; i <= i__1; ++i) {
/* L400: */
	d[i] *= u;
    }
    u = -(double)s;
    if (k > *m) {
	goto L470;
    }
    i__1 = *m;
    for (i = k; i <= i__1; ++i) {
/* L410: */
	b[i] = (double)0.;
    }
    i__1 = *n;
    for (l = j; l <= i__1; ++l) {
	t = d[l];
	a[j + l * a_dim1] = t;
	i__2 = *m;
	for (i = k; i <= i__2; ++i) {
/* L420: */
	    b[i] += t * a[i + l * a_dim1];
	}
    }
    i__2 = *n;
    for (l = j; l <= i__2; ++l) {
	t = d[l];
	i__1 = *m;
	for (i = k; i <= i__1; ++i) {
/* L430: */
	    a[i + l * a_dim1] -= t * b[i];
	}
    }
L440:
    h = k + 1;
    if (*ip == 0) {
	goto L460;
    }
    i__1 = *n;
    for (i = j; i <= i__1; ++i) {
/* L450: */
	p[i + j * p_dim1] = d[i];
    }
L460:
    d[j] = u;
    if (k < *m) {
	goto L490;
    }
    if (k > *m) {
	goto L620;
    }
    b[j] = a[*m + j * a_dim1];
    goto L330;
L470:
    i__1 = *n;
    for (i = j; i <= i__1; ++i) {
/* L480: */
	a[j + i * a_dim1] = d[i];
    }
    goto L440;
L490:
    i__1 = *m;
    for (i = h; i <= i__1; ++i) {
/* L500: */
	if (a[i + j * a_dim1] != (double)0.) {
	    goto L520;
	}
    }
    b[j] = a[k + j * a_dim1];
    a[k + j * a_dim1] = (double)0.;
    if (*iq == 0) {
	goto L330;
    }
    i__1 = *m;
    for (i = k; i <= i__1; ++i) {
/* L510: */
	q[i + j * q_dim1] = (double)0.;
    }
    goto L330;
L520:
    t = (r__1 = a[k + j * a_dim1], fabs(r__1));
    if (t != (double)0.) {
	u = (double)1. / t;
    }
    r = (double)1.;
    i__1 = *m;
    for (l = i; l <= i__1; ++l) {
	s = (r__1 = a[l + j * a_dim1], fabs(r__1));
	if (s <= t) {
	    goto L530;
	}
	u = (double)1. / s;
/* Computing 2nd power */
	r__1 = t * u;
	r = r * (r__1 * r__1) + (double)1.;
	t = s;
	goto L540;
L530:
/* Computing 2nd power */
	r__1 = s * u;
	r += r__1 * r__1;
L540:
	;
    }
    s = t * sqrt(r);
    r = a[k + j * a_dim1];
    u = (double)1. / sqrt(s * (s + fabs(r)));
    if (r < (double)0.) {
	s = -(double)s;
    }
    b[j] = -(double)s;
    a[k + j * a_dim1] = u * (r + s);
    i__1 = *m;
    for (i = h; i <= i__1; ++i) {
/* L550: */
	a[i + j * a_dim1] *= u;
    }
    if (*iq == 0) {
	goto L570;
    }
    i__1 = *m;
    for (i = k; i <= i__1; ++i) {
/* L560: */
	q[i + j * q_dim1] = a[i + j * a_dim1];
    }
L570:
    i__1 = *n;
    for (l = k; l <= i__1; ++l) {
	t = (double)0.;
	i__2 = *m;
	for (i = k; i <= i__2; ++i) {
/* L580: */
	    t += a[i + j * a_dim1] * a[i + l * a_dim1];
	}
	a[k + l * a_dim1] -= t * a[k + j * a_dim1];
	d[l] = a[k + l * a_dim1];
	i__2 = *m;
	for (i = h; i <= i__2; ++i) {
/* L590: */
	    a[i + l * a_dim1] -= t * a[i + j * a_dim1];
	}
/* L600: */
    }
    goto L350;
L610:
    d[*n] = a[*n + *n * a_dim1];
    b[*n - 1] = a[*n - 1 + *n * a_dim1];
L620:
    if (jq == 0) {
	goto L650;
    }
    if (*n > *m) {
	goto L640;
    }
    if (*n == *m) {
	goto L630;
    }
    if (jq == 1) {
	hsr3_(&q[q_offset], lq, m, n);
    }
    if (jq == 2) {
	hsr4_(&q[q_offset], lq, m, n);
    }
    if (jq == 3) {
	hsr5_(&q[q_offset], lq, m, n);
    }
    goto L650;
L630:
    hsr2_(&q[q_offset], lq, m);
    goto L650;
L640:
    hsr1_(&q[q_offset], lq, m);
L650:
    if (jp == 0) {
	return 0;
    }
    if (*n <= *m) {
	goto L660;
    }
    if (jp == 1) {
	hsr3_(&p[p_offset], lp, n, m);
    }
    if (jp == 2) {
	hsr4_(&p[p_offset], lp, n, m);
    }
    if (jp == 3) {
	hsr5_(&p[p_offset], lp, n, m);
    }
    return 0;
L660:
    hsr1_(&p[p_offset], lp, n);
    return 0;
} /* bidag2_ */


static int
hsr1_(double *a, int *la, int *n)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i, j, k, l, m;
    static real s;

    /* Parameter adjustments */
    a_dim1 = *la;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    a[a_dim1 + 1] = (double)1.;
    if (*n == 1) {
	return 0;
    }
    a[(a_dim1 << 1) + 1] = (double)0.;
    if (*n > 2) {
	goto L10;
    }
    a[a_dim1 + 2] = (double)0.;
    a[(a_dim1 << 1) + 2] = (double)1.;
    return 0;
L10:
    l = *n - 2;
    m = *n - 1;
    k = *n;
    s = a[k + l * a_dim1];
    a[*n + *n * a_dim1] = (double)1. - s * a[*n + l * a_dim1];
    a[m + *n * a_dim1] = -(double)s * a[m + l * a_dim1];
L20:
    j = m;
    m = l;
    --l;
    if (l == 0) {
	goto L50;
    }
    s = (double)0.;
    i__1 = *n;
    for (i = j; i <= i__1; ++i) {
/* L30: */
	s += a[i + l * a_dim1] * a[i + k * a_dim1];
    }
    a[m + k * a_dim1] = -(double)s * a[m + l * a_dim1];
    i__1 = *n;
    for (i = j; i <= i__1; ++i) {
/* L40: */
	a[i + k * a_dim1] -= s * a[i + l * a_dim1];
    }
    goto L20;
L50:
    a[k * a_dim1 + 1] = (double)0.;
    --k;
    m = k;
    l = k - 1;
    s = -(double)a[m + l * a_dim1];
    i__1 = *n;
    for (i = k; i <= i__1; ++i) {
/* L60: */
	a[i + k * a_dim1] = s * a[i + l * a_dim1];
    }
    a[k + k * a_dim1] += (double)1.;
    if (l > 1) {
	goto L20;
    }
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
/* L70: */
	a[i + a_dim1] = (double)0.;
    }
    return 0;
} /* hsr1_ */

static int
hsr2_(double *a, int *la, int *n)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i, j, k, l, m;
    static real s;

    /* Parameter adjustments */
    a_dim1 = *la;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    if (*n > 1) {
	goto L10;
    }
    a[a_dim1 + 1] = (double)1.;
    return 0;
L10:
    l = *n - 2;
    m = *n - 1;
    k = *n;
    s = a[*n + m * a_dim1];
    a[*n + *n * a_dim1] = (double)1. - s * a[*n + m * a_dim1];
    a[m + *n * a_dim1] = -(double)s * a[m + m * a_dim1];
L20:
    j = m;
    --m;
    if (m == 0) {
	goto L50;
    }
    s = (double)0.;
    i__1 = *n;
    for (i = j; i <= i__1; ++i) {
/* L30: */
	s += a[i + m * a_dim1] * a[i + k * a_dim1];
    }
    a[m + k * a_dim1] = -(double)s * a[m + m * a_dim1];
    i__1 = *n;
    for (i = j; i <= i__1; ++i) {
/* L40: */
	a[i + k * a_dim1] -= s * a[i + m * a_dim1];
    }
    goto L20;
L50:
    --k;
    m = k;
    s = -(double)a[k + k * a_dim1];
    i__1 = *n;
    for (i = k; i <= i__1; ++i) {
/* L60: */
	a[i + k * a_dim1] = s * a[i + k * a_dim1];
    }
    a[k + k * a_dim1] += (double)1.;
    if (k > 1) {
	goto L20;
    }
    return 0;
} /* hsr2_ */

static int
hsr3_(double *a, int *la, int *m, int *n)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i, j, k, l;
    static real s;

    /* Parameter adjustments */
    a_dim1 = *la;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    k = *n;
    if (*m >= *n) {
	goto L10;
    }
    psi::outfile->Printf("ERROR: ARGUMENT M MUST BE .GE. N IN SUBROUTINE HSR3\n");
    abort();
L10:
    s = -(double)a[k + k * a_dim1];
    i__1 = *m;
    for (i = k; i <= i__1; ++i) {
/* L20: */
	a[i + k * a_dim1] = s * a[i + k * a_dim1];
    }
    a[k + k * a_dim1] += (double)1.;
    if (k == 1) {
	return 0;
    }
    l = k;
L30:
    j = l;
    --l;
    if (l == 0) {
	goto L60;
    }
    s = (double)0.;
    i__1 = *m;
    for (i = j; i <= i__1; ++i) {
/* L40: */
	s += a[i + l * a_dim1] * a[i + k * a_dim1];
    }
    a[l + k * a_dim1] = -(double)s * a[l + l * a_dim1];
    i__1 = *m;
    for (i = j; i <= i__1; ++i) {
/* L50: */
	a[i + k * a_dim1] -= s * a[i + l * a_dim1];
    }
    goto L30;
L60:
    --k;
    goto L10;
} /* hsr3_ */

static int
hsr4_(double *a, int *la, int *m, int *n)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, k, l;
    static real s;

    /* Parameter adjustments */
    a_dim1 = *la;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    k = *m;
    if (*m > *n) {
	goto L10;
    }
    psi::outfile->Printf("ERROR: ARGUMENT M MUST BE .GE. N IN SUBROUTINE HSR4\n");
    abort();
L10:
    s = -(double)a[k + *n * a_dim1];
    i__1 = *m;
    for (i = *n; i <= i__1; ++i) {
/* L20: */
	a[i + k * a_dim1] = s * a[i + *n * a_dim1];
    }
    a[k + k * a_dim1] += (double)1.;
    l = *n;
L30:
    j = l;
    --l;
    if (l == 0) {
	goto L60;
    }
    s = (double)0.;
    i__1 = *m;
    for (i = j; i <= i__1; ++i) {
/* L40: */
	s += a[i + l * a_dim1] * a[i + k * a_dim1];
    }
    a[l + k * a_dim1] = -(double)s * a[l + l * a_dim1];
    i__1 = *m;
    for (i = j; i <= i__1; ++i) {
/* L50: */
	a[i + k * a_dim1] -= s * a[i + l * a_dim1];
    }
    goto L30;
L60:
    --k;
    if (k > *n) {
	goto L10;
    }
    k = *m - *n;
    i__1 = k;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
/* L70: */
	    a[i + j * a_dim1] = a[i + (j + *n) * a_dim1];
	}
    }
    return 0;
} /* hsr4_ */

/*      ________________________________________________________ */
/*     |                                                        | */
/*     |    PACKAGE SUBROUTINES: HSR3                           | */
/*     |________________________________________________________| */

static int
hsr5_(double *a, int *la, int *m, int *n)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i, j, k, l;
    static real s;

    /* Parameter adjustments */
    a_dim1 = *la;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    if (*m >= *n) {
	goto L10;
    }
    psi::outfile->Printf("ERROR: ARGUMENT M MUST BE .GE. N IN SUBROUTINE HSR5\n");
    abort();
L10:
    if (*m == *n) {
	goto L90;
    }
    k = *m;
L20:
    s = -(double)a[k + *n * a_dim1];
    i__1 = *m;
    for (i = *n; i <= i__1; ++i) {
/* L30: */
	a[i + k * a_dim1] = s * a[i + *n * a_dim1];
    }
    a[k + k * a_dim1] += (double)1.;
    l = *n;
L40:
    j = l;
    --l;
    if (l == 0) {
	goto L70;
    }
    s = (double)0.;
    i__1 = *m;
    for (i = j; i <= i__1; ++i) {
/* L50: */
	s += a[i + l * a_dim1] * a[i + k * a_dim1];
    }
    a[l + k * a_dim1] = -(double)s * a[l + l * a_dim1];
    i__1 = *m;
    for (i = j; i <= i__1; ++i) {
/* L60: */
	a[i + k * a_dim1] -= s * a[i + l * a_dim1];
    }
    goto L40;
L70:
    --k;
    if (k > *n) {
	goto L20;
    }
L90:
    hsr3_(&a[a_offset], la, m, n);
    return 0;
} /* hsr5_ */

/*      ________________________________________________________ */
/*     |                                                        | */
/*     |  COMPUTE SINGULAR VALUE DECOMPOSITION OF A BIDIAGONAL  | */
/*     | MATRIX:  A = Q TIMES DIAGONAL MATRIX TIMES P TRANSPOSE | */
/*     |                                                        | */
/*     |    INPUT:                                              | */
/*     |                                                        | */
/*     |         D     --DIAGONAL                               | */
/*     |                                                        | */
/*     |         N     --DIMENSION OF BIDIAGONAL MATRIX         | */
/*     |                                                        | */
/*     |         U     --OFF DIAGONAL                           | */
/*     |                                                        | */
/*     |         IU    --= 0 IF U IS SUPERDIAGONAL AND = 1 IF   | */
/*     |                 U IS SUBDIAGONAL                       | */
/*     |                                                        | */
/*     |         Q     --EITHER A SEGMENT OF AN IDENTITY MATRIX | */
/*     |                 OR A SEGMENT OF THE MATRIX USED TO     | */
/*     |                 BIDIAGONALIZE THE COEFFICIENT MATRIX   | */
/*     |                                                        | */
/*     |         LQ    --LEADING (ROW) DIMENSION OF ARRAY Q     | */
/*     |                                                        | */
/*     |         MQ    --NUMBER OF MATRIX ELEMENTS CONTAINED IN | */
/*     |                 EACH COLUMN OF Q                       | */
/*     |                                                        | */
/*     |         IQ    --AN INTEGER WHICH INDICATES WHETHER     | */
/*     |                 COLUMNS OF Q TO BE PROCESSED (= 0 MEANS| */
/*     |                 NO AND = 1 MEANS YES)                  | */
/*     |                                                        | */
/*     |         LP    --LEADING (ROW) DIMENSION OF ARRAY P     | */
/*     |                                                        | */
/*     |         MP    --NUMBER OF MATRIX ELEMENTS CONTAINED IN | */
/*     |                 EACH COLUMN OF P                       | */
/*     |                                                        | */
/*     |         IP    --AN INTEGER DEFINED LIKE IQ             | */
/*     |                                                        | */
/*     |         E     --WORK ARRAY WITH AT LEAST N ELEMENTS    | */
/*     |                                                        | */
/*     |         F     --WORK ARRAY WITH AT LEAST N ELEMENTS    | */
/*     |                                                        | */
/*     |    OUTPUT:                                             | */
/*     |                                                        | */
/*     |         Q     --Q FACTOR IN THE SINGULAR VALUE DECOMP. | */
/*     |                                                        | */
/*     |         D     --SINGULAR VALUES IN DECREASING ORDER    | */
/*     |                                                        | */
/*     |         P     --P FACTOR IN THE SINGULAR VALUE DECOMP. | */
/*     |                                                        | */
/*     |    BUILTIN FUNCTIONS: ABS,AMAX1,SIGN,SQRT              | */
/*     |    PACKAGE SUBROUTINES: EIG3,FGV,HSR1-HSR5,SCL,SFT,    | */
/*     |                         SNG0,SNG1,SORT2                | */
/*     |________________________________________________________| */

static int
singb_(double *d, int *n, double *u, int *iu,
       double *q, int *lq, int *mq, int *iq,
       double *p, int *lp, int *mp, int *ip, double *e, double *f)
{
    int one = 1;

    /* System generated locals */
    integer q_dim1, q_offset, p_dim1, p_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static real b, c;
    static integer g, h, i, j, k, l, m;
    static real r, s, t, v, w, x, y, z;
    static integer j0, k2, l1;
    static real t0, t1, t2, t3;
    static integer id, jl, ll, jr, ns;

    /* Parameter adjustments */
    --f;
    --e;
    p_dim1 = *lp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    q_dim1 = *lq;
    q_offset = q_dim1 + 1;
    q -= q_offset;
    --u;
    --d;

    /* Function Body */
    if (*n > 1) {
	goto L10;
    }
    if (*iq == 1) {
	q[q_dim1 + 1] = (double)1.;
    }
    if (*ip == 1) {
	p[p_dim1 + 1] = (double)1.;
    }
    if (d[1] >= (double)0.) {
	return 0;
    }
    d[1] = -(double)d[1];
    if (*iq == 1) {
	q[q_dim1 + 1] = (double)-1.;
    }
    return 0;
L10:
    jl = *iq;
    if (jl == 0) {
	jl = 3;
    }
    if (jl != 3) {
	jl = 1;
    }
    jr = *ip;
    if (jr == 0) {
	jr = 3;
    }
    if (jr != 3) {
	jr = 2;
    }
    if (*iu == 0) {
	goto L20;
    }
    i = jl;
    jl = jr;
    jr = i;
L20:
    j = 0;
    l = *n - 1;
    k2 = *n - 2;
    i__1 = l;
    for (i = 1; i <= i__1; ++i) {
	e[i] = (double)1.;
	f[i] = (double)1.;
	if (u[i] == (double)0.) {
	    j = i;
	}
/* L30: */
	if (d[i] == (double)0.) {
	    j = i;
	}
    }
    e[*n] = (double)1.;
    f[*n] = (double)1.;
    b = (double)3.5527136788005009e-15;
    t = (double)1.;
L40:
    t *= (double).5;
    s = t + (double)1.;
    if (s > (double)1.) {
	goto L40;
    }
    t0 = (double)1. / (t + t);
/* Computing 2nd power */
    r__1 = t + t;
    t2 = r__1 * r__1;
    ns = *n * 50;
    l1 = 0;
    k = *n;
    ll = 0;
    goto L70;
L50:
    j = 0;
    i__1 = l;
    for (i = 1; i <= i__1; ++i) {
	if (u[i] == (double)0.) {
	    j = i;
	}
/* L60: */
	if (d[i] == (double)0.) {
	    j = i;
	}
    }
L70:
    if (j == 0) {
	goto L140;
    }
    if (u[j] == (double)0.) {
	goto L140;
    }
/*     ------------------------------- */
/*     |*** ZERO DIAGONAL ELEMENT ***| */
/*     ------------------------------- */
    i = j;
    v = u[j];
    u[j] = (double)0.;
    s = -(double)v * (r__1 = e[j], fabs(r__1));
/*     --------------------------- */
/*     |*** PROCESS LEFT SIDE ***| */
/*     --------------------------- */
L80:
    h = i;
    ++i;
    r = (r__1 = e[i], fabs(r__1)) * d[i];
    fgv_(&x, &y, &t, &r, &s, &e[j], &e[i]);
    if (t == (double)1.) {
	goto L90;
    }
    d[i] -= y * v;
    sng0_(&q[q_offset], lq, mq, &p[p_offset], lp, mp, &jl, &j, &i, &x, &y);
    if (i == k) {
	goto L100;
    }
    v = x * u[i];
    s = -(double)v * (r__1 = e[j], fabs(r__1));
    goto L80;
L90:
    d[i] = d[i] * y - v;
    sng1_(&q[q_offset], lq, mq, &p[p_offset], lp, mp, &jl, &j, &i, &x, &y);
    if (i == k) {
	goto L100;
    }
    v = u[i];
    u[i] = v * y;
    s = -(double)v * (r__1 = e[j], fabs(r__1));
    goto L80;
L100:
    if (j == 1) {
	goto L130;
    }
    i = j - 1;
    s = u[i];
    u[i] = (double)0.;
/*     ---------------------------- */
/*     |*** PROCESS RIGHT SIDE ***| */
/*     ---------------------------- */
L110:
    h = i;
    --i;
    r = d[h];
    fgv_(&x, &y, &t, &r, &s, &f[h], &f[j]);
    if (t == (double)1.) {
	goto L120;
    }
    d[h] = r + x * s;
    sng0_(&q[q_offset], lq, mq, &p[p_offset], lp, mp, &jr, &h, &j, &x, &y);
    if (h == 1) {
	goto L130;
    }
    s = -(double)y * u[i];
    goto L110;
L120:
    d[h] = x * r + s;
    sng1_(&q[q_offset], lq, mq, &p[p_offset], lp, mp, &jr, &h, &j, &x, &y);
    if (h == 1) {
	goto L130;
    }
    s = -(double)u[i];
    u[i] = x * u[i];
    goto L110;
L130:
    scl_(&d[1], &u[1], n, &q[q_offset], lq, mq, &p[p_offset], lp, mp, &e[1], &
	    f[1], &b, &one, &k, &jl, &jr);
L140:
    ++j;
    if (j == k) {
	goto L320;
    }
    s = (double)0.;
    t = (double)0.;
/*     ----------------------------- */
/*     |*** SET ERROR TOLERANCE ***| */
/*     ----------------------------- */
    i__1 = l;
    for (i = j; i <= i__1; ++i) {
	x = (r__1 = d[i] * e[i] * (d[i] * f[i]), fabs(r__1));
	y = (r__1 = u[i] * e[i] * (u[i] * f[i + 1]), fabs(r__1));
/* Computing MAX */
	r__1 = max(s,x);
	s = dmax(r__1,y);
/* L150: */
	t = t + x + y;
    }
    x = (r__1 = e[k] * d[k] * (f[k] * d[k]), fabs(r__1));
    s = dmax(s,x);
    t3 = s * t2;
    t += x;
    if (t == (double)0.) {
	t = (double)1.;
    }
    t1 = t0 / t;
    goto L280;
L160:
    ++ll;
    if (ll > ns) {
	goto L530;
    }
    if (l > j) {
	goto L170;
    }
    s = (double)0.;
    t = (double)0.;
    goto L180;
L170:
    s = u[k2];
    t = e[k2];
/*     ----------------------- */
/*     |*** COMPUTE SHIFT ***| */
/*     ----------------------- */
L180:
    sft_(&c, &d[k], &d[l], &u[l], &s, &e[k], &e[l], &t, &f[k], &f[l]);
    v = (r__1 = e[j], fabs(r__1));
    w = (r__1 = f[j], fabs(r__1));
    t = d[j] * w;
    r = t * (d[j] * v) - c;
    s = t * (u[j] * v);
    id = 1;
    j0 = j;
    i = j;
    h = j - 1;
L190:
    g = h;
    h = i;
    ++i;
    z = (r__1 = f[i], fabs(r__1));
/*     ---------------------------- */
/*     |*** PROCESS RIGHT SIDE ***| */
/*     ---------------------------- */
    fgv_(&x, &y, &t, &r, &s, &f[h], &f[i]);
    v = d[h];
    if (t == (double)1.) {
	goto L210;
    }
    if (h == j) {
	goto L200;
    }
    t = u[g] + x * s;
    u[g] = t;
    if ((r__1 = t * e[g] * (t * f[h]), fabs(r__1)) > t3) {
	goto L200;
    }
    j = h;
    id = 1;
L200:
    r = v + x * u[h];
    s = x * d[i];
    u[h] -= y * v;
    sng0_(&q[q_offset], lq, mq, &p[p_offset], lp, mp, &jr, &h, &i, &x, &y);
    goto L230;
L210:
    if (h == j) {
	goto L220;
    }
    t = x * u[g] + s;
    u[g] = t;
    if ((r__1 = t * e[g] * (t * f[h]), fabs(r__1)) > t3) {
	goto L220;
    }
    j = h;
    id = 1;
L220:
    r = x * v + u[h];
    s = d[i];
    u[h] = y * u[h] - v;
    d[i] = y * s;
    sng1_(&q[q_offset], lq, mq, &p[p_offset], lp, mp, &jr, &h, &i, &x, &y);
/*     --------------------------- */
/*     |*** PROCESS LEFT SIDE ***| */
/*     --------------------------- */
L230:
    fgv_(&x, &y, &t, &r, &s, &e[h], &e[i]);
    if (t == (double)1.) {
	goto L250;
    }
    t = r + x * s;
    d[h] = t;
    if ((r__1 = t * e[h] * (t * f[h]), fabs(r__1)) > t3) {
	goto L240;
    }
    id = 0;
    j = h;
L240:
    r = u[h] + x * d[i];
    d[i] -= y * u[h];
    u[h] = r;
    sng0_(&q[q_offset], lq, mq, &p[p_offset], lp, mp, &jl, &h, &i, &x, &y);
    if (i == k) {
	goto L270;
    }
    s = x * u[i];
    goto L190;
L250:
    t = s + x * r;
    d[h] = t;
    if ((r__1 = t * e[h] * (t * f[h]), fabs(r__1)) > t3) {
	goto L260;
    }
    id = 0;
    j = h;
L260:
    r = d[i] + x * u[h];
    d[i] = y * d[i] - u[h];
    u[h] = r;
    sng1_(&q[q_offset], lq, mq, &p[p_offset], lp, mp, &jl, &h, &i, &x, &y);
    if (i == k) {
	goto L270;
    }
    s = u[i];
    u[i] = s * y;
    goto L190;
L270:
    scl_(&d[1], &u[1], n, &q[q_offset], lq, mq, &p[p_offset], lp, mp, &e[1], &
	    f[1], &b, &j0, &k, &jl, &jr);
    if (id == 0) {
	goto L70;
    }
L280:
    w = e[l];
    x = e[k];
    y = f[l];
    z = f[k];
    r = (r__1 = x * d[k] * (z * d[k]), fabs(r__1));
    s = (r__1 = w * d[l] * (y * d[l]), fabs(r__1));
    t = (r__1 = x * u[l] * (y * u[l]), fabs(r__1));
/*     ------------------------------ */
/*     |*** TEST FOR CONVERGENCE ***| */
/*     ------------------------------ */
    if (s * t1 * (t * t1) > (double)1.) {
	goto L160;
    }
    ++l1;
    if (l1 > 40) {
	goto L290;
    }
    if (t == (double)0.) {
	goto L290;
    }
    r += t;
    if (s / r * (t / r) > t2) {
	goto L160;
    }
L290:
    l1 = 0;
    if (s > r) {
	goto L310;
    }
    r = -(double)d[k] * fabs(x);
    s = u[l] * fabs(w);
    fgv_(&w, &y, &t, &r, &s, &e[l], &e[k]);
    x = e[k];
    if (t == (double)1.) {
	goto L300;
    }
    d[k] -= y * u[l];
    sng0_(&q[q_offset], lq, mq, &p[p_offset], lp, mp, &jl, &l, &k, &w, &y);
    goto L310;
L300:
    d[l] *= w;
    d[k] = y * d[k] - u[l];
    sng1_(&q[q_offset], lq, mq, &p[p_offset], lp, mp, &jl, &l, &k, &w, &y);
L310:
    r__1 = sqrt((fabs(x)));
    r__2 = sqrt((fabs(z)));
    t = r_sign(&r__1, &x) * d[k] * r_sign(&r__2, &z);
    if (t < (double)0.) {
	e[k] = -(double)e[k];
    }
    d[k] = fabs(t);
    k = l;
    l = k2;
    --k2;
    if (k > j) {
	goto L280;
    }
L320:
    x = e[k];
    z = f[k];
    r__1 = sqrt((fabs(x)));
    r__2 = sqrt((fabs(z)));
    t = r_sign(&r__1, &x) * d[k] * r_sign(&r__2, &z);
    if (t < (double)0.) {
	e[k] = -(double)e[k];
    }
    d[k] = fabs(t);
    k = l;
    l = k2;
    --k2;
    if (k > 1) {
	goto L50;
    }
    if (k == 0) {
	goto L330;
    }
    goto L320;
/*     ------------------------- */
/*     |*** RESCALE P AND Q ***| */
/*     ------------------------- */
L330:
    switch ((int)jl) {
	case 1:  goto L340;
	case 2:  goto L360;
	case 3:  goto L380;
    }
L340:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	t = e[j];
	r__1 = sqrt((fabs(t)));
	t = r_sign(&r__1, &t);
	i__2 = *mq;
	for (i = 1; i <= i__2; ++i) {
/* L350: */
	    q[i + j * q_dim1] *= t;
	}
    }
    goto L380;
L360:
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
	t = e[j];
	r__1 = sqrt((fabs(t)));
	t = r_sign(&r__1, &t);
	i__1 = *mp;
	for (i = 1; i <= i__1; ++i) {
/* L370: */
	    p[i + j * p_dim1] *= t;
	}
    }
L380:
    switch ((int)jr) {
	case 1:  goto L390;
	case 2:  goto L410;
	case 3:  goto L430;
    }
L390:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	t = f[j];
	r__1 = sqrt((fabs(t)));
	t = r_sign(&r__1, &t);
	i__2 = *mq;
	for (i = 1; i <= i__2; ++i) {
/* L400: */
	    q[i + j * q_dim1] *= t;
	}
    }
    goto L430;
L410:
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
	t = f[j];
	r__1 = sqrt((fabs(t)));
	t = r_sign(&r__1, &t);
	i__1 = *mp;
	for (i = 1; i <= i__1; ++i) {
/* L420: */
	    p[i + j * p_dim1] *= t;
	}
    }
/*     -------------------------------- */
/*     |*** REORDER THE EIGENPAIRS ***| */
/*     -------------------------------- */
L430:
    sort2_(&d[1], &e[1], &f[1], n);
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	j = e[i];
/* L440: */
	f[j] = (real) i;
    }
    m = *n + 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = m - j;
	k = e[l];
	i = f[j];
	e[i] = (real) k;
	f[k] = (real) i;
	t = d[j];
	d[j] = d[k];
	d[k] = t;
	if (jl == 1) {
	    goto L450;
	}
	if (jr != 1) {
	    goto L480;
	}
L450:
	s = (double)0.;
	i__2 = *mq;
	for (i = 1; i <= i__2; ++i) {
	    t = q[i + k * q_dim1];
	    q[i + k * q_dim1] = q[i + j * q_dim1];
	    q[i + j * q_dim1] = t;
/* L460: */
	    s += t * t;
	}
	s = (double)1. / sqrt(s);
	i__2 = *mq;
	for (i = 1; i <= i__2; ++i) {
/* L470: */
	    q[i + j * q_dim1] = s * q[i + j * q_dim1];
	}
L480:
	if (jr == 2) {
	    goto L490;
	}
	if (jl != 2) {
	    goto L520;
	}
L490:
	s = (double)0.;
	i__2 = *mp;
	for (i = 1; i <= i__2; ++i) {
	    t = p[i + k * p_dim1];
	    p[i + k * p_dim1] = p[i + j * p_dim1];
	    p[i + j * p_dim1] = t;
/* L500: */
	    s += t * t;
	}
	s = (double)1. / sqrt(s);
	i__2 = *mp;
	for (i = 1; i <= i__2; ++i) {
/* L510: */
	    p[i + j * p_dim1] = s * p[i + j * p_dim1];
	}
L520:
	;
    }
    e[1] = (real) (*n);
    return 0;
L530:
    k = *n - k + 1;
    psi::outfile->Printf("SINCE THE STOPPING CRITERION NOT SATISFIED\n");
    psi::outfile->Printf("AFTER %d ITERATIONS, WE STOP WHILE COMPUTING\n", ns);
    psi::outfile->Printf("EIGENVALUE NUMBER %d", k);
    e[1] = (real) k;
    return 0;
} /* singb_ */

/* % */
static int
fgv_(double *x, double *y, double *s, double *p, double *q,
     double *a, double *b)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    static real c, r, t;

    if (fabs(*p) > fabs(*q)) {
	goto L10;
    }
    if (*q == (double)0.) {
	goto L110;
    }
    r = *a / *b;
    *s = *p / *q;
    t = fabs(r) * *s * *s;
    if (t < (double)1.) {
	goto L70;
    }
    t /= t + (double)1.;
    r = *b / *a;
    *s = *q / *p;
    goto L20;
L10:
    r = *b / *a;
    *s = *q / *p;
    t = fabs(r) * *s * *s;
    if (t > (double)1.) {
	goto L60;
    }
    t = (double)1. / (t + (double)1.);
L20:
    r__1 = *a * t;
    *a = r_sign(&r__1, p);
    if (r_sign(&r, p) == r) {
	goto L40;
    }
    *b = -(double)(r__1 = *b * t, fabs(r__1));
    goto L50;
L40:
    *b = (r__1 = *b * t, fabs(r__1));
L50:
    *y = *s;
    *x = fabs(r) * *s;
    *s = (double)0.;
    return 0;
L60:
    t /= t + (double)1.;
    r = *a / *b;
    *s = *p / *q;
    goto L80;
L70:
    t = (double)1. / (t + (double)1.);
L80:
    c = *a;
    r__1 = *b * t;
    *a = r_sign(&r__1, q);
    if (r_sign(&r, q) == r) {
	goto L90;
    }
    *b = -(double)(r__1 = c * t, fabs(r__1));
    goto L100;
L90:
    *b = (r__1 = c * t, fabs(r__1));
L100:
    *y = *s;
    *x = fabs(r) * *s;
    *s = (double)1.;
    return 0;
L110:
    *x = (double)0.;
    *y = (double)0.;
    *s = (double)0.;
    return 0;
} /* fgv_ */

/* % */
static int
sng0_(double *q, int *lq, int *m, double *p, int *lp, int *n,
      int *l, int *j, int *k, double *x, double *y)
{
    /* System generated locals */
    integer q_dim1, q_offset, p_dim1, p_offset, i__1;

    /* Local variables */
    static integer i;
    static real s, t;

    /* Parameter adjustments */
    p_dim1 = *lp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    q_dim1 = *lq;
    q_offset = q_dim1 + 1;
    q -= q_offset;

    /* Function Body */
    switch ((int)*l) {
	case 1:  goto L10;
	case 2:  goto L30;
	case 3:  goto L50;
    }
L10:
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	t = q[i + *j * q_dim1];
	s = q[i + *k * q_dim1];
	q[i + *j * q_dim1] = t + *x * s;
/* L20: */
	q[i + *k * q_dim1] = s - *y * t;
    }
    return 0;
L30:
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	t = p[i + *j * p_dim1];
	s = p[i + *k * p_dim1];
	p[i + *j * p_dim1] = t + *x * s;
/* L40: */
	p[i + *k * p_dim1] = s - *y * t;
    }
L50:
    return 0;
} /* sng0_ */

/* % */
static int
sng1_(double *q, int *lq, int *m, double *p,
      int *lp, int *n, int *l, int *j, int *k, double *x, double *y)
{
    /* System generated locals */
    integer q_dim1, q_offset, p_dim1, p_offset, i__1;

    /* Local variables */
    static integer i;
    static real s, t;

    /* Parameter adjustments */
    p_dim1 = *lp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    q_dim1 = *lq;
    q_offset = q_dim1 + 1;
    q -= q_offset;

    /* Function Body */
    switch ((int)*l) {
	case 1:  goto L10;
	case 2:  goto L30;
	case 3:  goto L50;
    }
L10:
    i__1 = *m;
    for (i = 1; i <= i__1; ++i) {
	t = q[i + *j * q_dim1];
	s = q[i + *k * q_dim1];
	q[i + *j * q_dim1] = *x * t + s;
/* L20: */
	q[i + *k * q_dim1] = *y * s - t;
    }
    return 0;
L30:
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	t = p[i + *j * p_dim1];
	s = p[i + *k * p_dim1];
	p[i + *j * p_dim1] = *x * t + s;
/* L40: */
	p[i + *k * p_dim1] = *y * s - t;
    }
L50:
    return 0;
} /* sng1_ */

/* % */
static int
sft_(double *s, double *a, double *b, double *c,
     double *d, double *e2, double *e1, double *e0, double *f2, double *f1)
{
    static real w, x, y, z, g0, g1, g2, h1, h2;

    g0 = fabs(*e0);
    g1 = fabs(*e1);
    g2 = fabs(*e2);
    h1 = fabs(*f1);
    h2 = fabs(*f2);
    w = *a * g2 * (*a * h2) + *c * g1 * (*c * h2);
    x = *b * g1 * (*b * h1) + *d * g0 * (*d * h1);
    y = *b * g1 * (*c * h1);
    z = *b * g1 * (*c * h2);
    eig3_(s, s, &x, &w, &y, &z);
    return 0;
} /* sft_ */

/* % */
static int
scl_(double *d, double *u, int */*n*/, double *q,
     int *lq, int *mq, double *p, int *lp, int *mp,
     double *e, double *f, double *b,
     int *j, int *k, int *jl,
     int *jr)
{
    /* System generated locals */
    integer p_dim1, p_offset, q_dim1, q_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer h, i, l, m;
    static real r, s, t, v;

    /* Parameter adjustments */
    --f;
    --e;
    p_dim1 = *lp;
    p_offset = p_dim1 + 1;
    p -= p_offset;
    q_dim1 = *lq;
    q_offset = q_dim1 + 1;
    q -= q_offset;
    --u;
    --d;

    /* Function Body */
    t = (double)1.;
    i__1 = *k;
    for (i = *j; i <= i__1; ++i) {
	if ((r__1 = e[i], fabs(r__1)) < t) {
	    t = (r__2 = e[i], fabs(r__2));
	}
/* L10: */
	if ((r__1 = f[i], fabs(r__1)) < t) {
	    t = (r__2 = f[i], fabs(r__2));
	}
    }
    if (t > *b) {
	return 0;
    }
/*     ---------------------------- */
/*     |*** RESCALE THE MATRIX ***| */
/*     ---------------------------- */
    r = e[*j];
    r__1 = sqrt((fabs(r)));
    r = r_sign(&r__1, &r);
    e[*j] = (double)1.;
    s = f[*j];
    r__1 = sqrt((fabs(s)));
    s = r_sign(&r__1, &s);
    f[*j] = (double)1.;
    d[*j] = r * d[*j] * s;
    t = r;
    if (*jl == 1) {
	goto L20;
    }
    if (*jr != 1) {
	goto L40;
    }
    t = s;
L20:
    i__1 = *mq;
    for (i = 1; i <= i__1; ++i) {
/* L30: */
	q[i + *j * q_dim1] *= t;
    }
L40:
    t = s;
    if (*jr == 2) {
	goto L50;
    }
    if (*jl != 2) {
	goto L70;
    }
    t = r;
L50:
    i__1 = *mp;
    for (i = 1; i <= i__1; ++i) {
/* L60: */
	p[i + *j * p_dim1] *= t;
    }
L70:
    l = *j + 1;
    h = *j;
    i__1 = *k;
    for (m = l; m <= i__1; ++m) {
	t = e[m];
	r__1 = sqrt((fabs(t)));
	t = r_sign(&r__1, &t);
	e[m] = (double)1.;
	s = f[m];
	r__1 = sqrt((fabs(s)));
	s = r_sign(&r__1, &s);
	f[m] = (double)1.;
	d[m] = s * d[m] * t;
	u[h] = r * u[h] * s;
	h = m;
	r = t;
	v = t;
	if (*jl == 1) {
	    goto L80;
	}
	if (*jr != 1) {
	    goto L100;
	}
	v = s;
L80:
	i__2 = *mq;
	for (i = 1; i <= i__2; ++i) {
/* L90: */
	    q[i + m * q_dim1] *= v;
	}
L100:
	v = s;
	if (*jr == 2) {
	    goto L110;
	}
	if (*jl != 2) {
	    goto L130;
	}
	v = t;
L110:
	i__2 = *mp;
	for (i = 1; i <= i__2; ++i) {
/* L120: */
	    p[i + m * p_dim1] *= v;
	}
L130:
	;
    }
    return 0;
} /* scl_ */

static int
eig3_(double *ea, double *eb, double *a, double *b, double *y, double *z)
{
    /* Local variables */
    static real c, s, t;

    t = (*b - *a) * (double).5;
    c = sqrt((fabs(*y))) * sqrt((fabs(*z)));
    if (fabs(t) > fabs(c)) {
	goto L30;
    }
    if (c != (double)0.) {
	goto L10;
    }
    *ea = *a;
    *eb = *b;
    return 0;
L10:
    t /= fabs(c);
    s = fabs(c) / (fabs(t) + sqrt(t * t + (double)1.));
    if (t < (double)0.) {
	goto L20;
    }
    *ea = *a - s;
    *eb = *b + s;
    return 0;
L20:
    *ea = *a + s;
    *eb = *b - s;
    return 0;
L30:
    t = fabs(c) / t;
    s = t * fabs(c) / (sqrt(t * t + (double)1.) + (double)1.);
    *ea = *a - s;
    *eb = *b + s;
    return 0;
} /* eig3_ */

/*      ________________________________________________________ */
/*     |                                                        | */
/*     |            SORT AN ARRAY IN INCREASING ORDER           | */
/*     |                                                        | */
/*     |    INPUT:                                              | */
/*     |                                                        | */
/*     |         X     --ARRAY OF NUMBERS                       | */
/*     |                                                        | */
/*     |         W     --WORKING ARRAY (LENGTH  AT LEAST N)     | */
/*     |                                                        | */
/*     |         N     --NUMBER OF ARRAY ELEMENTS TO SORT       | */
/*     |                                                        | */
/*     |    OUTPUT:                                             | */
/*     |                                                        | */
/*     |         X     --ORIGINAL ARRAY                         | */
/*     |                                                        | */
/*     |         Y     --INDICES OF X GIVING INCREASING ORDER   | */
/*     |________________________________________________________| */

static int
sort2_(double *x, double *y, double *w, int *n)
{
    static integer i, j, k, l, m, p, q;
    static real s, t;

    /* Parameter adjustments */
    --w;
    --y;
    --x;

    /* Function Body */
    i = 1;
L10:
    k = i;
L20:
    j = i;
    y[i] = (real) i;
    ++i;
    if (j == *n) {
	goto L30;
    }
    if (x[i] >= x[j]) {
	goto L20;
    }
    w[k] = (real) i;
    goto L10;
L30:
    if (k == 1) {
	return 0;
    }
    w[k] = (real) (*n + 1);
L40:
    m = 1;
    l = 1;
L50:
    i = l;
    if (i > *n) {
	goto L120;
    }
    p = y[i];
    s = x[p];
    j = w[i];
    k = j;
    if (j > *n) {
	goto L100;
    }
    q = y[j];
    t = x[q];
    l = w[j];
    y[i] = (real) l;
L60:
    if (s > t) {
	goto L70;
    }
    w[m] = (real) p;
    ++m;
    ++i;
    if (i == k) {
	goto L80;
    }
    p = y[i];
    s = x[p];
    goto L60;
L70:
    w[m] = (real) q;
    ++m;
    ++j;
    if (j == l) {
	goto L110;
    }
    q = y[j];
    t = x[q];
    goto L60;
L80:
    w[m] = (real) q;
    k = m + l - j;
    i = j - m;
L90:
    ++m;
    if (m == k) {
	goto L50;
    }
    w[m] = y[m + i];
    goto L90;
L100:
    y[i] = (real) j;
    l = j;
L110:
    w[m] = (real) p;
    k = m + k - i;
    i -= m;
    goto L90;
L120:
    i = 1;
L130:
    k = i;
    j = y[i];
L140:
    y[i] = w[i];
    ++i;
    if (i < j) {
	goto L140;
    }
    w[k] = (real) i;
    if (i <= *n) {
	goto L130;
    }
    if (k > 1) {
	goto L40;
    }
    return 0;
} /* sort2_ */
