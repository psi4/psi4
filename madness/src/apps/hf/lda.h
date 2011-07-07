/*
  This file is part of MADNESS.
  
  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:
  
  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367
  
  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
  
  $Id: lda.h 1981 2010-10-12 03:01:55Z rjharrison $
*/
/*
 * lda.h
 *
 *  Created on: Nov 10, 2008
 *      Author: wsttiger
 */

#ifndef LDA_H_
#define LDA_H_

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <world/world.h>
#include <math.h>
#include <madness_config.h>

typedef double doublereal;
typedef MADNESS_FORINT integer;
typedef int logical;

#define max(a,b) ((a)>=(b)?(a):(b))

/* static double pow_dd(doublereal* a, doublereal* b) { */
/*     return pow(*a,*b); */
/* } */

using namespace madness;

/* lda.f -- translated by f2c (version 20050501).
   You must link the resulting object file with libf2c:
  on Microsoft Windows system, link with libf2c.lib;
  on Linux or Unix systems, link with .../path/to/libf2c.a -lm
  or, if you install libf2c.a in a standard place, with -lf2c -lm
  -- in that order, at the end of the command line, as in
    cc *.o -lf2c -lm
  Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

    http://www.netlib.org/f2c/libf2c.zip
*/


/* Table of constant values */

#include <cmath>

static inline double pow(const double* a, const double* b) {
    return pow(*a, *b);
}

static double c_b2 = .333333333333333333333333333333333;
static double c_b7 = .333333333333333333333333333333;
static double c_b8 = .5;
static double c_b14 = 1.333333333333333333333333333333;


inline /* Subroutine */ int x_rks_s__(const double *r__, double *f, double *
  dfdra)
{

    /* Local variables */
    double ra13;


/*     This subroutine evaluates the spin polarised exchange functional */
/*     in the Local Density Approximation [1], and the corresponding */
/*     potential. Often this functional is referred to as the Dirac */
/*     functional [2] or Slater functional. */

/*     [1] F. Bloch, Zeitschrift fuer Physik, Vol. 57 (1929) 545. */

/*     [2] P.A.M. Dirac, Proceedings of the Cambridge Philosophical */
/*         Society, Vol. 26 (1930) 376. */

/*     Parameters: */

/*     r     the total electron density */
/*     f     On return the functional value */
/*     dfdra On return the derivative of f with respect to alpha electron */
/*           density */


/*     Ax = -3/4*(6/pi)**(1/3) */
/*     Bx = -(6/pi)**(1/3) */
/*     C  = (1/2)**(1/3) */




    ra13 = pow(r__, &c_b2) * .793700525984099737375852819636154;
    *f = *r__ * -.930525736349100025002010218071667 * ra13;
    *dfdra = ra13 * -1.24070098179880003333601362409556;

    return 0;
} /* x_rks_s__ */

/* ----------------------------------------------------------------------- */
inline /* Subroutine */ int x_uks_s__(double *ra, double *rb, double *f,
  double *dfdra, double *dfdrb)
{
    /* Local variables */
    double ra13, rb13;


/*     This subroutine evaluates the spin polarised exchange functional */
/*     in the Local Density Approximation [1], and the corresponding */
/*     potential. Often this functional is referred to as the Dirac */
/*     functional [2] or Slater functional. */

/*     [1] F. Bloch, Zeitschrift fuer Physik, Vol. 57 (1929) 545. */

/*     [2] P.A.M. Dirac, Proceedings of the Cambridge Philosophical */
/*         Society, Vol. 26 (1930) 376. */

/*     Parameters: */

/*     ra    the alpha electron density */
/*     rb    the beta  electron density */
/*     f     On return the functional value */
/*     dfdra On return the derivative of f with respect to ra */
/*     dfdrb On return the derivative of f with respect to rb */


/*     Ax = -3/4*(6/pi)**(1/3) */
/*     Bx = -(6/pi)**(1/3) */




    ra13 = pow(ra, &c_b2);
    rb13 = pow(rb, &c_b2);
    *f = (*ra * ra13 + *rb * rb13) * -.930525736349100025002010218071667;
    *dfdra = ra13 * -1.24070098179880003333601362409556;
    *dfdrb = rb13 * -1.24070098179880003333601362409556;

    return 0;
} /* x_uks_s__ */

inline /* Subroutine */ int c_rks_vwn5__(const double *r__, double *f, double *
  dfdra)
{
    /* Local variables */
    double a2, b2, c2, d2, i1, i2, i3, p1, p2, p3, p4, t4, t5, t6,
      t7, iv, alpha_rho13__, iv2, pp1, pp2, inv, srho, srho13;


/*     This subroutine evaluates the Vosko, Wilk and Nusair correlation */
/*     functional number 5 [1] for the closed shell case, with the */
/*     parametrisation as given in table 5. */

/*     The original code was obtained from Dr. Phillip Young, */
/*     with corrections from Dr. Paul Sherwood. */

/*     [1] S.H. Vosko, L. Wilk, and M. Nusair */
/*         "Accurate spin-dependent electron liquid correlation energies */
/*          for local spin density calculations: a critical analysis", */
/*         Can.J.Phys, Vol. 58 (1980) 1200-1211. */

/*     Parameters: */

/*     r      the total electron density */
/*     f      On return the functional value */
/*     dfdra  On return the derivative of f with respect to the alpha */
/*            electron density */




/* VWN interpolation parameters */

/* paramagnetic */
    a2 = .0621814;
    b2 = 3.72744;
    c2 = 12.9352;
    d2 = -.10498;

/* t4 = (1/(4/3)*pi)**(1/3) */
    t4 = .620350490899399531;

/* t5 = 0.5/(2**(1/3)-1) */
    t5 = 1.92366105093153617;

/* t6 = 2/(3*(2**(1/3)-1)) */
    t6 = 2.56488140124204822;

/* t7 = 2.25*(2**(1/3)-1) */
    t7 = .584822362263464735;

/* Paramagnetic interpolation constants */

    p1 = 6.1519908197590798;
    p2 = a2 * .5;
    p3 = 9.6902277115443745e-4;
    p4 = .038783294878113009;

/* closed shell case */
    srho = *r__;
    srho13 = pow(&srho, &c_b7);
    alpha_rho13__ = pow(&c_b8, &c_b7) * srho;
    iv2 = t4 / srho13;
    iv = sqrt(iv2);

/* paramagnetic */
    inv = 1. / (iv2 + b2 * iv + c2);
    i1 = log(iv2 * inv);
    i2 = log((iv - d2) * (iv - d2) * inv);
/* corrected b1->b2 ps Apr98 */
    i3 = atan(p1 / (iv * 2. + b2));
    pp1 = p2 * i1 + p3 * i2 + p4 * i3;
    pp2 = a2 * (1. / iv - iv * inv * (b2 / (iv - d2) + 1.));

    *f = pp1 * srho;
    *dfdra = pp1 - iv * .166666666666666666666666666666 * pp2;

    return 0;
} /* c_rks_vwn5__ */

/* ----------------------------------------------------------------------- */
inline /* Subroutine */ int c_uks_vwn5__(double *ra, double *rb, double *
  f, double *dfdra, double *dfdrb)
{
    /* System generated locals */
    double d__1, d__2;

    /* Local variables */
    double v, beta_rho13__, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3,
       c3, d3, f1, f2, f3, p1, p2, p3, s1, t1, t2, s2, t4, t5, t6, t7,
      s3, s4, p4, f4, i1, i2, i3, iv, alpha_rho13__, ff1, ff2, iv2, pp1,
       pp2, ss1, ss2, tau, inv, vwn1, vwn2, dtau, zeta, srho, zeta3,
      zeta4, srho13, inter1, inter2;


/*     This subroutine evaluates the Vosko, Wilk and Nusair correlation */
/*     functional number 5 [1], with the parametrisation as given in */
/*     table 5. */

/*     The original code was obtained from Dr. Phillip Young, */
/*     with corrections from Dr. Paul Sherwood. */

/*     [1] S.H. Vosko, L. Wilk, and M. Nusair */
/*         "Accurate spin-dependent electron liquid correlation energies */
/*          for local spin density calculations: a critical analysis", */
/*         Can.J.Phys, Vol. 58 (1980) 1200-1211. */

/*     Parameters: */

/*     ra     the alpha-electron density */
/*     rb     the beta-electron density */
/*     f      On return the functional value */
/*     dfdra  On return the derivative of f with respect to ra */
/*     dfdrb  On return the derivative of f with respect to rb */



/*     tn13 = 2**(1/3) */
/*     tn43 = 2**(4/3) */

/* VWN interpolation parameters */

/* spin stiffness */
    a1 = -.0337737278807791058;
    b1 = 1.13107;
    c1 = 13.0045;
    d1 = -.0047584;
/* paramagnetic */
    a2 = .0621814;
    b2 = 3.72744;
    c2 = 12.9352;
    d2 = -.10498;
/* ferromagnetic */
/* try cadpac/nwchem value (.5*a2) */
    a3 = .0310907;
    b3 = 7.06042;
    c3 = 18.0578;
    d3 = -.325;

/* t4 = (1/(4/3)*pi)**(1/3) */
    t4 = .620350490899399531;

/* t5 = 0.5/(2**(1/3)-1) */
    t5 = 1.92366105093153617;

/* t6 = 2/(3*(2**(1/3)-1)) */
    t6 = 2.56488140124204822;

/* t7 = 2.25*(2**(1/3)-1) */
    t7 = .584822362263464735;

/* Spin stiffness interpolation constants */

    s1 = 7.12310891781811772;
    s2 = a1 * .5;
    s3 = -6.9917323507644313e-6;
    s4 = -.0053650918488835769;

/* Paramagnetic interpolation constants */

    p1 = 6.1519908197590798;
    p2 = a2 * .5;
    p3 = 9.6902277115443745e-4;
    p4 = .038783294878113009;

/* Ferromagnetic interpolation constants */

    f1 = 4.73092690956011364;
    f2 = a3 * .5;

/*      F3 = -0.244185082989490298d-02 *0.5d0 */
/*      F4 = -0.570212323620622186d-01 *0.5d0 */

/*  try nwchem values */

    f3 = .00224786709554261133;
    f4 = .0524913931697809227;

/* Interpolation intervals */

    inter1 = .99999999989999999;
    inter2 = -.99999999989999999;

/* open shell case */
    alpha_rho13__ = pow(ra, &c_b7);
    beta_rho13__ = pow(rb, &c_b7);
    srho = *ra + *rb;
    srho13 = pow(&srho, &c_b7);
    iv2 = t4 / srho13;
    iv = sqrt(iv2);

/* spin-stiffness */
    inv = 1. / (iv2 + b1 * iv + c1);
    i1 = log(iv2 * inv);
    i2 = log((iv - d1) * (iv - d1) * inv);
    i3 = atan(s1 / (iv * 2. + b1));
    ss1 = s2 * i1 + s3 * i2 + s4 * i3;
    ss2 = a1 * (1. / iv - iv * inv * (b1 / (iv - d1) + 1.));

/* paramagnetic */
    inv = 1. / (iv2 + b2 * iv + c2);
    i1 = log(iv2 * inv);
    i2 = log((iv - d2) * (iv - d2) * inv);
/* corrected b1->b2 ps Apr98 */
    i3 = atan(p1 / (iv * 2. + b2));
    pp1 = p2 * i1 + p3 * i2 + p4 * i3;
    pp2 = a2 * (1. / iv - iv * inv * (b2 / (iv - d2) + 1.));

/* ferromagnetic */
    inv = 1. / (iv2 + b3 * iv + c3);
    i1 = log(iv2 * inv);
    i2 = log((iv - d3) * (iv - d3) * inv);
    i3 = atan(f1 / (iv * 2. + b3));
    ff1 = f2 * i1 + f3 * i2 + f4 * i3;
    ff2 = a3 * (1. / iv - iv * inv * (b3 / (iv - d3) + 1.));

/* polarisation function */

    zeta = (*ra - *rb) / srho;
    zeta3 = zeta * zeta * zeta;
    zeta4 = zeta3 * zeta;
    if (zeta > inter1) {
  vwn1 = t5 * .51984209978974638;
  vwn2 = t6 * 1.25992104989487316476721060727823;
    } else if (zeta < inter2) {
  vwn1 = t5 * .51984209978974638;
  vwn2 = t6 * -1.25992104989487316476721060727823;
    } else {
  d__1 = zeta + 1.;
  d__2 = 1. - zeta;
  vwn1 = (pow(&d__1, &c_b14) + pow(&d__2, &c_b14) - 2.) * t5;
  d__1 = zeta + 1.;
  d__2 = 1. - zeta;
  vwn2 = (pow(&d__1, &c_b7) - pow(&d__2, &c_b7)) * t6;
    }
    ss1 *= t7;
    ss2 *= t7;
    tau = ff1 - pp1 - ss1;
    dtau = ff2 - pp2 - ss2;

    v = pp1 + vwn1 * (ss1 + tau * zeta4);
    *f = v * srho;

    t1 = v - iv * .166666666666666666666666666667 * (pp2 + vwn1 * (ss2 + dtau
      * zeta4));
    t2 = vwn2 * (ss1 + tau * zeta4) + vwn1 * 4. * tau * zeta3;
    *dfdra = t1 + t2 * (1. - zeta);
    *dfdrb = t1 - t2 * (zeta + 1.);

    return 0;
} /* c_uks_vwn5__ */

////***************************************************************************
//int rks_x_lda__ (integer *ideriv, integer *npt, doublereal *rhoa1,
//                  doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa,
//                  doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa,
//                  doublereal *v2sigmaaa2);
////***************************************************************************
//
////***************************************************************************
//int rks_c_vwn5__ (integer *ideriv, integer *npt, doublereal *rhoa1,
//                  doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa,
//                  doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa,
//                  doublereal *v2sigmaaa2);
////***************************************************************************

inline /* Subroutine */ int rks_c_vwn5__(integer *ideriv, integer *npt, doublereal *
  rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa,
  doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa,
  doublereal *v2sigmaaa2)
{
  // WSTHORNTON
  static doublereal c_b2 = .16666666666666666;
  static doublereal c_b3 = .33333333333333331;
   /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal), atan(
      doublereal);

    /* Local variables */
    static integer i__;
    static doublereal s1, t1, t2, t4, t6, t7, t10, t20, t11, t22, t13, t23,
      t16, t17, t25, t19, t26, t28, t29, t32, t33, t37, t38, t40, t41,
      t43, t51, t52, t27, t34, t39, t46, t47, t48, t49, t53, t55, t56,
      t58, t60, t63, t66, t67, t68, t69, t77, t79, t80, t81, t88, t92,
      t94, t102, t103, t105, t107, t125, t134, t138, rho;


/*     S.H. Vosko, L. Wilk, and M. Nusair */
/*     Accurate spin-dependent electron liquid correlation energies for */
/*     local spin density calculations: a critical analysis */
/*     Can. J. Phys. 58 (1980) 1200-1211 */


/*     CITATION: */

/*     Functionals were obtained from the Density Functional Repository */
/*     as developed and distributed by the Quantum Chemistry Group, */
/*     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD */
/*     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or */
/*     Paul Sherwood for further information. */

/*     COPYRIGHT: */

/*     Users may incorporate the source code into software packages and */
/*     redistribute the source code provided the source code is not */
/*     changed in anyway and is properly cited in any documentation or */
/*     publication related to its use. */

/*     ACKNOWLEDGEMENT: */

/*     The source code was generated using Maple 8 through a modified */
/*     version of the dfauto script published in: */

/*        R. Strange, F.R. Manby, P.J. Knowles */
/*        Automatic code generation in density functional theory */
/*        Comp. Phys. Comm. 136 (2001) 310-318. */

    /* Parameter adjustments */
    --v2sigmaaa2;
    --v2rhoasigmaaa;
    --v2rhoa2;
    --vsigmaaa;
    --vrhoa;
    --zk;
    --sigmaaa1;
    --rhoa1;

    /* Function Body */
    if (*ideriv == 0) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = 1 / rho;
    t2 = pow_dd(&t1, &c_b3);
    t4 = pow_dd(&t1, &c_b2);
    t7 = 1 / (t2 * .6203504908994 + t4 * 2.935818660072219 +
      12.9352);
    t10 = log(t2 * .6203504908994 * t7);
    t16 = atan(6.15199081975908 / (t4 * 1.575246635799487 +
      3.72744));
/* Computing 2nd power */
    d__1 = t4 * .7876233178997433 + .10498;
    t20 = d__1 * d__1;
    t22 = log(t20 * t7);
    zk[i__] = rho * (t10 * .0310907 + t16 * .03878329487811301 +
      t22 * 9.690227711544374e-4);
      } else {
/* rho */
    zk[i__] = 0.;
      }
/* rho */
  }
    } else if (*ideriv == 1) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = 1 / rho;
    t2 = pow_dd(&t1, &c_b3);
    t4 = pow_dd(&t1, &c_b2);
    t6 = t2 * .6203504908994 + t4 * 2.935818660072219 + 12.9352;
    t7 = 1 / t6;
    t10 = log(t2 * .6203504908994 * t7);
    t11 = t10 * .0310907;
    t13 = t4 * 1.575246635799487 + 3.72744;
    t16 = atan(6.15199081975908 / t13);
    t17 = t16 * .03878329487811301;
    t19 = t4 * .7876233178997433 + .10498;
/* Computing 2nd power */
    d__1 = t19;
    t20 = d__1 * d__1;
    t22 = log(t20 * t7);
    t23 = t22 * 9.690227711544374e-4;
    zk[i__] = rho * (t11 + t17 + t23);
/* Computing 2nd power */
    d__1 = t2;
    t25 = d__1 * d__1;
    t26 = 1 / t25;
/* Computing 2nd power */
    d__1 = rho;
    t28 = d__1 * d__1;
    t29 = 1 / t28;
/* Computing 2nd power */
    d__1 = t6;
    t32 = d__1 * d__1;
    t33 = 1 / t32;
/* Computing 2nd power */
    d__1 = t4;
    t37 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = t37;
    t38 = d__1 * d__1;
    t40 = 1 / t38 / t4;
    t41 = t40 * t29;
    t43 = t26 * -.2067834969664667 * t29 - t41 *
      .4893031100120365;
/* Computing 2nd power */
    d__1 = t13;
    t51 = d__1 * d__1;
    t52 = 1 / t51;
    vrhoa[i__] = t11 + t17 + t23 + rho * ((t26 *
      -.2067834969664667 * t7 * t29 - t2 * .6203504908994 *
      t33 * t43) * .05011795824473985 / t2 * t6 + t52 *
      .0626408570946439 * t40 * t29 / (t52 * 37.8469910464
      + 1.) + (t19 * -.2625411059665811 * t7 * t41 - t20 *
      1. * t33 * t43) * 9.690227711544374e-4 / t20 * t6);
      } else {
/* rho */
    zk[i__] = 0.;
    vrhoa[i__] = 0.;
      }
/* rho */
  }
    } else if (*ideriv == 2) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = 1 / rho;
    t2 = pow_dd(&t1, &c_b3);
    t4 = pow_dd(&t1, &c_b2);
    t6 = t2 * .6203504908994 + t4 * 2.935818660072219 + 12.9352;
    t7 = 1 / t6;
    t10 = log(t2 * .6203504908994 * t7);
    t11 = t10 * .0310907;
    t13 = t4 * 1.575246635799487 + 3.72744;
    t16 = atan(6.15199081975908 / t13);
    t17 = t16 * .03878329487811301;
    t19 = t4 * .7876233178997433 + .10498;
/* Computing 2nd power */
    d__1 = t19;
    t20 = d__1 * d__1;
    t22 = log(t20 * t7);
    t23 = t22 * 9.690227711544374e-4;
    zk[i__] = rho * (t11 + t17 + t23);
/* Computing 2nd power */
    d__1 = t2;
    t25 = d__1 * d__1;
    t26 = 1 / t25;
    t27 = t26 * t7;
/* Computing 2nd power */
    d__1 = rho;
    t28 = d__1 * d__1;
    t29 = 1 / t28;
/* Computing 2nd power */
    d__1 = t6;
    t32 = d__1 * d__1;
    t33 = 1 / t32;
    t34 = t2 * t33;
/* Computing 2nd power */
    d__1 = t4;
    t37 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = t37;
    t38 = d__1 * d__1;
    t39 = t38 * t4;
    t40 = 1 / t39;
    t41 = t40 * t29;
    t43 = t26 * -.2067834969664667 * t29 - t41 *
      .4893031100120365;
    t46 = t27 * -.2067834969664667 * t29 - t34 * .6203504908994 *
      t43;
    t47 = 1 / t2;
    t48 = t46 * t47;
    t49 = t48 * t6;
/* Computing 2nd power */
    d__1 = t13;
    t51 = d__1 * d__1;
    t52 = 1 / t51;
    t53 = t52 * t40;
    t55 = t52 * 37.8469910464 + 1.;
    t56 = 1 / t55;
    t58 = t53 * t29 * t56;
    t60 = t19 * t7;
    t63 = t20 * t33;
    t66 = t60 * -.2625411059665811 * t41 - t63 * 1. * t43;
    t67 = 1 / t20;
    t68 = t66 * t67;
    t69 = t68 * t6;
    vrhoa[i__] = t11 + t17 + t23 + rho * (t49 *
      .05011795824473985 + t58 * .0626408570946439 + t69 *
      9.690227711544374e-4);
    t77 = 1 / t25 / t1;
/* Computing 2nd power */
    d__1 = t28;
    t79 = d__1 * d__1;
    t80 = 1 / t79;
    t81 = t77 * t7 * t80;
    t88 = 1 / t28 / rho;
    t92 = 1 / t32 / t6;
/* Computing 2nd power */
    d__1 = t43;
    t94 = d__1 * d__1;
    t102 = 1 / t39 / t1;
    t103 = t102 * t80;
    t105 = t40 * t88;
    t107 = t77 * -.1378556646443111 * t80 + t26 *
      .4135669939329333 * t88 - t103 * .4077525916766971 +
      t105 * .978606220024073;
    t125 = t80 * t56;
/* Computing 2nd power */
    d__1 = t51;
    t134 = d__1 * d__1;
/* Computing 2nd power */
    d__1 = t55;
    t138 = d__1 * d__1;
    s1 = t49 * .2004718329789594 + t58 * .2505634283785756;
    v2rhoa2[i__] = s1 + t69 * .00387609108461775 + rho * 2. * ((
      t81 * -.1378556646443111 + t26 * .4135669939329333 *
      t33 * t29 * t43 + t27 * .4135669939329333 * t88 + t2 *
       1.2407009817988 * t92 * t94 - t34 * .6203504908994 *
      t107) * .05011795824473985 * t47 * t6 + t46 *
      .01670598608157995 / t2 / t1 * t6 * t29 + t48 *
      .05011795824473985 * t43 + .03289159980064473 / t51 /
      t13 * t77 * t125 + t52 * .05220071424553658 * t102 *
      t125 - t53 * .1252817141892878 * t88 * t56 -
      1.244848083156773 / t134 / t13 * t77 * t80 / t138 + (
      t81 * .03446391616107778 + t19 * .5250822119331622 *
      t33 * t41 * t43 - t60 * .2187842549721509 * t103 +
      t60 * .5250822119331622 * t105 + t20 * 2. * t92 * t94
      - t63 * 1. * t107) * 9.690227711544374e-4 * t67 * t6
      + t66 * 2.544083100456872e-4 / t20 / t19 * t6 * t40 *
      t29 + t68 * 9.690227711544374e-4 * t43);
      } else {
/* rho */
    zk[i__] = 0.;
    vrhoa[i__] = 0.;
    v2rhoa2[i__] = 0.;
      }
/* rho */
  }
    }
/* ideriv */
    return 0;
} /* rks_c_vwn5__ */

inline /* Subroutine */ int rks_x_lda__(integer *ideriv, integer *npt, doublereal *
  rhoa1, doublereal *sigmaaa1, doublereal *zk, doublereal *vrhoa,
  doublereal *vsigmaaa, doublereal *v2rhoa2, doublereal *v2rhoasigmaaa,
  doublereal *v2sigmaaa2)
{
  // WSTHORNTON
  static doublereal c_b2 = .33333333333333331;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__;
    static doublereal t1, t5, rho;


/*     P.A.M. Dirac */
/*     Proceedings of the Cambridge Philosophical Society, 26 (1930) 376 */


/*     CITATION: */

/*     Functionals were obtained from the Density Functional Repository */
/*     as developed and distributed by the Quantum Chemistry Group, */
/*     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD */
/*     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or */
/*     Paul Sherwood for further information. */

/*     COPYRIGHT: */

/*     Users may incorporate the source code into software packages and */
/*     redistribute the source code provided the source code is not */
/*     changed in anyway and is properly cited in any documentation or */
/*     publication related to its use. */

/*     ACKNOWLEDGEMENT: */

/*     The source code was generated using Maple 8 through a modified */
/*     version of the dfauto script published in: */

/*        R. Strange, F.R. Manby, P.J. Knowles */
/*        Automatic code generation in density functional theory */
/*        Comp. Phys. Comm. 136 (2001) 310-318. */

    /* Parameter adjustments */
    --v2sigmaaa2;
    --v2rhoasigmaaa;
    --v2rhoa2;
    --vsigmaaa;
    --vrhoa;
    --zk;
    --sigmaaa1;
    --rhoa1;

    /* Function Body */
    if (*ideriv == 0) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = pow_dd(&rho, &c_b2);
    zk[i__] = t1 * -.7385587663820224 * rho;
      } else {
/* rho */
    zk[i__] = 0.;
      }
/* rho */
  }
    } else if (*ideriv == 1) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = pow_dd(&rho, &c_b2);
    zk[i__] = t1 * -.7385587663820224 * rho;
    vrhoa[i__] = t1 * -.9847450218426965;
      } else {
/* rho */
    zk[i__] = 0.;
    vrhoa[i__] = 0.;
      }
/* rho */
  }
    } else if (*ideriv == 2) {
  i__1 = *npt;
  for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
      d__1 = 0., d__2 = rhoa1[i__];
      rho = max(d__1,d__2);
      if (rho > 1e-20) {
    t1 = pow_dd(&rho, &c_b2);
    zk[i__] = t1 * -.7385587663820224 * rho;
    vrhoa[i__] = t1 * -.9847450218426965;
/* Computing 2nd power */
    d__1 = t1;
    t5 = d__1 * d__1;
    v2rhoa2[i__] = -.6564966812284644 / t5;
      } else {
/* rho */
    zk[i__] = 0.;
    vrhoa[i__] = 0.;
    v2rhoa2[i__] = 0.;
      }
/* rho */
  }
    }
/* ideriv */
    return 0;
} /* rks_x_lda__ */

const double THRESH_RHO = 1e-8;
const double THRESH_GRHO = 1e-20;

inline void wst_munge_grho(int npoint, double *rho, double *grho) {
    for (int i=0; i<npoint; i++) {
        if (rho[i]<THRESH_RHO) rho[i] = THRESH_RHO;
        if ((rho[i] <=THRESH_RHO) ||
            (grho[i] < THRESH_GRHO)) grho[i] = THRESH_GRHO;
    }
}

inline void wst_munge_rho(int npoint, double *rho) {
    for (int i=0; i<npoint; i++) {
        if (rho[i]<THRESH_RHO) rho[i] = THRESH_RHO;
    }
}

//***************************************************************************
inline void xc_rks_generic_lda(Tensor<double> rho_alpha,           ///< Alpha-spin density at each grid point
                          Tensor<double> f,                         ///< Value of functional at each grid point
                          Tensor<double> df_drho)                   ///< Derivative of functional w.r.t. rho_alpha
  {
    MADNESS_ASSERT(rho_alpha.iscontiguous());
    MADNESS_ASSERT(f.iscontiguous());
    MADNESS_ASSERT(df_drho.iscontiguous());

    rho_alpha = rho_alpha.flat();
    f = f.flat();
    df_drho = df_drho.flat();

      integer ideriv = 2;
      integer npt = rho_alpha.dim(0);

      Tensor<double> gamma_alpha(npt);
      Tensor<double> tf(npt);
      Tensor<double> tdf_drho(npt);
      Tensor<double> tdf_dgamma(npt);
      Tensor<double> td2f_drho2(npt);
      Tensor<double> td2f_drhodgamma(npt);
      Tensor<double> td2f_dgamma2(npt);

      wst_munge_rho(npt, rho_alpha.ptr());

      f.fill(0.0);
      df_drho.fill(0.0);

      int returnvalue = ::rks_x_lda__(&ideriv, &npt, rho_alpha.ptr(), gamma_alpha.ptr(),
               tf.ptr(),
               tdf_drho.ptr(), tdf_dgamma.ptr(),
               td2f_drho2.ptr(), td2f_drhodgamma.ptr(), td2f_dgamma2.ptr());

      f.gaxpy(1.0, tf, 1.0);
      df_drho.gaxpy(1.0, tdf_drho, 1.0);

      tf.fill(0.0);
      tdf_drho.fill(0.0);

      returnvalue = ::rks_c_vwn5__(&ideriv, &npt, rho_alpha.ptr(), gamma_alpha.ptr(),
                tf.ptr(),
                tdf_drho.ptr(), tdf_dgamma.ptr(),
                td2f_drho2.ptr(), td2f_drhodgamma.ptr(), td2f_dgamma2.ptr());

      f.gaxpy(1.0, tf, 1.0);
      df_drho.gaxpy(1.0, tdf_drho, 1.0);

  }
  //***************************************************************************

  //***************************************************************************
  template <int NDIM>
 inline void dft_xc_lda_V(const Key<NDIM>& key, Tensor<double>& t)
  {
    Tensor<double> enefunc = copy(t);
    Tensor<double> V = copy(t);
    ::xc_rks_generic_lda(t, enefunc, V);
    t(___) = V(___);
  }
  //***************************************************************************

  //***************************************************************************
  template <int NDIM>
 inline void dft_xc_lda_ene(const Key<NDIM>& key, Tensor<double>& t)
  {
    Tensor<double> V = copy(t);
    Tensor<double> enefunc = copy(t);
    ::xc_rks_generic_lda(t, enefunc, V);
    t(___) = enefunc(___);
  }
  //***************************************************************************

//  //***************************************************************************
//  static double munge(double r) {
//      if (r < 1e-12) r = 1e-12;
//      return r;
//  }
//  //***************************************************************************
//
//  //***************************************************************************
//  static void ldaop(const Key<3>& key, Tensor<double>& t) {
//      UNARY_OPTIMIZED_ITERATOR(double, t, double r=munge(2.0* *_p0); double q; double dq1; double dq2;x_rks_s__(&r, &q, &dq1);c_rks_vwn5__(&r, &q, &dq2); *_p0 = dq1+dq2);
//  }
//  //***************************************************************************
//
//  //***************************************************************************
//  static void ldaeop(const Key<3>& key, Tensor<double>& t) {
//      UNARY_OPTIMIZED_ITERATOR(double, t, double r=munge(2.0* *_p0); double q1; double q2; double dq;x_rks_s__(&r, &q1, &dq);c_rks_vwn5__(&r, &q2, &dq); *_p0 = q1+q2);
//  }
//  //***************************************************************************



#endif /* LDA_H_ */
