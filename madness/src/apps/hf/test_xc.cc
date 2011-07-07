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
  
  $Id$
*/
#include <mra/mra.h>
#include <iostream>

#include "libxc.h"

using std::cout;
using std::endl;

using namespace madness;

const double PI = 3.1415926535897932384;

typedef Vector<double,3> coordT;

//static inline double pow(const double* a, const double* b) {
//    return pow(*a, *b);
//}

//static double c_b2 = .333333333333333333333333333333333;
//static double c_b7 = .333333333333333333333333333333;
//static double c_b8 = .5;
//static double c_b14 = 1.333333333333333333333333333333;

//inline /* Subroutine */ int x_rks_s__(const double *r__, double *f, double *
//  dfdra)
//{
//
//    /* Local variables */
//    static double ra13;
//
//
///*     This subroutine evaluates the spin polarised exchange functional */
///*     in the Local Density Approximation [1], and the corresponding */
///*     potential. Often this functional is referred to as the Dirac */
///*     functional [2] or Slater functional. */
//
///*     [1] F. Bloch, Zeitschrift fuer Physik, Vol. 57 (1929) 545. */
//
///*     [2] P.A.M. Dirac, Proceedings of the Cambridge Philosophical */
///*         Society, Vol. 26 (1930) 376. */
//
///*     Parameters: */
//
///*     r     the total electron density */
///*     f     On return the functional value */
///*     dfdra On return the derivative of f with respect to alpha electron */
///*           density */
//
//
///*     Ax = -3/4*(6/pi)**(1/3) */
///*     Bx = -(6/pi)**(1/3) */
///*     C  = (1/2)**(1/3) */
//
//
//
//
//    ra13 = pow(r__, &c_b2) * .793700525984099737375852819636154;
//    *f = *r__ * -.930525736349100025002010218071667 * ra13;
//    *dfdra = ra13 * -1.24070098179880003333601362409556;
//
//    return 0;
//} /* x_rks_s__ */
//
///* ----------------------------------------------------------------------- */
//inline /* Subroutine */ int x_uks_s__(double *ra, double *rb, double *f,
//  double *dfdra, double *dfdrb)
//{
//    /* Local variables */
//    static double ra13, rb13;
//
//
///*     This subroutine evaluates the spin polarised exchange functional */
///*     in the Local Density Approximation [1], and the corresponding */
///*     potential. Often this functional is referred to as the Dirac */
///*     functional [2] or Slater functional. */
//
///*     [1] F. Bloch, Zeitschrift fuer Physik, Vol. 57 (1929) 545. */
//
///*     [2] P.A.M. Dirac, Proceedings of the Cambridge Philosophical */
///*         Society, Vol. 26 (1930) 376. */
//
///*     Parameters: */
//
///*     ra    the alpha electron density */
///*     rb    the beta  electron density */
///*     f     On return the functional value */
///*     dfdra On return the derivative of f with respect to ra */
///*     dfdrb On return the derivative of f with respect to rb */
//
//
///*     Ax = -3/4*(6/pi)**(1/3) */
///*     Bx = -(6/pi)**(1/3) */
//
//
//
//
//    ra13 = pow(ra, &c_b2);
//    rb13 = pow(rb, &c_b2);
//    *f = (*ra * ra13 + *rb * rb13) * -.930525736349100025002010218071667;
//    *dfdra = ra13 * -1.24070098179880003333601362409556;
//    *dfdrb = rb13 * -1.24070098179880003333601362409556;
//
//    return 0;
//} /* x_uks_s__ */
//
//inline /* Subroutine */ int c_rks_vwn5__(const double *r__, double *f, double *
//  dfdra)
//{
//    /* Local variables */
//    static double a2, b2, c2, d2, i1, i2, i3, p1, p2, p3, p4, t4, t5, t6,
//      t7, iv, alpha_rho13__, iv2, pp1, pp2, inv, srho, srho13;
//
//
///*     This subroutine evaluates the Vosko, Wilk and Nusair correlation */
///*     functional number 5 [1] for the closed shell case, with the */
///*     parametrisation as given in table 5. */
//
///*     The original code was obtained from Dr. Phillip Young, */
///*     with corrections from Dr. Paul Sherwood. */
//
///*     [1] S.H. Vosko, L. Wilk, and M. Nusair */
///*         "Accurate spin-dependent electron liquid correlation energies */
///*          for local spin density calculations: a critical analysis", */
///*         Can.J.Phys, Vol. 58 (1980) 1200-1211. */
//
///*     Parameters: */
//
///*     r      the total electron density */
///*     f      On return the functional value */
///*     dfdra  On return the derivative of f with respect to the alpha */
///*            electron density */
//
//
//
//
///* VWN interpolation parameters */
//
///* paramagnetic */
//    a2 = .0621814;
//    b2 = 3.72744;
//    c2 = 12.9352;
//    d2 = -.10498;
//
///* t4 = (1/(4/3)*pi)**(1/3) */
//    t4 = .620350490899399531;
//
///* t5 = 0.5/(2**(1/3)-1) */
//    t5 = 1.92366105093153617;
//
///* t6 = 2/(3*(2**(1/3)-1)) */
//    t6 = 2.56488140124204822;
//
///* t7 = 2.25*(2**(1/3)-1) */
//    t7 = .584822362263464735;
//
///* Paramagnetic interpolation constants */
//
//    p1 = 6.1519908197590798;
//    p2 = a2 * .5;
//    p3 = 9.6902277115443745e-4;
//    p4 = .038783294878113009;
//
///* closed shell case */
//    srho = *r__;
//    srho13 = pow(&srho, &c_b7);
//    alpha_rho13__ = pow(&c_b8, &c_b7) * srho;
//    iv2 = t4 / srho13;
//    iv = sqrt(iv2);
//
///* paramagnetic */
//    inv = 1. / (iv2 + b2 * iv + c2);
//    i1 = log(iv2 * inv);
//    i2 = log((iv - d2) * (iv - d2) * inv);
///* corrected b1->b2 ps Apr98 */
//    i3 = atan(p1 / (iv * 2. + b2));
//    pp1 = p2 * i1 + p3 * i2 + p4 * i3;
//    pp2 = a2 * (1. / iv - iv * inv * (b2 / (iv - d2) + 1.));
//
//    *f = pp1 * srho;
//    *dfdra = pp1 - iv * .166666666666666666666666666666 * pp2;
//
//    return 0;
//} /* c_rks_vwn5__ */
//
///* ----------------------------------------------------------------------- */
//inline /* Subroutine */ int c_uks_vwn5__(double *ra, double *rb, double *
//  f, double *dfdra, double *dfdrb)
//{
//    /* System generated locals */
//    double d__1, d__2;
//
//    /* Local variables */
//    static double v, beta_rho13__, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3,
//       c3, d3, f1, f2, f3, p1, p2, p3, s1, t1, t2, s2, t4, t5, t6, t7,
//      s3, s4, p4, f4, i1, i2, i3, iv, alpha_rho13__, ff1, ff2, iv2, pp1,
//       pp2, ss1, ss2, tau, inv, vwn1, vwn2, dtau, zeta, srho, zeta3,
//      zeta4, srho13, inter1, inter2;
//
//
///*     This subroutine evaluates the Vosko, Wilk and Nusair correlation */
///*     functional number 5 [1], with the parametrisation as given in */
///*     table 5. */
//
///*     The original code was obtained from Dr. Phillip Young, */
///*     with corrections from Dr. Paul Sherwood. */
//
///*     [1] S.H. Vosko, L. Wilk, and M. Nusair */
///*         "Accurate spin-dependent electron liquid correlation energies */
///*          for local spin density calculations: a critical analysis", */
///*         Can.J.Phys, Vol. 58 (1980) 1200-1211. */
//
///*     Parameters: */
//
///*     ra     the alpha-electron density */
///*     rb     the beta-electron density */
///*     f      On return the functional value */
///*     dfdra  On return the derivative of f with respect to ra */
///*     dfdrb  On return the derivative of f with respect to rb */
//
//
//
///*     tn13 = 2**(1/3) */
///*     tn43 = 2**(4/3) */
//
///* VWN interpolation parameters */
//
///* spin stiffness */
//    a1 = -.0337737278807791058;
//    b1 = 1.13107;
//    c1 = 13.0045;
//    d1 = -.0047584;
///* paramagnetic */
//    a2 = .0621814;
//    b2 = 3.72744;
//    c2 = 12.9352;
//    d2 = -.10498;
///* ferromagnetic */
///* try cadpac/nwchem value (.5*a2) */
//    a3 = .0310907;
//    b3 = 7.06042;
//    c3 = 18.0578;
//    d3 = -.325;
//
///* t4 = (1/(4/3)*pi)**(1/3) */
//    t4 = .620350490899399531;
//
///* t5 = 0.5/(2**(1/3)-1) */
//    t5 = 1.92366105093153617;
//
///* t6 = 2/(3*(2**(1/3)-1)) */
//    t6 = 2.56488140124204822;
//
///* t7 = 2.25*(2**(1/3)-1) */
//    t7 = .584822362263464735;
//
///* Spin stiffness interpolation constants */
//
//    s1 = 7.12310891781811772;
//    s2 = a1 * .5;
//    s3 = -6.9917323507644313e-6;
//    s4 = -.0053650918488835769;
//
///* Paramagnetic interpolation constants */
//
//    p1 = 6.1519908197590798;
//    p2 = a2 * .5;
//    p3 = 9.6902277115443745e-4;
//    p4 = .038783294878113009;
//
///* Ferromagnetic interpolation constants */
//
//    f1 = 4.73092690956011364;
//    f2 = a3 * .5;
//
///*      F3 = -0.244185082989490298d-02 *0.5d0 */
///*      F4 = -0.570212323620622186d-01 *0.5d0 */
//
///*  try nwchem values */
//
//    f3 = .00224786709554261133;
//    f4 = .0524913931697809227;
//
///* Interpolation intervals */
//
//    inter1 = .99999999989999999;
//    inter2 = -.99999999989999999;
//
///* open shell case */
//    alpha_rho13__ = pow(ra, &c_b7);
//    beta_rho13__ = pow(rb, &c_b7);
//    srho = *ra + *rb;
//    srho13 = pow(&srho, &c_b7);
//    iv2 = t4 / srho13;
//    iv = sqrt(iv2);
//
///* spin-stiffness */
//    inv = 1. / (iv2 + b1 * iv + c1);
//    i1 = log(iv2 * inv);
//    i2 = log((iv - d1) * (iv - d1) * inv);
//    i3 = atan(s1 / (iv * 2. + b1));
//    ss1 = s2 * i1 + s3 * i2 + s4 * i3;
//    ss2 = a1 * (1. / iv - iv * inv * (b1 / (iv - d1) + 1.));
//
///* paramagnetic */
//    inv = 1. / (iv2 + b2 * iv + c2);
//    i1 = log(iv2 * inv);
//    i2 = log((iv - d2) * (iv - d2) * inv);
///* corrected b1->b2 ps Apr98 */
//    i3 = atan(p1 / (iv * 2. + b2));
//    pp1 = p2 * i1 + p3 * i2 + p4 * i3;
//    pp2 = a2 * (1. / iv - iv * inv * (b2 / (iv - d2) + 1.));
//
///* ferromagnetic */
//    inv = 1. / (iv2 + b3 * iv + c3);
//    i1 = log(iv2 * inv);
//    i2 = log((iv - d3) * (iv - d3) * inv);
//    i3 = atan(f1 / (iv * 2. + b3));
//    ff1 = f2 * i1 + f3 * i2 + f4 * i3;
//    ff2 = a3 * (1. / iv - iv * inv * (b3 / (iv - d3) + 1.));
//
///* polarisation function */
//
//    zeta = (*ra - *rb) / srho;
//    zeta3 = zeta * zeta * zeta;
//    zeta4 = zeta3 * zeta;
//    if (zeta > inter1) {
//  vwn1 = t5 * .51984209978974638;
//  vwn2 = t6 * 1.25992104989487316476721060727823;
//    } else if (zeta < inter2) {
//  vwn1 = t5 * .51984209978974638;
//  vwn2 = t6 * -1.25992104989487316476721060727823;
//    } else {
//  d__1 = zeta + 1.;
//  d__2 = 1. - zeta;
//  vwn1 = (pow(&d__1, &c_b14) + pow(&d__2, &c_b14) - 2.) * t5;
//  d__1 = zeta + 1.;
//  d__2 = 1. - zeta;
//  vwn2 = (pow(&d__1, &c_b7) - pow(&d__2, &c_b7)) * t6;
//    }
//    ss1 *= t7;
//    ss2 *= t7;
//    tau = ff1 - pp1 - ss1;
//    dtau = ff2 - pp2 - ss2;
//
//    v = pp1 + vwn1 * (ss1 + tau * zeta4);
//    *f = v * srho;
//
//    t1 = v - iv * .166666666666666666666666666667 * (pp2 + vwn1 * (ss2 + dtau
//      * zeta4));
//    t2 = vwn2 * (ss1 + tau * zeta4) + vwn1 * 4. * tau * zeta3;
//    *dfdra = t1 + t2 * (1. - zeta);
//    *dfdrb = t1 - t2 * (zeta + 1.);
//
//    return 0;
//} /* c_uks_vwn5__ */

//*****************************************************************************
static double gauss_func3(const coordT& r)
{
  const double x=r[0], y=r[1], z=r[2];
  return exp(-3.0*(x*x + y*y + z*z));
}
//*****************************************************************************

//*****************************************************************************
static double gauss_func1(const coordT& r)
{
  const double x=r[0], y=r[1], z=r[2];
  return exp(-1.0*(x*x + y*y + z*z));
}
//*****************************************************************************

//*****************************************************************************
void test_xc1(World& world)
{
  XC(lda_type) xc_c_func;
  XC(lda_type) xc_x_func;
  XC(lda_type) xc_c_func2;
  XC(lda_type) xc_x_func2;

  xc_lda_init(&xc_c_func, XC_LDA_C_VWN,XC_UNPOLARIZED);
  xc_lda_x_init(&xc_x_func, XC_UNPOLARIZED, 3, 0);
  xc_lda_init(&xc_c_func2, XC_LDA_C_VWN,XC_POLARIZED);
  xc_lda_x_init(&xc_x_func2, XC_POLARIZED, 3, 0);

  for (int k = -10; k < 10; k++)
  {
    double r = pow(3.5, k);
    double tmpr[2];
    tmpr[0] = 0.5 * r;
    tmpr[1] = 0.5 * r;
    double libxc_v1 = 0.0, libxc_v2 = 0.0, libxc_e1 = 0.0, libxc_e2 = 0.0;
    double da_v1 = 0.0, da_v2 = 0.0, da_e1 = 0.0, da_e2 = 0.0;
    double libxc_v1_sp[2], libxc_v2_sp[2], libxc_e1_sp[2], libxc_e2_sp[2];
    x_rks_s__(&r, &da_e1, &da_v1);
    c_rks_vwn5__(&r, &da_e2, &da_v2);
    xc_lda_vxc(&xc_x_func, &r, &libxc_e1, &libxc_v1);
    xc_lda_vxc(&xc_c_func, &r, &libxc_e2, &libxc_v2);
    xc_lda_vxc(&xc_x_func2, &tmpr[0], &libxc_e1_sp[0], &libxc_v1_sp[0]);
    xc_lda_vxc(&xc_c_func2, &tmpr[0], &libxc_e2_sp[0], &libxc_v2_sp[0]);

    printf("%.5e     %8e     %8e     %8e\n", r, libxc_e1+libxc_e2, da_e1+da_e2, libxc_e1_sp[0]+libxc_e2_sp[0]);
  }
}
//*****************************************************************************

////***************************************************************************
//static double munge(double r) {
//  if (r < 1e-12) r = 1e-12;
//  return r;
//}
////***************************************************************************

//***************************************************************************
static void da_ldaop(const Key<3>& key, Tensor<double>& t) {
  UNARY_OPTIMIZED_ITERATOR(double, t, double r=munge(2.0* *_p0); double q; double dq1; double dq2;x_rks_s__(&r, &q, &dq1);c_rks_vwn5__(&r, &q, &dq2); *_p0 = dq1+dq2);
}
//***************************************************************************

////***************************************************************************
//static void libxc_ldaop(const Key<3>& key, Tensor<double>& t) {
//  XC(lda_type) xc_c_func;
//  XC(lda_type) xc_x_func;
//  xc_lda_init(&xc_c_func, XC_LDA_C_VWN,XC_UNPOLARIZED);
//  xc_lda_x_init(&xc_x_func, XC_UNPOLARIZED, 3, 0);
//  UNARY_OPTIMIZED_ITERATOR(double, t, double r=munge(2.0* *_p0); double q; double dq1; double dq2;xc_lda_vxc(&xc_x_func, &r, &q, &dq1);xc_lda_vxc(&xc_c_func, &r, &q, &dq2); *_p0 = dq1+dq2);
//}
////***************************************************************************

//*****************************************************************************
void multiply_op(const Key<3>& key, Tensor<double> tcube,
                 Tensor<double> lcube,
                 Tensor<double> rcube)
{
  TERNARY_OPTIMIZED_ITERATOR(double, tcube, double, lcube, double, rcube, *_p0 = *_p1 * *_p2;);
}
//*****************************************************************************

//*****************************************************************************
void add_op(const Key<3>& key, Tensor<double> tcube,
                 Tensor<double> lcube,
                 Tensor<double> rcube)
{
  TERNARY_OPTIMIZED_ITERATOR(double, tcube, double, lcube, double, rcube, *_p0 = *_p1 + *_p2;);
}
//*****************************************************************************

//*****************************************************************************
void test_xc2(World& world)
{
  // Function defaults
  int funck = 8;
  double thresh = 1e-6;
  double L = 10.0;
  FunctionDefaults<3>::set_k(funck);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);

  Function<double,3> f = FunctionFactory<double,3>(world).f(gauss_func3);
  Function<double,3> f2 = FunctionFactory<double,3>(world).f(gauss_func1);
  Function<double,3> fda = copy(f);
  Function<double,3> flibxc = copy(f);
  fda.unaryop(&da_ldaop);
  flibxc.unaryop(&libxc_ldaop);
  f.scale(0.5);
  Function<double,3> flibxc_sp = binary_op(f, f, &::libxc_ldaop_sp);
  Function<double,3> prod = binary_op(f, f2, &::multiply_op);
  Function<double,3> sum = binary_op(f, f2, &::add_op);

  if (world.rank() == 0) printf("\n");
  double bstep = L / 100.0;
  fda.reconstruct();
  flibxc.reconstruct();
  for (int i=0; i<101; i++)
  {
    coordT p(-L/2 + i*bstep);
    if (world.rank() == 0) printf("%5.2f%18.8e%15.8f%15.8f%15.8f\n", p[0], f(p), fda(p), flibxc(p), flibxc_sp(p));
//    if (world.rank() == 0) printf("%.2f\t\t%.8e\t%.8f\t%.8f\t%.8f\n", p[0], f(p), fda(p), flibxc(p), flibxc_sp(p));
//    if (world.rank() == 0) printf("%.2f\t\t%.8f\t%.8f\n", p[0], f(p)+f2(p), sum(p));
  }
  if (world.rank() == 0) printf("\n");

}
//*****************************************************************************

#define TO_STRING(s) TO_STRING2(s)
#define TO_STRING2(s) #s

//*****************************************************************************
int main(int argc, char** argv)
{
  initialize(argc, argv);
  World world(MPI::COMM_WORLD);
  if (world.rank() == 0)
  {
    print("");
    print("--------------------------------------------");
    print("   MADNESS", " multiresolution testsuite");
    print("--------------------------------------------");
    print("");
    print("   number of processors ...", world.size());
    print("    processor frequency ...", cpu_frequency());
    print("            host system ...", TO_STRING(HOST_SYSTEM));
    print("             byte order ...", TO_STRING(MADNESS_BYTE_ORDER));
    print("          configured by ...", MADNESS_CONFIGURATION_USER);
    print("          configured on ...", MADNESS_CONFIGURATION_HOST);
    print("          configured at ...", MADNESS_CONFIGURATION_DATE);
    print("                    CXX ...", MADNESS_CONFIGURATION_CXX);
    print("               CXXFLAGS ...", MADNESS_CONFIGURATION_CXXFLAGS);
#ifdef WORLD_WATCHDOG
    print("               watchdog ...", WATCHDOG_BARK_INTERVAL,
        WATCHDOG_TIMEOUT);
#endif
#ifdef OPTERON_TUNE
    print("             tuning for ...", "opteron");
#elif defined(CORE_DUO_TUNE)
    print("             tuning for ...", "core duo");
#else
    print("             tuning for ...", "core2");
#endif
#ifdef BOUNDS_CHECKING
    print(" tensor bounds checking ...", "enabled");
#endif
#ifdef TENSOR_INSTANCE_COUNT
    print("  tensor instance count ...", "enabled");
#endif
    print(" ");
  }

  try
  {
    printf("WSTHORNTON: Starting up the world ... \n");

    startup(world,argc,argv);
    if (world.rank() == 0) print("Initial tensor instance count", BaseTensor::get_instance_count());
    test_xc2(world);
  }
  catch (const MPI::Exception& e)
  {
    //        print(e);
    error("caught an MPI exception");
  }
  catch (const madness::MadnessException& e)
  {
    print(e);
    error("caught a MADNESS exception");
  }
  catch (const madness::TensorException& e)
  {
    print(e);
    error("caught a Tensor exception");
  }
  catch (const char* s)
  {
    print(s);
    error("caught a string exception");
  }
  catch (const std::string& s)
  {
    print(s);
    error("caught a string (class) exception");
  }
  catch (const std::exception& e)
  {
    print(e.what());
    error("caught an STL exception");
  }
  catch (...)
  {
    error("caught unhandled exception");
  }

  if (world.rank() == 0)
    print("entering final fence");
  world.gop.fence();
  if (world.rank() == 0)
    print("done with final fence");
  if (world.rank() == 0)
    print("Final tensor instance count", BaseTensor::get_instance_count());

  finalize();
  return 0;
}
//*****************************************************************************

