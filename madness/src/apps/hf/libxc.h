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
  
  $Id: libxc.h 1822 2010-02-11 20:10:27Z wsttiger@gmail.com $
*/
/*
 * libxc.h
 *
 *  Created on: Nov 23, 2008
 *      Author: wsttiger
 */

#ifndef LIBXC_H_
#define LIBXC_H_

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <world/world.h>
//#include "xc.h"
#include "lda.h"

using namespace madness;

//***************************************************************************
static double munge(double r) {
  if (r < 1e-15) r = 2e-15;
  return r;
}
//***************************************************************************

//***************************************************************************
template <typename T>
inline static void ldaop(const Key<3>& key, Tensor<T>& t) {
    UNARY_OPTIMIZED_ITERATOR(T, t, double r=munge(2.0* *_p0); double q; double dq1; double dq2;x_rks_s__(&r, &q, &dq1);c_rks_vwn5__(&r, &q, &dq2); *_p0 = dq1+dq2);
}
//***************************************************************************

//***************************************************************************
template <typename T>
inline static void ldaeop(const Key<3>& key, Tensor<T>& t) {
    UNARY_OPTIMIZED_ITERATOR(T, t, double r=munge(2.0* *_p0); double q1; double q2; double dq;x_rks_s__(&r, &q1, &dq);c_rks_vwn5__(&r, &q2, &dq); *_p0 = q1+q2);
}
//***************************************************************************

////***************************************************************************
//template <typename T>
//inline static void libxc_ldaop(const Key<3>& key, Tensor<T>& t) {
//  XC(lda_type) xc_c_func;
//  XC(lda_type) xc_x_func;
//  xc_lda_init(&xc_c_func, XC_LDA_C_VWN,XC_UNPOLARIZED);
//  xc_lda_x_init(&xc_x_func, XC_UNPOLARIZED, 3, 0);
//  UNARY_OPTIMIZED_ITERATOR(T, t, double r=munge(2.0* *_p0); double q; double dq1; double dq2;
//                           xc_lda_vxc(&xc_x_func, &r, &q, &dq1); xc_lda_vxc(&xc_c_func, &r, &q, &dq2);
//                           *_p0 = dq1+dq2);
//}
////***************************************************************************

////***************************************************************************
//template <typename T>
//inline static void libxc_ldaop_sp(const Key<3>& key, Tensor<T>& t, Tensor<T>& a, Tensor<T>& b)
//{
//  XC(lda_type) xc_c_func;
//  XC(lda_type) xc_x_func;
//  xc_lda_init(&xc_c_func, XC_LDA_C_VWN,XC_POLARIZED);
//  xc_lda_x_init(&xc_x_func, XC_POLARIZED, 3, 0);
//  TERNARY_OPTIMIZED_ITERATOR(T, t, T, a, T, b, double r[2]; r[0] = munge(*_p1);
//                             r[1] = munge(*_p2); double q[2]; double dq1[2]; double dq2[2];
//                             xc_lda_vxc(&xc_x_func, &r[0], &q[0], &dq1[0]); xc_lda_vxc(&xc_c_func, &r[0], &q[0], &dq2[0]);
//                             *_p0 = dq1[0]+dq2[0]);
//}
////***************************************************************************

////***************************************************************************
//template <typename T>
//inline static void libxc_ldaeop_sp(const Key<3>& key, Tensor<T>& t, Tensor<T>& a, Tensor<T>& b)
//{
//  XC(lda_type) xc_c_func;
//  XC(lda_type) xc_x_func;
//  xc_lda_init(&xc_c_func, XC_LDA_C_VWN,XC_POLARIZED);
//  xc_lda_x_init(&xc_x_func, XC_POLARIZED, 3, 0);
//  TERNARY_OPTIMIZED_ITERATOR(T, t, T, a, T, b, double r[2]; r[0] = munge(*_p1);
//                             r[1] = munge(*_p2); double q1[2]; double q2[2]; double dq[2];
//                             xc_lda_vxc(&xc_x_func, &r[0], &q1[0], &dq[0]); xc_lda_vxc(&xc_c_func, &r[0], &q2[0], &dq[0]);
//                             *_p0 = q1[0]+q2[0]);
//}
////***************************************************************************

////***************************************************************************
//inline static void libxc_ldaeop_sp(const Key<3>& key, Tensor<double>& t) {
//  XC(lda_type) xc_c_func;
//  XC(lda_type) xc_x_func;
//  xc_lda_init(&xc_c_func, XC_LDA_C_VWN,XC_UNPOLARIZED);
//  xc_lda_x_init(&xc_x_func, XC_UNPOLARIZED, 3, 0);
//  UNARY_OPTIMIZED_ITERATOR(double, t, double r=munge(2.0* *_p0); double q1; double q2; double dq; xc_lda_vxc(&xc_x_func, &r, &q1, &dq); xc_lda_vxc(&xc_c_func, &r, &q2, &dq); *_p0 = q1+q2);
//}
////***************************************************************************

//const double THRESH_RHO = 1e-8;
//const double THRESH_GRHO = 1e-20;
//
////***************************************************************************
//inline void wst_munge_grho(int npoint, double *rho, double *grho) {
//    for (int i=0; i<npoint; i++) {
//        if (rho[i]<THRESH_RHO) rho[i] = THRESH_RHO;
//        if ((rho[i] <=THRESH_RHO) ||
//            (grho[i] < THRESH_GRHO)) grho[i] = THRESH_GRHO;
//    }
//}
////***************************************************************************
//
////***************************************************************************
//inline void wst_munge_rho(int npoint, double *rho) {
//    for (int i=0; i<npoint; i++) {
//        if (rho[i]<THRESH_RHO) rho[i] = THRESH_RHO;
//    }
//}
////***************************************************************************
//
////***************************************************************************
//inline void xc_generic_lda(Tensor<double> rho_alpha,           ///< Alpha-spin density at each grid point
//                          Tensor<double> f,                         ///< Value of functional at each grid point
//                          Tensor<double> df_drho,                   ///< Derivative of functional w.r.t. rho_alpha
//                          bool spinpol)
//    {
//    MADNESS_ASSERT(rho_alpha.iscontiguous());
//    MADNESS_ASSERT(f.iscontiguous());
//    MADNESS_ASSERT(df_drho.iscontiguous());
//
//    rho_alpha = rho_alpha.flat();
//    f = f.flat();
//    df_drho = df_drho.flat();
//
//    XC(lda_type) xc_c_func;
//    XC(lda_type) xc_x_func;
//
//    int npt = rho_alpha.dim(0);
//
//    Tensor<double> tf(npt);
//    Tensor<double> tdf_drho(npt);
//    double* rhoptr = rho_alpha.ptr();
//    double* tfptr = tf.ptr();
//    double* tdf_drhoptr = tdf_drho.ptr();
//
//    tf.fill(0.0);
//    tdf_drho.fill(0.0);
//    f.fill(0.0);
//    df_drho.fill(0.0);
//
//    wst_munge_rho(npt, rhoptr);
//
//    xc_lda_init(&xc_c_func, XC_LDA_C_VWN,XC_UNPOLARIZED);
//    for (int i = 0; i < npt; i++)
//    {
//      xc_lda_vxc(&xc_c_func, &rhoptr[i], &tfptr[i], &tdf_drhoptr[i]);
//    }
//
//    f.gaxpy(1.0, tf, 1.0);
//    df_drho.gaxpy(1.0, tdf_drho, 1.0);
//
//    tf.fill(0.0);
//    tdf_drho.fill(0.0);
//
//    xc_lda_x_init(&xc_x_func, XC_UNPOLARIZED, 3, 0);
//    for (int i = 0; i < npt; i++)
//    {
//      xc_lda_vxc(&xc_x_func, &rhoptr[i], &tfptr[i], &tdf_drhoptr[i]);
//    }
//
//    f.gaxpy(1.0, tf, 1.0);
//    df_drho.gaxpy(1.0, tdf_drho, 1.0);
//}
//  //***************************************************************************
//
//  //***************************************************************************
//  template <int NDIM>
// inline void xc_lda_V(const Key<NDIM>& key, Tensor<double>& t)
//  {
//    Tensor<double> enefunc = copy(t);
//    Tensor<double> V = copy(t);
//    ::xc_generic_lda(t, enefunc, V, false);
//    t(___) = V(___);
//  }
//  //***************************************************************************
//
//  //***************************************************************************
//  template <int NDIM>
// inline void xc_lda_ene(const Key<NDIM>& key, Tensor<double>& t)
//  {
//    Tensor<double> V = copy(t);
//    Tensor<double> enefunc = copy(t);
//    ::xc_generic_lda(t, enefunc, V, false);
//    t(___) = enefunc(___);
//  }
//  //***************************************************************************



#endif /* LIBXC_H_ */
