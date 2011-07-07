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
#ifndef POPERATOR_H_
#define POPERATOR_H_

#include <constants.h>
#include "mra/operator.h"

#define WST_PI madness::constants::pi
//#define WST_PI 3.14159265358979323846264338328

namespace madness
{
  //***************************************************************************
  const double acut1e_6 = 0.25; //0.6450626287524907;
  //***************************************************************************

//  //***************************************************************************
//  template<typename Q, int NDIM>
//  SeparatedConvolution<Q, NDIM> PeriodicHFExchangeOp(World& world, long k,
//  double lo, double eps, Tensor<double> L, Vector<double,3> q)
//  {
//    // bsh_fit generates representation for 1/4Pir but we want 1/r
//    // so have to scale eps by 1/4Pi
//    Tensor<double> coeff, expnt;
//    //if (mu==0) eps /= 4.0*pi;
//    bsh_fit(0.0, lo, 100.0*L[0], eps/(4.0 * WST_PI), &coeff, &expnt, false); //eps /(4.0*pi)
//
//    // Scale coefficients according to the dimensionality and add to the list of operators
//    std::vector< std::shared_ptr< Convolution1D<Q> > > ops;
//    for (int i=0; i < coeff.dim(0); ++i)
//    {
//      if (expnt[i]*L[0]*L[0] > acut1e_6)
//      {
//        double c = pow(4 * WST_PI * coeff[i], 1.0/double(NDIM));
//        ops.push_back(std::shared_ptr< Convolution1D<Q> >(new GaussianConvolution1D<Q>(k,
//        	  c*L[0], expnt[i]*L[0]*L[0], 1.0, 0, true, q)));
//      }
//    }
//
//    return SeparatedConvolution<Q, NDIM>(world, ops);
//  }
//  //***************************************************************************

  //***************************************************************************
  template<typename Q, int NDIM>
  SeparatedConvolution<Q, NDIM> PeriodicCoulombOp(World& world, long k, double lo, double eps, Tensor<double> L)
  {
    // bsh_fit generates representation for 1/4Pir but we want 1/r
    // so have to scale eps by 1/4Pi
    Tensor<double> coeff, expnt;
    //if (mu==0) eps /= 4.0*pi;
    bsh_fit(0.0, lo, 100.0*L[0], eps/(4.0 * WST_PI), &coeff, &expnt, false); //eps /(4.0*pi)

    // Scale coefficients according to the dimensionality and add to the list of operators
    std::vector< std::shared_ptr< Convolution1D<Q> > > ops;
    for (int i=0; i < coeff.dim(0); ++i)
    {
      if (expnt[i]*L[0]*L[0] > acut1e_6)
      {
        double c = pow(4 * WST_PI * coeff[i], 1.0/double(NDIM));
        ops.push_back(std::shared_ptr< Convolution1D<Q> >(new GaussianConvolution1D<Q>(k, c*L[0], expnt[i]*L[0]*L[0], 1.0, 0, true)));
      }
    }

    return SeparatedConvolution<Q, NDIM>(world, ops);
  }
  //***************************************************************************

  //***************************************************************************
  template<typename Q, int NDIM>
  SeparatedConvolution<Q, NDIM>* PeriodicCoulombOpPtr(World& world, long k, double lo, double eps, Tensor<double> L)
  {
    // bsh_fit generates representation for 1/4Pir but we want 1/r
    // so have to scale eps by 1/4Pi
    Tensor<double> coeff, expnt;
    //if (mu==0) eps /= 4.0*pi;
    bsh_fit(0.0, lo, 100.0*L[0], eps/(4.0 * WST_PI), &coeff, &expnt, false); //eps /(4.0*pi)

    // Scale coefficients according to the dimensionality and add to the list of operators
    std::vector< std::shared_ptr< Convolution1D<Q> > > ops;
    for (int i=0; i < coeff.dim(0); ++i)
    {
      if (expnt[i]*L[0]*L[0] > acut1e_6)
      {
        double c = pow(4 * WST_PI * coeff[i], 1.0/double(NDIM));
        ops.push_back(std::shared_ptr< Convolution1D<Q> >(new GaussianConvolution1D<double>(k, c*L[0], expnt[i]*L[0]*L[0], 1.0, 0, true)));
      }
    }

    return new SeparatedConvolution<Q, NDIM>(world, ops);
  }
  //***************************************************************************

  //***************************************************************************
  template<typename Q, int NDIM>
  SeparatedConvolution<Q, NDIM> PeriodicBSHOp(World& world, double mu, long k, double lo, double eps, Tensor<double> L)
  {
    // bsh_fit generates representation for 1/4Pir but we want 1/r
    // so have to scale eps by 1/4Pi
    Tensor<double> coeff, expnt;
    //if (mu==0) eps /= 4.0*pi;
    bsh_fit(mu, lo, 10.0*L[0], eps, &coeff, &expnt, false); //eps /(4.0*pi)

    // Scale coefficients according to the dimensionality and add to the list of operators
    std::vector< std::shared_ptr< Convolution1D<Q> > > ops;
    for (int i=0; (i < coeff.dim(0)); ++i)
    {
      double c = pow(coeff[i], 1.0/double(NDIM));
      ops.push_back(std::shared_ptr< Convolution1D<Q> >(new GaussianConvolution1D<double>(k, c*L[0], expnt[i]*L[0]*L[0], 1.0, 0, true)));
    }

    return SeparatedConvolution<Q, NDIM>(world, ops);
  }
  //***************************************************************************

  //***************************************************************************
  template<typename Q, int NDIM>
  SeparatedConvolution<Q, NDIM>* PeriodicBSHOpPtr(World& world, double mu, long k, double lo, double eps, Tensor<double> L)
  {
    // bsh_fit generates representation for 1/4Pir but we want 1/r
    // so have to scale eps by 1/4Pi
    Tensor<double> coeff, expnt;
    //if (mu==0) eps /= 4.0*pi;
    bsh_fit(mu, lo, 10.0*L[0], eps, &coeff, &expnt, false); //eps /(4.0*pi)

    // Scale coefficients according to the dimensionality and add to the list of operators
    std::vector< std::shared_ptr< Convolution1D<Q> > > ops;
    for (int i=0; (i < coeff.dim(0)); ++i)
    {
      double c = pow(coeff[i], 1.0/double(NDIM));
      ops.push_back(std::shared_ptr< Convolution1D<Q> >(new GaussianConvolution1D<double>(k, c*L[0], expnt[i]*L[0]*L[0], 1.0, 0, true)));
    }

    return new SeparatedConvolution<Q, NDIM>(world, ops);
  }
  //***************************************************************************

};

#endif
