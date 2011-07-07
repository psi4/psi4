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
#ifndef UTIL_H_
#define UTIL_H_

#include <mra/mra.h>
#include <world/world.h>

namespace madness {
//  void printfunc(const World& world, Function<double,3> f, int npts)
//  {
//    Tensor<double> LL = FunctionDefaults<3>::get_cell_width();
//    double L = LL[0];
//    double bstep = L / npts;
//    f.reconstruct();
//    for (int i = 0; i <= npts; i++)
//    {
//      Vector<double,3> p(-L/2 + i * bstep);
//      if (world.rank() == 0) printf("%.2f\t\t%.8f\n", p[0], f(p));
//    }
//    if (world.rank() == 0) printf("\n");
//  }

//  void printfunc(const World& world, Function<double,3> f1, Function<double,3> f2, int npts)
//  {
//    Tensor<double> LL = FunctionDefaults<3>::get_cell_width();
//    double L = LL[0];
//    double bstep = L / npts;
//    f1.reconstruct();
//    f2.reconstruct();
//    for (int i = 0; i <= npts; i++)
//    {
//      Vector<double,3> p(-L/2 + i * bstep);
//      if (world.rank() == 0) printf("%.2f\t\t%.8f\t%.8f\n", p[0], f1(p), f2(p));
//    }
//    if (world.rank() == 0) printf("\n");
//  }
}
//
//#include <mra/mra.h>
//#include <world/world.h>
//#include <vector>
//
//namespace madness
//{
//  class OnesFunctor :
//  public FunctionFunctorInterface<double,3>
//  {
//  private:
//
//  public:
//    //*************************************************************************
//    OnesFunctor()
//    {
//    }
//    //*************************************************************************
//
//    //*************************************************************************
//    virtual ~OnesFunctor() {}
//    //*************************************************************************
//
//    //*************************************************************************
//    double operator()(const coordT& x) const
//    {
//      return 1.0;
//    }
//    //*************************************************************************
//  };
//
//  //***************************************************************************
//  class ZerosFunctor :
//  public FunctionFunctorInterface<double,3>
//  {
//  private:
//
//  public:
//    //*************************************************************************
//    ZerosFunctor()
//    {
//    }
//    //*************************************************************************
//
//    //*************************************************************************
//    virtual ~ZerosFunctor() {}
//    //*************************************************************************
//
//    //*************************************************************************
//    double operator()(const coordT& x) const
//    {
//      return 0.0;
//    }
//    //*************************************************************************
//  };
//  //***************************************************************************
//}

#endif
