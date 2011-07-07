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


  $Id: miketest.cc 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/


#include <iostream>
#include <cstdio>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "tensor.h"
using madness::Tensor;


int main() {

    int k = 10;
    int twok = 2*k;
    double ops = 2*3*twok*twok*twok*twok;
    long times = 10000;
    double million = 1e6;
    double used;
    double mops;
    double start;

    Tensor<double> x = Tensor<double>(2*k,2*k,2*k);
    Tensor<double> r = Tensor<double>(2*k,2*k,2*k);
    Tensor<double> w = Tensor<double>(2*k,2*k,2*k);
    Tensor<double> c = Tensor<double>(2*k,2*k);

    start = std::clock();
    for (long i=0; i<times; ++i) {
        r = transform(x,c);
    }
    used = (double)(std::clock()-start)/(double)CLOCKS_PER_SEC;
    mops = ((double)times*ops)/(used*million);
    std::cout << "TRANSFORM MOPS=" << mops << "   "
    << used << "   "<< times << "   "<< ops << "   "<< million << std::endl;


    start = std::clock();
    for (long i=0; i<times; ++i) {
        fast_transform(x,c,r,w);
    }
    used = (double)(std::clock()-start)/(double)CLOCKS_PER_SEC;
    mops = ((double)times*ops)/(used*million);
    std::cout << "TRANSFORM MOPS=" << mops << "   "
    << used << "   "<< times << "   "<< ops << "   "<< million << std::endl;

    start = std::clock();
    for (long i=0; i<times; ++i) {
        r = transform3d(x,c);
    }
    used = (double)(std::clock()-start)/(double)CLOCKS_PER_SEC;
    mops = ((double)times*ops)/(used*million);
    std::cout << "TRANSFORM MOPS=" << mops << "   "
    << used << "   "<< times << "   "<< ops << "   "<< million << std::endl;

}
