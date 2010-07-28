// shellrot.cc - renamed to shellrotation.cc
//
// Modified for PSI4 programming style.
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <cstdio>

#include <libmints/integral.h>
#include <libmints/pointgrp.h>
#include <libciomr/libciomr.h>
#include "shellrotation.h"

using namespace psi;

namespace psi {
    extern FILE *outfile;
}

ShellRotation::ShellRotation(int n)
    : n_(n), am_(0), r_(0)
{
    if (n_) {
        r_ = new double*[n_];
        for (int i=0; i<n_; ++i)
            r_[i] = new double[n_];
    }
}

ShellRotation::ShellRotation(const ShellRotation& other)
    : n_(0), am_(0), r_(0)
{
    *this = other;
}

ShellRotation::ShellRotation(int a, SymmetryOperation& so, const shared_ptr<IntegralFactory>& ints, int pure) :
    n_(0), am_(0), r_(0)
{
    init(a, so, ints);
}

ShellRotation::~ShellRotation()
{
    done();
}

ShellRotation& ShellRotation::operator=(const ShellRotation& other)
{
    done();

    n_ = other.n_;
    am_ = other.am_;

    if (n_ && other.r_) {
        r_ = new double*[n_];
        for (int i=0; i<n_; ++i) {
            r_[i] = new double[n_];
            memcpy(r_[i], other.r_[i], sizeof(double)*n_);
        }
    }

    return *this;
}

void ShellRotation::done()
{
    if (r_) {
        for (int i=0; i < n_; i++) {
            if (r_[i]) delete[] r_[i];
        }
        delete[] r_;
        r_=0;
    }
    n_=0;
}

void ShellRotation::init(int a, SymmetryOperation& so, const shared_ptr<IntegralFactory>& ints)
{
    done();

    am_ = a;

    if (a == 0) {
        n_ = 1;
        r_ = new double*[1];
        r_[0] = new double[1];
        r_[0][0] = 1.0;
        return;
    }

    CartesianIter* ip = ints->cartesian_iter(a);
    RedundantCartesianIter* jp = ints->redundant_cartesian_iter(a);

    CartesianIter& I = *ip;
    RedundantCartesianIter& J = *jp;
    int lI[3];
    int k, iI;

    n_ = I.n();
    r_ = new double*[n_];

    for (I.start(); I; I.next()) {
        r_[I.bfn()] = new double[n_];
        memset(r_[I.bfn()],0,sizeof(double)*n_);

        for (J.start(); J; J.next()) {
            double tmp = 1.0;

            for (k=0; k<3; ++k) {
                lI[k] = I.l(k);
            }

            for (k=0; k<am_; ++k) {
                for (iI=0; lI[iI]==0; iI++);
                lI[iI]--;
                double contrib = so(J.axis(k),iI);
                tmp *= contrib;
            }

            r_[I.bfn()][J.bfn()] += tmp;
        }
    }

    delete ip;
    delete jp;
}

ShellRotation ShellRotation::operate(const ShellRotation& rot) const
{
    if (n_ != rot.n_) {
        throw PSIEXCEPTION("ShellRotation::operate(): dimensions don't match.");
    }

    ShellRotation ret(n_);
    ret.am_ = am_;

    for (int i=0; i < n_; i++) {
        for (int j=0; j < n_; j++) {
            double t=0;
            for (int k=0; k < n_; k++)
                t += rot.r_[i][k] * r_[k][j];
            ret.r_[i][j] = t;
        }
    }

    return ret;
}

ShellRotation ShellRotation::transform(const ShellRotation& rot) const
{
    int i,j,k;

    if (rot.n_ != n_) {
        throw PSIEXCEPTION("ShellRotation::transform(): dimensions don't match.");
    }

    ShellRotation ret(n_), foo(n_);
    ret.am_ = foo.am_ = am_;

    // foo = r * d
    for (i=0; i < n_; i++) {
        for (j=0; j < n_; j++) {
            double t=0;
            for (k=0; k < n_; k++)
                t += rot.r_[i][k] * r_[k][j];
            foo.r_[i][j] = t;
        }
    }

    // ret = (r*d)*r~ = foo*r~
    for (i=0; i < n_; i++) {
        for (j=0; j < n_; j++) {
            double t=0;
            for (k=0; k < n_; k++)
                t += foo.r_[i][k]*rot.r_[j][k];
            ret.r_[i][j]=t;
        }
    }

    return ret;
}

double ShellRotation::trace() const
{
    double t=0;
    for (int i=0; i < n_; i++)
        t += r_[i][i];
    return t;
}

void ShellRotation::print() const
{
    print_mat(r_, n_, n_, outfile);
}
