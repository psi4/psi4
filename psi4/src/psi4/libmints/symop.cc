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

//
// symop.cc
//
// Additional modifications made by Justin Turney <jturney@ccqc.uga.edu>
// for use in PSI4.
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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
// the Free Software Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//


#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libparallel/ParallelPrinter.h"
#include "psi4/psi4-dec.h"

#include <cmath>

using namespace std;
using namespace psi;

////////////////////////////////////////////////////////////////////////

SymmetryOperation::SymmetryOperation() : bits_(0)
{
    zero();
}

SymmetryOperation::SymmetryOperation(const SymmetryOperation &so)
    : bits_(so.bits_)
{
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            d[i][j] = so.d[i][j];
        }
    }
}

SymmetryOperation::~SymmetryOperation()
{
}

SymmetryOperation & SymmetryOperation::operator = (SymmetryOperation const & so) // Assignment operator
{
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            d[i][j] = so.d[i][j];
        }
    }
    bits_ = so.bits_;
    return *this;
}

void SymmetryOperation::analyze_d()
{
#define EQUAL(x, y) (fabs((x) - (y)) < 1.0e-5)

    if (EQUAL(d[0][0], 1.0) && EQUAL(d[1][1], 1.0) && EQUAL(d[2][2], 1.0))
        bits_ = SymmOps::E;
    else if (EQUAL(d[0][0], 1.0) && EQUAL(d[1][1], -1.0) && EQUAL(d[2][2], -1.0))
        bits_ = SymmOps::C2_x;
    else if (EQUAL(d[0][0], -1.0) && EQUAL(d[1][1], 1.0) && EQUAL(d[2][2], -1.0))
        bits_ = SymmOps::C2_y;
    else if (EQUAL(d[0][0], -1.0) && EQUAL(d[1][1], -1.0) && EQUAL(d[2][2], 1.0))
        bits_ = SymmOps::C2_z;
    else if (EQUAL(d[0][0], 1.0) && EQUAL(d[1][1], 1.0) && EQUAL(d[2][2], -1.0))
        bits_ = SymmOps::Sigma_xy;
    else if (EQUAL(d[0][0], 1.0) && EQUAL(d[1][1], -1.0) && EQUAL(d[2][2], 1.0))
        bits_ = SymmOps::Sigma_xz;
    else if (EQUAL(d[0][0], -1.0) && EQUAL(d[1][1], 1.0) && EQUAL(d[2][2], 1.0))
        bits_ = SymmOps::Sigma_yz;
    else if (EQUAL(d[0][0], -1.0) && EQUAL(d[1][1], -1.0) && EQUAL(d[2][2], -1.0))
        bits_ = SymmOps::i;
}

SymmetryOperation
SymmetryOperation::operate(const SymmetryOperation& r) const
{
    SymmetryOperation ret;
    for (int i=0; i < 3; i++)
        for (int j=0; j < 3; j++) {
            double t=0;
            for (int k=0; k < 3; k++)
                t += r.d[i][k]*d[k][j];
            ret.d[i][j] = t;
        }
    ret.analyze_d();
    return ret;
}

SymmetryOperation
SymmetryOperation::transform(const SymmetryOperation& r) const
{
    int i,j,k;
    SymmetryOperation ret,foo;

    // foo = r * d
    for (i=0; i < 3; i++) {
        for (j=0; j < 3; j++) {
            double t=0;
            for (k=0; k < 3; k++)
                t += r.d[i][k] * d[k][j];
            foo.d[i][j] = t;
        }
    }

    // ret = (r*d)*r~ = foo*r~
    for (i=0; i < 3; i++) {
        for (j=0; j < 3; j++) {
            double t=0;
            for (k=0; k < 3; k++)
                t += foo.d[i][k]*r.d[j][k];
            ret.d[i][j]=t;
        }
    }

    ret.analyze_d();
    return ret;
}

// Clockwise rotation by 2pi/n degrees
void
SymmetryOperation::rotation(int n)
{
    double theta = (n) ? 2.0*M_PI/n : 2.0*M_PI;
    rotation(theta);
}

// Clockwise rotation by theta degrees
void
SymmetryOperation::rotation(double theta)
{
    zero();

    double ctheta = cos(theta);
    double stheta = sin(theta);

    d[0][0] = ctheta;
    d[0][1] = stheta;
    d[1][0] = -stheta;
    d[1][1] = ctheta;
    d[2][2] = 1.0;

    analyze_d();
}

void
SymmetryOperation::transpose()
{
    for (int i=1; i<3; i++) {
        for (int j=0; j<i; j++) {
            double tmp = d[i][j];
            d[i][j] = d[j][i];
            d[j][i] = tmp;
        }
    }
    analyze_d();
}

void
SymmetryOperation::print(std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   printer->Printf( "        1          2          3\n");
    printer->Printf( "  1  ");
    printer->Printf( "%10.7f ", d[0][0]);
    printer->Printf( "%10.7f ", d[0][1]);
    printer->Printf( "%10.7f \n", d[0][2]);
    printer->Printf( "  2  ");
    printer->Printf( "%10.7f ", d[1][0]);
    printer->Printf( "%10.7f ", d[1][1]);
    printer->Printf( "%10.7f \n", d[1][2]);
    printer->Printf( "  3  ");
    printer->Printf( "%10.7f ", d[2][0]);
    printer->Printf( "%10.7f ", d[2][1]);
    printer->Printf( "%10.7f \n", d[2][2]);
    outfile->Printf( "bits_ = %d\n", bits_);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
