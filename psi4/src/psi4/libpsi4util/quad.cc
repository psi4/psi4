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

#include "quad.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "psi4/libparallel/ParallelPrinter.h"
using namespace psi;

namespace psi {

Quadrature::Quadrature(int npoints)
{

    npoints_ = npoints;
    index_ = 0;
    w_ = new double[npoints];
    t_ = new double[npoints];
}

Quadrature::~Quadrature()
{
    delete[] w_;
    delete[] t_;
}

ChebyshevIIQuadrature::ChebyshevIIQuadrature(int npoints, double t0) :
  Quadrature(npoints) , center_(t0)
{
    // Compute Becke-style mapping
    // (Seems to span the space better)
    double x,temp;
    double INVLN2 = 1.0/log(2.0);
    for (int tau = 1; tau<=npoints_; tau++) {
        x = cos(tau/(npoints_+1.0)*M_PI);
        t_[tau-1] = t0*(1.0-x)/(1.0+x);
        temp = sin(tau/(npoints_+1.0)*M_PI);
        w_[tau-1] = 2.0*M_PI/(npoints_+1)*temp*temp*t0/((1.0+x)*(1.0+x)*
          sqrt(1.0-x*x));
    }
}

void ChebyshevIIQuadrature::print(std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   printer->Printf( "  Chebyshev Type II Quadrature of %d Points\n", npoints_);
    printer->Printf( "        for integration on [0, \\infty)\n");
    printer->Printf( "           Center %14.10E\n", center_);
    printer->Printf( "\n");
    printer->Printf( "  Index       Point         Weight\n");
    for (int k = 0; k < npoints_; k++)
        printer->Printf( "   %3d      %8.3E    %8.3E\n", k+1, t_[k], w_[k]);
}

}
