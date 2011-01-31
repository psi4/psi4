#include "quad.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>

using namespace psi;

namespace psi {

Quadrature::Quadrature(int npoints) {

    npoints_ = npoints;
    index_ = 0;
    w_ = new double[npoints]; 
    t_ = new double[npoints]; 
}

Quadrature::~Quadrature() {
		delete[] w_;
		delete[] t_;
}

ChebyshevIIQuadrature::ChebyshevIIQuadrature(int npoints, double t0) :
  Quadrature(npoints) {
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

void ChebyshevIIQuadrature::print(FILE* out)
{
    fprintf(out, "  Chebyshev Type II Quadrature of %d Points\n", npoints_);
    fprintf(out, "        for integration on [0, \\infty)\n");
    fprintf(out, "\n");
    fprintf(out, "  Index     Point     Weight\n");
    for (int k = 0; k < npoints_; k++)
        fprintf(out, "   %d      %8.3E    %8.3E\n", k+1, t_[k], w_[k]);
}

}
