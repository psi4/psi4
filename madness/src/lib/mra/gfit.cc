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
#include <linalg/tensor_lapack.h>

using std::cout;
using std::endl;
using std::max;
using std::min;

/// \file gfit.cc
/// \brief Implements generation of separated operator representations using Gaussians

/// \ingroup function


namespace madness {
    using std::abs;

    static const double pi = 3.14159265358979323846264338328;


    // Returns in q[0..4] int(r^2(n+1)*exp(-alpha*r^2),r=0..R) n=-1,0,1,2
    static void gaussian_spherical_moments(double alpha, double R, Tensor<double>& q) {
        q[0] = -(-0.1e1 + exp(-alpha * R*R)) / alpha / 0.2e1;
        q[1] = (-0.2e1 * R * pow(alpha, 0.3e1 / 0.2e1) + sqrt(pi) * erf(R * sqrt(alpha)) * alpha * exp(alpha * R*R)) * pow(alpha, -0.5e1 / 0.2e1) * exp(-alpha * R*R) / 0.4e1;
        q[2] = -(-0.1e1 + exp(-alpha * R*R) + exp(-alpha * R*R) * alpha * R*R) * pow(alpha, -0.2e1) / 0.2e1;
        q[3] = -(-0.3e1 * sqrt(pi) * erf(R * sqrt(alpha)) * pow(alpha, 0.2e1) * exp(alpha * R*R) + 0.6e1 * R * pow(alpha, 0.5e1 / 0.2e1) + 0.4e1 * pow(R, 0.3e1) * pow(alpha, 0.7e1 / 0.2e1)) * pow(alpha, -0.9e1 / 0.2e1) * exp(-alpha * R*R) / 0.8e1;
    }

    // Returns in q[0..4] int(r^2(n+1)*exp(-mu*r)/(4*Pi*r),r=0..R) n=-1,0,1,2
    static void bsh_spherical_moments(double mu, double R, Tensor<double>& q) {
        if (mu == 0.0) {
            q[0] = R / pi  / 0.4e1;
            q[1] = pow(R, 0.2e1) / pi / 0.8e1;
            q[2] = pow(R, 0.3e1) / pi / 0.12e2;
            q[3] = pow(R, 0.4e1) / pi / 0.16e2;
        }
        else  {
            q[0] = (exp(mu * R) - 0.1e1) / mu * exp(-mu * R) / pi / 0.4e1;
            q[1] = -(-exp(mu * R) + 0.1e1 + mu * R) * pow(mu, -0.2e1) / pi * exp(-mu * R) / 0.4e1;
            q[2] = -(-0.2e1 * exp(mu * R) + 0.2e1 + 0.2e1 * mu * R + R*R * pow(mu, 0.2e1)) *
                   pow(mu, -0.3e1) / pi * exp(-mu * R) / 0.4e1;
            q[3] = -(-0.6e1 * exp(mu * R) + 0.6e1 + 0.6e1 * mu * R + 0.3e1 * R*R * pow(mu, 0.2e1) +
                     pow(R, 0.3e1) * pow(mu, 0.3e1)) * pow(mu, -0.4e1) / pi * exp(-mu * R) / 0.4e1;
        }
    }


    /// Form separated expansion for BSH (and Coulomb, mu=0) via trapezoidal quadrature.

    /// For Coulomb, prune, assuming hi = 1.0 (was 2sqrt(3))
    /// See Harrison et al in LNCS for details
    ///
    /// exp(-mu*r) / (4*pi*r) = sum(k) c[k]*exp(-t[k]*r^2)  + O(eps) for lo<=r<=hi
    void bsh_fit_old(double mu, double lo, double hi, double eps,
                     Tensor<double> *pcoeff, Tensor<double> *pexpnt, bool prnt) {


        if (mu > 0) {
            // Restrict hi according to the exponential decay
            double r = -log(4*pi*0.01*eps);
            r = -log(r * 4*pi*0.01*eps);
            if (hi > r) {
                //print("restrictING HI", hi, r);
                hi = r;
            }
        }

        double T;
        double slo, shi;

        if (eps >= 1e-2) T = 5;
        else if (eps >= 1e-4) T = 10;
        else if (eps >= 1e-6) T = 14;
        else if (eps >= 1e-8) T = 18;
        else if (eps >= 1e-10) T = 22;
        else if (eps >= 1e-12) T = 26;
        else T = 30;

        if (mu > 0) {
            slo = -0.5*log(4.0*T/(mu*mu));
        }
        else {
            slo = log(eps/hi) - 1.0;
        }

        shi = 0.5*log(T/(lo*lo));
        double h = 1.0/(.2-.50*log10(eps)); // was 0.47

        long npt = long((shi-slo)/h);
        h = (shi-slo)/(npt+1.0);

        //if (prnt)
        //cout << "slo " << slo << " shi " << shi << " npt " << npt << " h " << h << endl;

        Tensor<double> coeff(npt), expnt(npt);

        for (int i=0; i<npt; ++i) {
            double s = slo + h*(npt-i);	// i+1
            coeff[i] = h*2.0/sqrt(pi)*exp(-mu*mu*exp(-2.0*s)/4.0)*exp(s);
            coeff[i] = coeff[i]/(4.0*pi);
            expnt[i] = exp(2.0*s);
        }

        // Prune small exponents from Coulomb fit.  Evaluate a gaussian at
        // the range midpoint, and replace it there with the next most
        // diffuse gaussian.  Then examine the resulting error at the two
        // end points ... if this error is less than the desired
        // precision, can discard the diffuse gaussian.

        if (mu == 0.0) {
            double mid = sqrt(0.5*(hi*hi + lo*lo));
            long i;
            for (i=npt-1; i>0; --i) {
                double cnew = coeff[i]*exp(-(expnt[i]-expnt[i-1])*mid*mid);
                double errlo = coeff[i]*exp(-expnt[i]*lo*lo) -
                               cnew*exp(-expnt[i-1]*lo*lo);
                double errhi = coeff[i]*exp(-expnt[i]*hi*hi) -
                               cnew*exp(-expnt[i-1]*hi*hi);
                if (max(abs(errlo),abs(errhi)) > 0.03*eps) break;
                npt--;
                coeff[i-1] = coeff[i-1] + cnew;
            }
            coeff = coeff(Slice(0,npt-1));
            expnt = expnt(Slice(0,npt-1));
        }

        if (prnt) {
            for (int i=0; i<npt; ++i)
                cout << i << " " << coeff[i] << " " << expnt[i] << endl;

            long npt = 300;
            //double hi = 1.0;
            //if (mu) hi = min(1.0,30.0/mu);
            cout << "       x         value   abserr   relerr" << endl;
            cout << "  ------------  ------- -------- -------- " << endl;
            double step = exp(log(hi/lo)/(npt+1));
            for (int i=0; i<=npt; ++i) {
                double r = lo*(pow(step,i+0.5));
                double exact = exp(-mu*r)/r/4.0/pi;
                double test = 0.0;
                for (int j=0; j<coeff.dim(0); ++j)
                    test += coeff[j]*exp(-r*r*expnt[j]);
                double err = 0.0;
                if (exact) err = (exact-test)/exact;
                printf("  %.6e %8.1e %8.1e %8.1e\n",r, exact, exact-test, err);
            }
        }
        *pcoeff = coeff;
        *pexpnt = expnt;
    }

    void bsh_fit(double mu, double lo, double hi, double eps,
                 Tensor<double> *pcoeff, Tensor<double> *pexpnt, bool prnt) {

        if (mu > 0) {
            // Restrict hi according to the exponential decay
            double r = -log(4*pi*0.01*eps);
            r = -log(r * 4*pi*0.01*eps);
            if (hi > r) hi = r;
        }

        double T;
        double slo, shi;

        if (eps >= 1e-2) T = 5;
        else if (eps >= 1e-4) T = 10;
        else if (eps >= 1e-6) T = 14;
        else if (eps >= 1e-8) T = 18;
        else if (eps >= 1e-10) T = 22;
        else if (eps >= 1e-12) T = 26;
        else T = 30;

        if (mu > 0) {
            slo = -0.5*log(4.0*T/(mu*mu));
        }
        else {
            slo = log(eps/hi) - 1.0;
        }
        shi = 0.5*log(T/(lo*lo));

        // Resolution required for quadrature over s
        double h = 1.0/(0.2-.50*log10(eps)); // was 0.5 was 0.47

        // Truncate the number of binary digits in h's mantissa
        // so that rounding does not occur when performing
        // manipulations to determine the quadrature points and
        // to limit the number of distinct values in case of
        // multiple precisions being used at the same time.
        h = floor(64.0*h)/64.0;


        // Round shi/lo up/down to an integral multiple of quadrature points
        shi = ceil(shi/h)*h;
        slo = floor(slo/h)*h;

        long npt = long((shi-slo)/h+0.5);

        //if (prnt)
        //cout << "slo " << slo << " shi " << shi << " npt " << npt << " h " << h << " " << eps << endl;

        Tensor<double> coeff(npt), expnt(npt);

        for (int i=0; i<npt; ++i) {
            double s = slo + h*(npt-i);	// i+1
            coeff[i] = h*2.0/sqrt(pi)*exp(-mu*mu*exp(-2.0*s)/4.0)*exp(s);
            coeff[i] = coeff[i]/(4.0*pi);
            expnt[i] = exp(2.0*s);
        }

        // Prune large exponents from the fit ... never necessary due to construction

        // Prune small exponents from Coulomb fit.  Evaluate a gaussian at
        // the range midpoint, and replace it there with the next most
        // diffuse gaussian.  Then examine the resulting error at the two
        // end points ... if this error is less than the desired
        // precision, can discard the diffuse gaussian.

        if (mu == 0.0) {
            double mid = lo + (hi-lo)*0.5;
            long i;
            for (i=npt-1; i>0; --i) {
                double cnew = coeff[i]*exp(-(expnt[i]-expnt[i-1])*mid*mid);
                double errlo = coeff[i]*exp(-expnt[i]*lo*lo) -
                               cnew*exp(-expnt[i-1]*lo*lo);
                double errhi = coeff[i]*exp(-expnt[i]*hi*hi) -
                               cnew*exp(-expnt[i-1]*hi*hi);
                if (max(abs(errlo),abs(errhi)) > 0.03*eps) break;
                npt--;
                coeff[i-1] = coeff[i-1] + cnew;
            }
            coeff = coeff(Slice(0,npt-1));
            expnt = expnt(Slice(0,npt-1));
        }

        // Modify the coeffs of the largest exponents to satisfy the moment conditions
        //
        // SETTING NMOM>1 TURNS OUT TO BE A BAD IDEA (AS CURRENTLY IMPLEMENTED)
        // [It is accurate and efficient for a one-shot application but it seems to
        //  introduce fine-scale noise that amplifies during iterative solution of
        //  the SCF and DFT equations ... the symptom is negative coeffs in the fit]
        //
        // SET NMOM=0 or 1 (1 recommended) unless you are doing a one-shot application
        //
        // Determine the effective range of the four largest exponents and compute
        // moments of the exact and remainder of the fit.  Then adjust the coeffs
        // to reproduce the exact moments in that volume.
        //
        // If nmom!=4 we assume that we will eliminate n=-1 which is stored first
        // in the moment list
        //
        // <r^i|gj> cj = <r^i|exact-remainder>
        const long nmom = 1;
        if (nmom > 0) {
            Tensor<double> q(4), qg(4);
            double range = sqrt(-log(1e-6)/expnt[nmom-1]);
            if (prnt) print("exponent(nmom-1)",expnt[nmom-1],"has range", range);

            bsh_spherical_moments(mu, range, q);
            Tensor<double> M(nmom,nmom);
            for (int i=nmom; i<npt; ++i) {
                Tensor<double> qt(4);
                gaussian_spherical_moments(expnt[i], range, qt);
                qg += qt*coeff[i];
            }
            if (nmom != 4) {
                q = q(Slice(1,nmom));
                qg = qg(Slice(1,nmom));
            }
            if (prnt) {
                print("moments", q);
                print("moments", qg);
            }
            q = q - qg;
            for (int j=0; j<nmom; ++j) {
                Tensor<double> qt(4);
                gaussian_spherical_moments(expnt[j], range, qt);
                if (nmom != 4) qt = qt(Slice(1,nmom));
                for (int i=0; i<nmom; ++i) {
                    M(i,j) = qt[i];
                }
            }
            Tensor<double> ncoeff;
            gesv(M, q, ncoeff);
            if (prnt) {
                print("M\n",M);
                print("old coeffs", coeff(Slice(0,nmom-1)));
                print("new coeffs", ncoeff);
            }

            coeff(Slice(0,nmom-1)) = ncoeff;
        }

        if (prnt) {
            for (int i=0; i<npt; ++i)
                cout << i << " " << coeff[i] << " " << expnt[i] << endl;

            long npt = 300;
            //double hi = 1.0;
            //if (mu) hi = min(1.0,30.0/mu);
            cout << "       x         value   abserr   relerr" << endl;
            cout << "  ------------  ------- -------- -------- " << endl;
            double step = exp(log(hi/lo)/(npt+1));
            for (int i=0; i<=npt; ++i) {
                double r = lo*(pow(step,i+0.5));
                double exact = exp(-mu*r)/r/4.0/pi;
                double test = 0.0;
                for (int j=0; j<coeff.dim(0); ++j)
                    test += coeff[j]*exp(-r*r*expnt[j]);
                double err = 0.0;
                if (exact) err = (exact-test)/exact;
                printf("  %.6e %8.1e %8.1e %8.1e\n",r, exact, exact-test, err);
            }
        }
        *pcoeff = coeff;
        *pexpnt = expnt;
    }

    void bsh_fit_ndim(int ndim, double mu, double lo, double hi, double eps,
                      Tensor<double> *pcoeff, Tensor<double> *pexpnt, bool prnt) {

        if (mu > 0) {
            // Restrict hi according to the exponential decay
            double r = -log(4*pi*0.01*eps);
            r = -log(r * 4*pi*0.01*eps);
            if (hi > r) hi = r;
        }


        // Determine range of quadrature by estimating when
        // kernel drops to exp(-100)

        double slo, shi;
        if (mu > 0) {
            slo = -0.5*log(4.0*100.0/(mu*mu));
            slo = -0.5*log(4.0*(slo*ndim - 2.0*slo + 100.0)/(mu*mu));
        }
        else {
            slo = log(eps/hi) - 1.0;
        }
        shi = 0.5*log(100.0/(lo*lo));

        // Resolution required for quadrature over s
        double h = 1.0/(0.2-.50*log10(eps)); // was 0.5 was 0.47

        // Truncate the number of binary digits in h's mantissa
        // so that rounding does not occur when performing
        // manipulations to determine the quadrature points and
        // to limit the number of distinct values in case of
        // multiple precisions being used at the same time.
        h = floor(64.0*h)/64.0;


        // Round shi/lo up/down to an integral multiple of quadrature points
        shi = ceil(shi/h)*h;
        slo = floor(slo/h)*h;

        long npt = long((shi-slo)/h+0.5);

        if (prnt)
            cout << "bsh: mu " << mu << " lo " << lo << " hi " << hi << " eps " << eps << " slo " << slo << " shi " << shi << " npt " << npt << " h " << h << endl;


        // Compute expansion pruning small coeffs and large exponents
        Tensor<double> coeff(npt), expnt(npt);
        int nnpt=0;
        for (int i=0; i<npt; ++i) {
            double s = slo + h*(npt-i);	// i+1
            double c = exp(-0.25*mu*mu*exp(-2.0*s)+(ndim-2)*s)*0.5/pow(pi,0.5*ndim);
            double p = exp(2.0*s);
            c = c*h;
            if (c*exp(-p*lo*lo) > eps) {
                coeff(nnpt) = c;
                expnt(nnpt) = p;
                ++nnpt;
            }
        }
        npt = nnpt;


        // Prune small exponents from Coulomb fit.  Evaluate a gaussian at
        // the range midpoint, and replace it there with the next most
        // diffuse gaussian.  Then examine the resulting error at the two
        // end points ... if this error is less than the desired
        // precision, can discard the diffuse gaussian.

        if (mu == 0.0) {
            double mid = lo + (hi-lo)*0.5;
            long i;
            for (i=npt-1; i>0; --i) {
                double cnew = coeff[i]*exp(-(expnt[i]-expnt[i-1])*mid*mid);
                double errlo = coeff[i]*exp(-expnt[i]*lo*lo) -
                               cnew*exp(-expnt[i-1]*lo*lo);
                double errhi = coeff[i]*exp(-expnt[i]*hi*hi) -
                               cnew*exp(-expnt[i-1]*hi*hi);
                if (max(abs(errlo),abs(errhi)) > 0.03*eps) break;
                npt--;
                coeff[i-1] = coeff[i-1] + cnew;
            }
        }

        // Shrink array to correct size
        coeff = coeff(Slice(0,npt-1));
        expnt = expnt(Slice(0,npt-1));


        if (prnt) {
            for (int i=0; i<npt; ++i)
                cout << i << " " << coeff[i] << " " << expnt[i] << endl;

            long npt = 300;
            cout << "       x         value" << endl;
            cout << "  ------------  ---------------------" << endl;
            double step = exp(log(hi/lo)/(npt+1));
            for (int i=0; i<=npt; ++i) {
                double r = lo*(pow(step,i+0.5));
                double test = 0.0;
                for (int j=0; j<coeff.dim(0); ++j)
                    test += coeff[j]*exp(-r*r*expnt[j]);
                printf("  %.6e %20.10e\n",r, test);
            }
        }

        *pcoeff = coeff;
        *pexpnt = expnt;
    }

    void truncate_periodic_expansion(Tensor<double>& c, Tensor<double>& e,
      double L, bool discardG0)
    {
      double tcut = 0.25/L/L;

      if (discardG0)
      {
        // Relies on expnts being in decreasing order
        for (int i=0; i<e.dim(0); ++i) {
          if (e(i) < tcut) {
            c = c(Slice(0,i));
            e = e(Slice(0,i));
            break;
          }
        }
      }
      else
      {
        // Relies on expnts being in decreasing order
        int icut = -1;
        for (int i=0; i<e.dim(0); ++i) {
          if (e(i) < tcut) {
            icut = i;
            break;
          }
        }
        if (icut > 0) {
          for (int i=icut+1; i<e.dim(0); ++i) {
            c(icut) += c(i);
          }
          c = c(Slice(0,icut));
          e = e(Slice(0,icut));
        }
      }
    }
}
