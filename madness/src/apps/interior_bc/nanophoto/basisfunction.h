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
  
  $Id: test_problems.h 1856 2010-04-06 14:03:52Z mgr522 $
*/

/** \file basisfunction.h
    \brief

*/

#ifndef MADNESS_INTERIOR_BC_BASISFUNCTION_H__INCLUDED
#define MADNESS_INTERIOR_BC_BASISFUNCTION_H__INCLUDED

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <string>

using namespace madness;

/** \brief Abstract Gaussian basis function.  Overriden for angular momenta. */
class GaussianBF : public FunctionFunctorInterface<double, 3> {
    private:
        GaussianBF() {}

    protected:
        std::vector<double> contract_coeffs;
        std::vector<double> exponents;
        Vector<double, 3> center;
        int n;

    public:
        /// \brief Sets up the basis function data.
        GaussianBF(std::vector<double> &contract_coeffs,
                   std::vector<double> &exponents,
                   const Vector<double, 3> &c)
            : contract_coeffs(contract_coeffs), exponents(exponents),
              n(0) {

            center[0] = c[0];
            center[1] = c[1];
            center[2] = c[2];
        }

        virtual ~GaussianBF() {}

        std::vector<Vector<double, 3> > special_points() const {
            std::vector<Vector<double, 3> > pts(0);
            pts.push_back(center);
            return pts;
        }

        Level special_level() { return 7; }
};

/** \brief An s basis function. */
class SBF : public GaussianBF {
    protected:
        std::vector<double> norms;

    public:
        SBF(std::vector<double> &contract_coeffs,
            std::vector<double> &exponents,
				Vector<double, 3> &center)
            : GaussianBF(contract_coeffs, exponents, center),
              norms(contract_coeffs.size()) {

            int i;

            // print a warning if there are more coefficients than exponents
            n = contract_coeffs.size();
            i = exponents.size();
            if(n != i) {
                print("Error: different number of contraction coefficients " \
                    "and exponents.");
                // use as many as possible...
                if(i < n)
                    n = i;
            }

            for(i = 0; i < n; ++i) {
                norms[i] = pow(8.0*exponents[i]*exponents[i]*exponents[i] /
                    (constants::pi*constants::pi*constants::pi), 0.25);
            }
        }

        double operator() (const Vector<double, 3> &pt) const {
            double r2, ret;

            r2 = (center[0] - pt[0])*(center[0] - pt[0]) +
                 (center[1] - pt[1])*(center[1] - pt[1]) +
                 (center[2] - pt[2])*(center[2] - pt[2]);

            ret = 0.0;
            for(int i = 0; i < n; ++i)
                ret += contract_coeffs[i] * norms[i] *
                    exp(-exponents[i] * r2);

            return ret;
        }
};

/** \brief A p basis function. */
class PBF : public GaussianBF {
    protected:
        std::vector<double> norms;
        int dir;

    public:
        static const int X = 0, Y = 1, Z = 2;

        PBF(std::vector<double> &contract_coeffs,
            std::vector<double> &exponents,
				Vector<double, 3> &center, int dir)
            : GaussianBF(contract_coeffs, exponents, center),
              norms(contract_coeffs.size()), dir(dir) {

            int i;

            // print a warning if there are more coefficients than exponents
            n = contract_coeffs.size();
            i = exponents.size();
            if(n != i) {
                print("Error: different number of contraction coefficients " \
                    "and exponents.");
                // use as many as possible...
                if(i < n)
                    n = i;
            }

            for(i = 0; i < n; ++i) {
                norms[i] = pow(128.0*exponents[i]*exponents[i]*exponents[i]
                    *exponents[i]*exponents[i] /
                    (constants::pi*constants::pi*constants::pi), 0.25);
            }
        }

        double operator() (const Vector<double, 3> &pt) const {
            double r2, ret;
            Vector<double, 3> d;

            d[0] = pt[0] - center[0];
            d[1] = pt[1] - center[1];
            d[2] = pt[2] - center[2];

            r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];

            ret = 0.0;
            for(int i = 0; i < n; ++i)
                ret += contract_coeffs[i] * norms[i] *
                    exp(-exponents[i] * r2);

            switch(dir) {
            case X:
                ret *= d[0];
                break;
            case Y:
                ret *= d[1];
                break;
            case Z:
                ret *= d[2];
                break;
            }

            return ret;
        }
};

/** \brief A d basis function. */
class DBF : public GaussianBF {
    protected:
        std::vector<double> norms;
        int dir;

    public:
        static const int XX = 0, XY = 1, XZ = 2, YY = 3, YZ = 4, ZZ = 5;

        DBF(std::vector<double> &contract_coeffs,
            std::vector<double> &exponents,
				Vector<double, 3> &center, int dir)
            : GaussianBF(contract_coeffs, exponents, center),
              norms(contract_coeffs.size()), dir(dir) {

            int i;

            // print a warning if there are more coefficients than exponents
            n = contract_coeffs.size();
            i = exponents.size();
            if(n != i) {
                print("Error: different number of contraction coefficients " \
                    "and exponents.");
                // use as many as possible...
                if(i < n)
                    n = i;
            }

            for(i = 0; i < n; ++i) {
                norms[i] = 2048.0*exponents[i]*exponents[i]*exponents[i]
                    *exponents[i]*exponents[i]*exponents[i]*exponents[i] /
                    (constants::pi*constants::pi*constants::pi);
                if(dir == XX || dir == YY || dir == ZZ)
                    norms[i] /= 9.0;
                norms[i] = pow(norms[i], 0.25);
            }
        }

        double operator() (const Vector<double, 3> &pt) const {
            double r2, ret;
            Vector<double, 3> d;

            d[0] = pt[0] - center[0];
            d[1] = pt[1] - center[1];
            d[2] = pt[2] - center[2];

            r2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];

            ret = 0.0;
            for(int i = 0; i < n; ++i)
                ret += contract_coeffs[i] * norms[i] *
                    exp(-exponents[i] * r2);

            switch(dir) {
            case XX:
                ret *= d[0]*d[0];
                break;
            case XY:
                ret *= d[0]*d[1];
                break;
            case XZ:
                ret *= d[0]*d[2];
                break;
            case YY:
                ret *= d[1]*d[1];
                break;
            case YZ:
                ret *= d[1]*d[2];
                break;
            case ZZ:
                ret *= d[2]*d[2];
                break;
            }

            return ret;
        }
};

#endif // MADNESS_INTERIOR_BC_BASISFUNCTION_H__INCLUDED
