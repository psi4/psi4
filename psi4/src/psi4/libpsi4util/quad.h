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

#ifndef QUAD_H
#define QUAD_H

#include "psi4/psi4-dec.h"

namespace psi {

class Quadrature {
    protected:
        /// Number of points in this quadrature rule
        int npoints_;
        /// Current index in the quadrature rule
        int index_;
        /// Set of points (arbitrary domain)
        double* t_;
        /// Set of weights (Cartesian points, no spherical r^2)
        double* w_;
    public:
        /// Constructor, allocates memory
        Quadrature(int npoints);
        /// Destructor, frees memory
        virtual ~Quadrature();

        /// Prints the Quadrature rule
        virtual void print(std::string = "outfile") = 0;

        /// Get the current quadrature weight
        double getWeight() { return w_[index_]; }

        /// Get the current quadrature point
        double getPoint() { return t_[index_]; }

        /// Move to the next point
        void nextPoint() { index_++; }

        /// Reset the quadrature
        void reset() { index_ = 0; }

        /// See if the quadrature is complete
        bool isDone() { return index_ >= npoints_; }
};

class ChebyshevIIQuadrature : public Quadrature {
    protected:
        /// The center of the quadrature
        double center_;

    public:
        /**!
        * Constructor, develops quadrature rule
        * \param npoints Number of points.
        * \param t0 Center point of distribution (1.0 is default)
        */
        ChebyshevIIQuadrature(int npoints, double t0 = 1.0);
        /**!
        * Destructor, does nothing
        */
        ~ChebyshevIIQuadrature() {}

        /// Prints the Quadrature rule
        void print(std::string OutFileRMR = "outfile");

};

}
#endif
