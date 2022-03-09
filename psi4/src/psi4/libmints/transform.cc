/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "wavefunction.h"
#include "matrix.h"
#include "integral.h"

#include "psi4/psi4-dec.h"

#include <cmath>

using namespace psi;

extern void solidharmonic(int l, Matrix& coefmat);

// The ordering here is arbitrary and doesn't have to match the
// basis set ordering
static inline int ncart(int l) { return (l >= 0) ? ((((l) + 2) * ((l) + 1)) >> 1) : 0; }
static inline int npure(int l) { return 2 * l + 1; }
static inline int icart(int a, int b, int c) { return (((((a + b + c + 1) << 1) - a) * (a + 1)) >> 1) - b - 1; }

void SphericalTransformComponent::init(int a, int b, int c, double coef, int cartindex, int pureindex) {
    a_ = a;
    b_ = b;
    c_ = c;
    coef_ = coef;
    cartindex_ = cartindex;
    pureindex_ = pureindex;
}

SphericalTransform::SphericalTransform() {}

SphericalTransform::SphericalTransform(int l, int subl) : l_(l) {
    if (subl == -1)
        subl_ = l;
    else
        subl_ = subl;

    init();
}

void SphericalTransform::init() {
    int cartdim = INT_NCART(l_);
    Matrix coefmat(cartdim, cartdim);
    coefmat.zero();

    // Compute the solid harmonic matrix elements
    solidharmonic(l_, coefmat);

    // Go through and grab the values.
    int pureindex = 0;

    for (int i = 1; i <= (l_ - subl_) / 2; ++i) pureindex += npure(subl_ + 2 * i);

    // There is npure and an INT_NPURE macro that do the same thing.
    for (int p = 0; p < npure(subl_); ++p) {
        for (int a = 0; a <= l_; ++a) {
            for (int b = 0; (a + b) <= l_; ++b) {
                int c = l_ - a - b;

                int cart1 = icart(a, b, c);
                int cart2 = INT_CARTINDEX(a + b + c, a, b);

                double coef = coefmat(cart1, p + pureindex);

                if (std::fabs(coef) > 1.0e-16) {
                    SphericalTransformComponent component;
                    component.init(a, b, c, coef, cart2, p);
                    components_.push_back(component);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

ISphericalTransform::ISphericalTransform() : SphericalTransform() {}

ISphericalTransform::ISphericalTransform(int l, int subl) : SphericalTransform(l, subl) {
    components_.clear();
    init();
}

void ISphericalTransform::init() {
    int cartdim = ncart(l_);
    Matrix coefmat(cartdim, cartdim);
    coefmat.zero();

    // Compute the solid harmonic matrix elements
    solidharmonic(l_, coefmat);

    // Invert and transpose the coefficient matrix
    coefmat.invert();
    coefmat.transpose_this();

    // Go through and grab the values.
    int pureindex = 0;

    for (int i = 1; i <= (l_ - subl_) / 2; ++i) pureindex += npure(subl_ + 2 * i);

    for (int p = 0; p < npure(subl_); ++p) {
        for (int a = 0; a <= l_; ++a) {
            for (int b = 0; (a + b) <= l_; ++b) {
                int c = l_ - a - b;

                int cart1 = icart(a, b, c);
                int cart2 = INT_CARTINDEX(a + b + c, a, b);

                double coef = coefmat(cart1, p + pureindex);

                if (std::fabs(coef) > 1.0e-16) {
                    SphericalTransformComponent component;
                    component.init(a, b, c, coef, cart2, p);
                    components_.push_back(component);
                }
            }
        }
    }
}
