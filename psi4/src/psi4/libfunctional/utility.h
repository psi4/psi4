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

#ifndef FUNCTIONAL_UTILITY
#define FUNCTIONAL_UTILITY

/// heaviside(x) = 1.0 if x >  0
//                 0.0 if x <= 0
// for matlab piecewise functions
inline double heaviside(double x)
{
    return (x > 0.0 ? 1.0 : 0.0);
}

/// dirac(...) = 0.0 for all x
// for matlab piecewise functions
inline double dirac(double x, ...)
{
    return 0.0;
}

double Ei(double x);
double expei(double x);
static const double expei_cutoff = 700.0;

#endif