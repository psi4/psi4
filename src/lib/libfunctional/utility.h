/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef FUNCTIONAL_UTILITY
#define FUNCTIONAL_UTILITY

#include <psiconfig.h>

#ifndef HAVE_FUNC_ERF
double erf(double x);
double erfc(double x);
#endif

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
