/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#define EXTERN
#include "globals.h"
namespace psi {
namespace ccdensity {

void V_build();
void Gijkl();
void Gabcd();
void Gibja();
void Gijka();
void Gciab();
void Gijab();

/* twopdm(): Computes all contributions to the two-particle density
** matrix for CC-like wave functions.
**
** Note that the contractions evaluated in the functions below
** actually build the bra-ket symmetrized two-particle density:
**
** Gamma'(pq,rs) = 1/2 [Gamma(pq,rs) + Gamma(rs,pq)],
**
** where Gamma(pq,rs) is the original, non-bra-ket-symmetric
** expression.  This is done to satisfy the
**
** TDC, July 2002
*/

void twopdm() {
    /*  V_build(); */
    Gijkl();
    Gabcd();
    Gijka();
    Gciab();
    Gibja();
    Gijab();
}

}  // namespace ccdensity
}  // namespace psi
