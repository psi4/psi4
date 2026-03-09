/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

void ltdensity_rohf(const struct TD_Params& S);
void ltdensity_uhf(const struct TD_Params& S);
void ltdensity_intermediates(const struct TD_Params& S);
void sort_ltd_rohf(const struct TD_Params& S);
void sort_ltd_uhf(const struct TD_Params& S);
void rtdensity(const struct TD_Params& S);
void sort_rtd_rohf(const struct TD_Params& S);
void sort_rtd_uhf(const struct TD_Params& S);

void tdensity(const struct TD_Params& S) {
    if (params.ref == 0 || params.ref == 1) {
        ltdensity_rohf(S);
        sort_ltd_rohf(S);
        rtdensity(S);
        sort_rtd_rohf(S);
    } else if (params.ref == 2) {
        ltdensity_uhf(S);
        sort_ltd_uhf(S);
        rtdensity(S);
        sort_rtd_uhf(S);
    }

    return;
}

}  // namespace ccdensity
}  // namespace psi
