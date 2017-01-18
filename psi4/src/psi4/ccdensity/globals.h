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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include "psi4/psifiles.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"

namespace psi {

namespace ccdensity {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

/* #define DEBUG_XI (1)*/

EXTERN struct MOInfo moinfo;
EXTERN struct Frozen frozen;
EXTERN struct Params params;
EXTERN struct RHO_Params *rho_params;
EXTERN struct TD_Params *td_params;
//EXTERN std::vector<struct XTD_Params> xtd_params;

}} // namespace psi::ccdensity
