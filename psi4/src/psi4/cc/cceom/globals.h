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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/

#ifndef _psi_src_bin_cceom_globals_h
#define _psi_src_bin_cceom_globals_h

#include "psi4/psifiles.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {

namespace cceom {

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN void check_sum(const char *term_lbl, int index, int irrep);

//#define EOM_DEBUG (0)

#define H_IRR (0)
#define MAX(I, J) ((I > J) ? I : J)
#define MIN(I, J) ((I < J) ? I : J)

EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct Eom_params eom_params;
EXTERN struct Local local;
EXTERN int ***dpd_dp;
}
}  // namespace psi

#endif // _psi_src_bin_cceom_globals_h
