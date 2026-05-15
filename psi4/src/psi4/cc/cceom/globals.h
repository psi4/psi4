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
void check_sum(const char *term_lbl, int index, int irrep);

//#define EOM_DEBUG (0)

// Use constexpr variables and inline functions instead of unsafe macros
inline constexpr int H_IRR = 0;
inline constexpr int MAX(int I, int J) { return (I > J) ? I : J; }
inline constexpr int MIN(int I, int J) { return (I < J) ? I : J; }

extern struct MOInfo moinfo;
extern struct Params params;
extern struct Eom_params eom_params;
extern struct Local local;
extern int ***dpd_dp;
}
}  // namespace psi

#endif // _psi_src_bin_cceom_globals_h
