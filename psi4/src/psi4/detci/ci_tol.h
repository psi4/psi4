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
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

#ifndef _psi_src_bin_detci_ci_tol_h
#define _psi_src_bin_detci_ci_tol_h

namespace psi { namespace detci {

#define S_MAX 0.999995 /* Max overlap for mitrush vects and collapsed vecs */
#define HD_MIN 1.0E-4  /* Minimum diagonal element of preconditioner       */
#define ZERO 1e-10
#define MPn_ZERO 1e-14
#define SA_NORM_TOL 1.0E-5 /* norm of schmidt orthogonalized vector */
#define MPn_NORM_TOL 1.0E-12

}} // namespace psi::detci

#endif // header guard