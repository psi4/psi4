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

#ifndef _DEFINES_H
#define	_DEFINES_H

#define GIT_ID "d12f233900069eb274854278e3aa1c733c34c9e6"

#define PSIF_DCFT_DPD 100
#define PSIF_DCFT_DENSITY 101
#define PRINT_ENERGY_COMPONENTS 0
#define ZERO 1.0E-16

#define ID(x) _ints->DPD_ID(x)

#ifndef INDEX
#define INDEX(i,j) (i > j ? i * (i + 1) / 2 + j : j + (j + 1) / 2 + i)
#endif

#endif	/* _DEFINES_H */