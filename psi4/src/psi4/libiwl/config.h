/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifndef _psi3_libiwl_config_h_
#define _psi3_libiwl_config_h_

namespace psi {

// On-disk label width. The bucket format packs four labels per integral as
// 16-bit signed integers, so any orbital index passed through IWL must fit in
// [INT16_MIN, INT16_MAX]. The write paths check this and throw on overflow.
using Label = short int;
using Value = double;

inline constexpr const char *IWL_KEY_BUF = "IWL Buffers";
inline constexpr const char *IWL_KEY_ONEL = "IWL One-electron matrix elements";

// Number of integrals packed into a single on-disk bucket. Bucket size in
// bytes is 2*sizeof(int) + 4*N*sizeof(Label) + N*sizeof(Value), i.e. ~28 kB at
// N = 2980. Kept at the historical value for on-disk compatibility.
inline constexpr int IWL_INTS_PER_BUF = 2980;
}

#endif
