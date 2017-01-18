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

/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "psi4/libpsi4util/libpsi4util.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void CCMRCC::build_Z_intermediates()
{
  // I am rewriting
  //   blas->append("Z_ijam[oov][o]{u} = #1234#  1/2 tau[oo][vv]{u} 2@2 <[vo]:[vv]>");
  // as
  blas->append("Z_ijam[oov][o]{u} = #1234#   tau[oo][vv]{u} 2@2 <[vo]|[vv]>");
  //

  blas->append("Z_iJaM[oOv][O]{u} = #1234#   tau[oO][vV]{u} 2@2 <[vo]|[vv]>");

  // I am rewriting
  //   blas->append("Z_iJAm[oOV][o]{u} = #1243# - tau[oO][vV]{u} 2@2 <[ov]|[vv]>");
  // as
  blas->append("Z_iJAm[oOV][o]{u} = #1234# - tau[oO][Vv]{u} 2@2 <[vo]|[vv]>");

  // I am rewriting
  // blas->append("Z_IJAM[OOV][O]{u} = #1234#  1/2 tau[OO][VV]{u} 2@2 <[vo]:[vv]>");
  // as
  blas->append("Z_IJAM[OOV][O]{u} = #1234#   tau[OO][VV]{u} 2@2 <[vo]|[vv]>");
  //
}

}} /* End Namespaces */
