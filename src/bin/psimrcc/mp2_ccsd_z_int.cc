/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <libpsi4util/libpsi4util.h>
#include <libmoinfo/libmoinfo.h>

#include "blas.h"
#include "mp2_ccsd.h"
#include "debugging.h"
#include "matrix.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void MP2_CCSD::build_Z_intermediates()
{
  blas->solve("Z_iJaM[aAa][O]{u} = #1234#   tau_oOvV[aA][vV]{u} 2@2 <[ao]|[vv]>");
  blas->solve("Z_iJAm[aAA][o]{u} = #1234# - tau_oOVv[aA][Vv]{u} 2@2 <[ao]|[vv]>");

  blas->solve("Z_iJaM[oAa][O]{u} = #1234#   tau_oOvV[oA][vV]{u} 2@2 <[ao]|[vv]>");
  blas->solve("Z_iJAm[oAA][o]{u} = #1234# - tau_oOVv[oA][Vv]{u} 2@2 <[ao]|[vv]>");

  blas->solve("Z_iJaM[aAv][O]{u} = #1234#   tau_oOvV[aA][vV]{u} 2@2 <[vo]|[vv]>");
}

}} /* End Namespace*/