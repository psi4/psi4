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

#include "psi4/liboptions/liboptions.h"

#include "mp2_ccsd.h"
#include "blas.h"
#include "debugging.h"

namespace psi{ namespace psimrcc{

void MP2_CCSD::build_t1_ia_amplitudes()
{
  blas->solve("t1_eqns[o][v]{u} = fock[o][v]{u}");
  blas->solve("t1_eqns[o][v]{u} += t1[o][v]{u} 2@2 F_ae[v][v]{u}");
  blas->solve("t1_eqns[o][v]{u} += - F_mi[o][o]{u} 1@1 t1[o][v]{u}");
  blas->solve("t1_eqns[o][v]{u} += #12# t2[ov][ov]{u} 2@1 F_me[ov]{u}");
  blas->solve("t1_eqns[o][v]{u} += #12# t2[ov][OV]{u} 2@1 F_ME[OV]{u}");

  blas->solve("t1_eqns[o][v]{u} += #12# - <[ov]|[ov]> 2@1 t1[ov]{u}");
  blas->solve("t1_eqns[o][v]{u} += #21#  ([ov]|[vo]) 1@1 t1[ov]{u}");
  blas->solve("t1_eqns[o][v]{u} += #21#  ([ov]|[vo]) 1@1 t1[OV]{u}");

  blas->solve("t1_eqns[o][v]{u} += 1/2 t2[o][ovv]{u} 2@2 <[v]:[ovv]>");
  blas->solve("t1_eqns[o][v]{u} +=     t2[o][OvV]{u} 2@2 <[v]|[ovv]>");

  blas->solve("t1_eqns[o][v]{u} += -1/2 <[o]:[voo]> 2@2 t2[v][voo]{u}");
  blas->solve("t1_eqns[o][v]{u} += - <[o]|[voo]> 2@2 t2[v][VoO]{u}");

  if(options_.get_str("MP2_CCSD_METHOD")=="I"){
    blas->reduce_spaces("t1_eqns[a][a]{u}","t1_eqns[o][v]{u}");
    blas->zero("t1_eqns[o][v]{u}");
    blas->expand_spaces("t1_eqns[a][a]{u}","t1_eqns[o][v]{u}");
  }

  blas->solve("t1_delta[o][v]{u} = t1_eqns[o][v]{u} / d1[o][v]{u} - t1[o][v]{u}");

  blas->solve("t1[o][v]{u} = t1_eqns[o][v]{u} / d1[o][v]{u}");
}

void MP2_CCSD::build_t1_IA_amplitudes()
{
//   START_TIMER(1,"Building the T1_IA Amplitudes");
//   blas->solve("t1[O][V]{u} = fock[O][V]{u} / d1[O][V]{u}");
//   END_TIMER(1);
  blas->solve("t1[O][V]{u} = t1[o][v]{u}");
}

}} /* End Namespaces */
