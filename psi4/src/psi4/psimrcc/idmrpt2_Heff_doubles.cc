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
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include "psi4/libpsi4util/libpsi4util.h"

#include "idmrpt2.h"
#include "blas.h"
#include "debugging.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

/**
 * \brief Computes the contribution (to be spin-factored)
 * \f[
 * <uv||xy>
 * + P(xy) \sum_{e} t_{uv}^{xe}(\mu) f_{ye}(\mu)
 * - P(uv) \sum_{m} t_{um}^{xy}(\mu) f_{mv}(\mu)
 * \f]
 * \f[
 * + \frac{1}{2} \sum_{mn} t_{mn}^{xy}(\mu) <mn||uv>
 * + \frac{1}{2} \sum_{ef} t_{uv}^{ef}(\mu) <xy||ef>
 * + P(uv)P(xy)  \sum_{me} t_{um}^{xe}(\mu) <my||ev>
 * \f]
 * \f[
 * + P(uv) \sum_{e} t_{u}^{e}(\mu) <xy||ev>
 * - P(xy) \sum_{m} t_{m}^{x}(\mu) <my||uv>
 * \f]
 */
void IDMRPT2::build_Heff_uvxy()
{
  START_TIMER(1,"Building the Heff_uvxy Matrix Elements");

  // Closed-shell
  blas->solve("Hijab[aa][aa]{c}  = <[aa]:[aa]>");

  blas->solve("Hijab[aa][aa]{c} += #3124# - t2_vvoo[v][aaa]{c} 1@2 fock[a][v]{c}");
  blas->solve("Hijab[aa][aa]{c} += #4123#   t2_vvoo[v][aaa]{c} 1@2 fock[a][v]{c}");

  blas->solve("Hijab[aa][aa]{c} += #1342#   t2_oovv[o][aaa]{c} 1@1 fock[o][a]{c}");
  blas->solve("Hijab[aa][aa]{c} += #2341# - t2_oovv[o][aaa]{c} 1@1 fock[o][a]{c}");

  blas->solve("Hijab[aa][aa]{c} += 1/2  <[oo]:[aa]> 1@1 t2_oovv[oo][aa]{c}");

  blas->solve("Hijab[aa][aa]{c} += 1/2 t2_oovv[aa][vv]{c} 2@2 <[aa]:[vv]>");

  blas->solve("Hijab[aa][aa]{c} += #1342#   t2_ovov[aa][ov]{c} 2@1 ([ov]:[aa])");
  blas->solve("Hijab[aa][aa]{c} += #1432# - t2_ovov[aa][ov]{c} 2@1 ([ov]:[aa])");
  blas->solve("Hijab[aa][aa]{c} += #2341# - t2_ovov[aa][ov]{c} 2@1 ([ov]:[aa])");
  blas->solve("Hijab[aa][aa]{c} += #2431#   t2_ovov[aa][ov]{c} 2@1 ([ov]:[aa])");

  blas->solve("Hijab[aa][aa]{c} += #1342#   t2_ovOV[aa][OV]{c} 2@1 ([ov]|[aa])");
  blas->solve("Hijab[aa][aa]{c} += #1432# - t2_ovOV[aa][OV]{c} 2@1 ([ov]|[aa])");
  blas->solve("Hijab[aa][aa]{c} += #2341# - t2_ovOV[aa][OV]{c} 2@1 ([ov]|[aa])");
  blas->solve("Hijab[aa][aa]{c} += #2431#   t2_ovOV[aa][OV]{c} 2@1 ([ov]|[aa])");

  blas->solve("Hijab[aa][aa]{c} += #1234#   t1_ov[a][v]{c} 2@1 <[v]:[aaa]>");
  blas->solve("Hijab[aa][aa]{c} += #2134# - t1_ov[a][v]{c} 2@1 <[v]:[aaa]>");

  blas->solve("Hijab[aa][aa]{c} += #3412# - t1_ov[o][a]{c} 1@1 <[o]:[aaa]>");
  blas->solve("Hijab[aa][aa]{c} += #4312#   t1_ov[o][a]{c} 1@1 <[o]:[aaa]>");


  // Open-shell
  blas->solve("Hijab[aa][aa]{o}  = <[aa]:[aa]>");

  blas->solve("Hijab[aa][aa]{o} += #3124# - t2_vvoo[v][aaa]{o} 1@2 fock[a][v]{o}");
  blas->solve("Hijab[aa][aa]{o} += #4123#   t2_vvoo[v][aaa]{o} 1@2 fock[a][v]{o}");

  blas->solve("Hijab[aa][aa]{o} += #1342#   t2_oovv[o][aaa]{o} 1@1 fock[o][a]{o}");
  blas->solve("Hijab[aa][aa]{o} += #2341# - t2_oovv[o][aaa]{o} 1@1 fock[o][a]{o}");

  blas->solve("Hijab[aa][aa]{o} += 1/2  <[oo]:[aa]> 1@1 t2_oovv[oo][aa]{o}");

  blas->solve("Hijab[aa][aa]{o} += 1/2 t2_oovv[aa][vv]{o} 2@2 <[aa]:[vv]>");

  blas->solve("Hijab[aa][aa]{o} += #1342#   t2_ovov[aa][ov]{o} 2@1 ([ov]:[aa])");
  blas->solve("Hijab[aa][aa]{o} += #1432# - t2_ovov[aa][ov]{o} 2@1 ([ov]:[aa])");
  blas->solve("Hijab[aa][aa]{o} += #2341# - t2_ovov[aa][ov]{o} 2@1 ([ov]:[aa])");
  blas->solve("Hijab[aa][aa]{o} += #2431#   t2_ovov[aa][ov]{o} 2@1 ([ov]:[aa])");

  blas->solve("Hijab[aa][aa]{o} += #1342#   t2_ovOV[aa][OV]{o} 2@1 ([ov]|[aa])");
  blas->solve("Hijab[aa][aa]{o} += #1432# - t2_ovOV[aa][OV]{o} 2@1 ([ov]|[aa])");
  blas->solve("Hijab[aa][aa]{o} += #2341# - t2_ovOV[aa][OV]{o} 2@1 ([ov]|[aa])");
  blas->solve("Hijab[aa][aa]{o} += #2431#   t2_ovOV[aa][OV]{o} 2@1 ([ov]|[aa])");

  blas->solve("Hijab[aa][aa]{o} += #1234#   t1_ov[a][v]{o} 2@1 <[v]:[aaa]>");
  blas->solve("Hijab[aa][aa]{o} += #2134# - t1_ov[a][v]{o} 2@1 <[v]:[aaa]>");

  blas->solve("Hijab[aa][aa]{o} += #3412# - t1_ov[o][a]{o} 1@1 <[o]:[aaa]>");
  blas->solve("Hijab[aa][aa]{o} += #4312#   t1_ov[o][a]{o} 1@1 <[o]:[aaa]>");

  DEBUGGING(3, blas->print("Hijab[aa][aa]{u}"); );

  END_TIMER(1);
}

void IDMRPT2::build_Heff_uVxY()
{
  START_TIMER(1,"Building the Heff_uVxY Matrix Elements");

  // Closed-shell
  blas->solve("HiJaB[aA][aA]{c}  = <[aa]|[aa]>");
  blas->solve("HiJaB[aA][aA]{c} += #3214# t2_VvOo[V][aAa]{c} 1@2 fock[a][v]{c}");
  blas->solve("HiJaB[aA][aA]{c} += #4123# t2_vVoO[v][AaA]{c} 1@2 fock[a][v]{c}");
  blas->solve("HiJaB[aA][aA]{c} += #1432# - t2_OoVv[O][aAa]{c} 1@1 fock[o][a]{c}");
  blas->solve("HiJaB[aA][aA]{c} += #2341# - t2_oOvV[o][AaA]{c} 1@1 fock[o][a]{c}");
  blas->solve("HiJaB[aA][aA]{c} += <[oo]|[aa]> 1@1 t2_oOvV[oO][aA]{c}");
  blas->solve("HiJaB[aA][aA]{c} += t2_oOvV[aA][vV]{c} 2@2 <[aa]|[vv]>");

  blas->solve("HiJaB[aA][aA]{c} += #1342# t2_ovov[aa][ov]{c} 2@1 ([ov]|[aa])");
  blas->solve("HiJaB[aA][aA]{c} += #1342# t2_ovOV[aa][OV]{c} 2@1 ([ov]:[aa])");
  blas->solve("HiJaB[aA][aA]{c} += #1423# - t2_oVOv[aA][Ov]{c} 2@2 <[aa]|[ov]>");
  blas->solve("HiJaB[aA][aA]{c} += #2314# - t2_oVOv[oV][Aa]{c} 1@2 <[aa]|[ov]>");
  blas->solve("HiJaB[aA][aA]{c} += #2431# t2_OVOV[AA][OV]{c} 2@1 ([ov]|[aa])");
  blas->solve("HiJaB[aA][aA]{c} += #2431# t2_ovOV[ov][AA]{c} 1@1 ([ov]:[aa])");

  blas->solve("HiJaB[aA][aA]{c} += #1234#   t1_ov[a][v]{c} 2@1 <[v]|[aaa]>");
  blas->solve("HiJaB[aA][aA]{c} += #2143#   t1_OV[A][V]{c} 2@1 <[v]|[aaa]>");
  blas->solve("HiJaB[aA][aA]{c} += #3412# - t1_ov[o][a]{c} 1@1 <[o]|[aaa]>");
  blas->solve("HiJaB[aA][aA]{c} += #4321# - t1_OV[O][A]{c} 1@1 <[o]|[aaa]>");

  // Open-shell
  blas->solve("HiJaB[aA][aA]{o}  = <[aa]|[aa]>");
  blas->solve("HiJaB[aA][aA]{o} += #3214# t2_VvOo[V][aAa]{o} 1@2 fock[A][V]{o}");
  blas->solve("HiJaB[aA][aA]{o} += #4123# t2_vVoO[v][AaA]{o} 1@2 fock[a][v]{o}");
  blas->solve("HiJaB[aA][aA]{o} += #1432# - t2_OoVv[O][aAa]{o} 1@1 fock[O][A]{o}");
  blas->solve("HiJaB[aA][aA]{o} += #2341# - t2_oOvV[o][AaA]{o} 1@1 fock[o][a]{o}");
  blas->solve("HiJaB[aA][aA]{o} += <[oo]|[aa]> 1@1 t2_oOvV[oO][aA]{o}");
  blas->solve("HiJaB[aA][aA]{o} += t2_oOvV[aA][vV]{o} 2@2 <[aa]|[vv]>");

  blas->solve("HiJaB[aA][aA]{o} += #1342# t2_ovov[aa][ov]{o} 2@1 ([ov]|[aa])");
  blas->solve("HiJaB[aA][aA]{o} += #1342# t2_ovOV[aa][OV]{o} 2@1 ([ov]:[aa])");
  blas->solve("HiJaB[aA][aA]{o} += #1423# - t2_oVOv[aA][Ov]{o} 2@2 <[aa]|[ov]>");
  blas->solve("HiJaB[aA][aA]{o} += #2314# - t2_oVOv[oV][Aa]{o} 1@2 <[aa]|[ov]>");
  blas->solve("HiJaB[aA][aA]{o} += #2431# t2_OVOV[AA][OV]{o} 2@1 ([ov]|[aa])");
  blas->solve("HiJaB[aA][aA]{o} += #2431# t2_ovOV[ov][AA]{o} 1@1 ([ov]:[aa])");

  blas->solve("HiJaB[aA][aA]{o} += #1234#   t1_ov[a][v]{o} 2@1 <[v]|[aaa]>");
  blas->solve("HiJaB[aA][aA]{o} += #2143#   t1_OV[A][V]{o} 2@1 <[v]|[aaa]>");
  blas->solve("HiJaB[aA][aA]{o} += #3412# - t1_ov[o][a]{o} 1@1 <[o]|[aaa]>");
  blas->solve("HiJaB[aA][aA]{o} += #4321# - t1_OV[O][A]{o} 1@1 <[o]|[aaa]>");

  DEBUGGING(3, blas->print("HiJaB[aA][aA]{u}"); );
  END_TIMER(1);
}

void IDMRPT2::build_Heff_UVXY()
{
  START_TIMER(1,"Building the Heff_UVXY Matrix Elements");

  // Closed-shell
  blas->solve("HIJAB[AA][AA]{c}  = Hijab[aa][aa]{c}");

  // Open-shell
  blas->solve("HIJAB[AA][AA]{o}  = <[aa]:[aa]>");
  blas->solve("HIJAB[AA][AA]{o} += #3124# - t2_VVOO[V][AAA]{o} 1@2 fock[A][V]{o}");
  blas->solve("HIJAB[AA][AA]{o} += #4123#   t2_VVOO[V][AAA]{o} 1@2 fock[A][V]{o}");
  blas->solve("HIJAB[AA][AA]{o} += #1342#   t2_OOVV[O][AAA]{o} 1@1 fock[O][A]{o}");
  blas->solve("HIJAB[AA][AA]{o} += #2341# - t2_OOVV[O][AAA]{o} 1@1 fock[O][A]{o}");
  blas->solve("HIJAB[AA][AA]{o} += 1/2  <[oo]:[aa]> 1@1 t2_OOVV[OO][AA]{o}");
  blas->solve("HIJAB[AA][AA]{o} += 1/2 t2_OOVV[AA][VV]{o} 2@2 <[aa]:[vv]>");
  blas->solve("HIJAB[AA][AA]{o} += #1342#   t2_OVOV[AA][OV]{o} 2@1 ([ov]:[aa])");
  blas->solve("HIJAB[AA][AA]{o} += #1432# - t2_OVOV[AA][OV]{o} 2@1 ([ov]:[aa])");
  blas->solve("HIJAB[AA][AA]{o} += #2341# - t2_OVOV[AA][OV]{o} 2@1 ([ov]:[aa])");
  blas->solve("HIJAB[AA][AA]{o} += #2431#   t2_OVOV[AA][OV]{o} 2@1 ([ov]:[aa])");
  blas->solve("HIJAB[AA][AA]{o} += #1342#   t2_ovOV[ov][AA]{o} 1@1 ([ov]|[aa])");
  blas->solve("HIJAB[AA][AA]{o} += #1432# - t2_ovOV[ov][AA]{o} 1@1 ([ov]|[aa])");
  blas->solve("HIJAB[AA][AA]{o} += #2341# - t2_ovOV[ov][AA]{o} 1@1 ([ov]|[aa])");
  blas->solve("HIJAB[AA][AA]{o} += #2431#   t2_ovOV[ov][AA]{o} 1@1 ([ov]|[aa])");
  blas->solve("HIJAB[AA][AA]{o} += #1234#   t1_OV[A][V]{o} 2@1 <[v]:[aaa]>");
  blas->solve("HIJAB[AA][AA]{o} += #2134# - t1_OV[A][V]{o} 2@1 <[v]:[aaa]>");
  blas->solve("HIJAB[AA][AA]{o} += #3412# - t1_OV[O][A]{o} 1@1 <[o]:[aaa]>");
  blas->solve("HIJAB[AA][AA]{o} += #4312#   t1_OV[O][A]{o} 1@1 <[o]:[aaa]>");

  DEBUGGING(3, blas->print("HIJAB[AA][AA]{u}"); );
  END_TIMER(1);
}

}} /* End Namespaces */
