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

#include "idmrpt2.h"
#include "blas.h"
#include "debugging.h"
#include "psi4/libmoinfo/libmoinfo.h"

namespace psi{

    namespace psimrcc{
    extern MOInfo *moinfo;

/**
 * \brief Computes the contribution
 * \f[
 * <iJ||aB>
 * + P(aB) \sum_{e} t_{iJ}^{ae}(\mu) \mathcal{F}_{Be}(\mu)
 * - P(iJ) \sum_{m} t_{im}^{aB}(\mu) \mathcal{F}_{mJ}(\mu)
 * \f]
 */
void IDMRPT2::build_t2_iJaB_amplitudes()
{
  START_TIMER(1,"Building the T2_iJaB Amplitudes");
  // Closed-shell
  blas->solve("t2_eqns[oO][vV]{c}  = <[oo]|[vv]>");
  blas->solve("t2_eqns[oO][vV]{c} += #3214# t2[V][vOo]{c} 1@2 F_ae[v][v]{c}");
  blas->solve("t2_eqns[oO][vV]{c} += #4123# t2[v][VoO]{c} 1@2 F_ae[v][v]{c}");
  blas->solve("t2_eqns[oO][vV]{c} += #1432# - t2[O][oVv]{c} 1@1 F_mi[o][o]{c}");
  blas->solve("t2_eqns[oO][vV]{c} += #2341# - t2[o][OvV]{c} 1@1 F_mi[o][o]{c}");

  // Open-shell
  blas->solve("t2_eqns[oO][vV]{o}  = <[oo]|[vv]>");
  blas->solve("t2_eqns[oO][vV]{o} += #3214# t2[V][vOo]{o} 1@2 F_AE[V][V]{o}");
  blas->solve("t2_eqns[oO][vV]{o} += #4123# t2[v][VoO]{o} 1@2 F_ae[v][v]{o}");
  blas->solve("t2_eqns[oO][vV]{o} += #1432# - t2[O][oVv]{o} 1@1 F_MI[O][O]{o}");
  blas->solve("t2_eqns[oO][vV]{o} += #2341# - t2[o][OvV]{o} 1@1 F_mi[o][o]{o}");
  END_TIMER(1);
}

/**
 * \brief Computes the contribution
 * \f[
 * <ij||ab>
 * + P(ab) \sum_{e} t_{ij}^{ae}(\mu) \mathcal{F}_{be}(\mu)
 * - P(ij) \sum_{m} t_{im}^{ab}(\mu) \mathcal{F}_{mj}(\mu)
 * \f]
 */
void IDMRPT2::build_t2_ijab_amplitudes()
{
  START_TIMER(1,"Building the T2_ijab Amplitudes");
  if(moinfo->get_ref_size(UniqueOpenShellRefs)==0){
    blas->solve("t2_eqns[oo][vv]{c}  = t2_eqns[oO][vV]{c}");
    blas->solve("t2_eqns[oo][vv]{c} += #2134# - t2_eqns[oO][vV]{c}");
  }else{
    // Closed-shell
    blas->solve("t2_eqns[oo][vv]{c}  = <[oo]:[vv]>");
    blas->solve("t2_eqns[oo][vv]{c} += #3124# - t2[v][voo]{c} 1@2 F_ae[v][v]{c}");
    blas->solve("t2_eqns[oo][vv]{c} += #4123#   t2[v][voo]{c} 1@2 F_ae[v][v]{c}");
    blas->solve("t2_eqns[oo][vv]{c} += #1342#   t2[o][ovv]{c} 1@1 F_mi[o][o]{c}");
    blas->solve("t2_eqns[oo][vv]{c} += #2341# - t2[o][ovv]{c} 1@1 F_mi[o][o]{c}");
  }

  // Open-shell
  blas->solve("t2_eqns[oo][vv]{o}  = <[oo]:[vv]>");
  blas->solve("t2_eqns[oo][vv]{o} += #3124# - t2[v][voo]{o} 1@2 F_ae[v][v]{o}");
  blas->solve("t2_eqns[oo][vv]{o} += #4123#   t2[v][voo]{o} 1@2 F_ae[v][v]{o}");
  blas->solve("t2_eqns[oo][vv]{o} += #1342#   t2[o][ovv]{o} 1@1 F_mi[o][o]{o}");
  blas->solve("t2_eqns[oo][vv]{o} += #2341# - t2[o][ovv]{o} 1@1 F_mi[o][o]{o}");
  END_TIMER(1);
}

/**
 * \brief Computes the contribution
 * \f[
 * <IJ||AB>
 * + P(AB) \sum_{e} t_{IJ}^{Ae}(\mu) \mathcal{F}_{Be}(\mu)
 * - P(IJ) \sum_{m} t_{Im}^{AB}(\mu) \mathcal{F}_{mJ}(\mu)
 * \f]
 */
void IDMRPT2::build_t2_IJAB_amplitudes()
{
  START_TIMER(1,"Building the T2_IJAB Amplitudes");
  // Closed-shell
  blas->solve("t2_eqns[OO][VV]{c}  = t2_eqns[oo][vv]{c}");

  // Open-shell
  blas->solve("t2_eqns[OO][VV]{o}  = <[oo]:[vv]>");
  blas->solve("t2_eqns[OO][VV]{o} += #3124# - t2[V][VOO]{o} 1@2 F_AE[V][V]{o}");
  blas->solve("t2_eqns[OO][VV]{o} += #4123#   t2[V][VOO]{o} 1@2 F_AE[V][V]{o}");
  blas->solve("t2_eqns[OO][VV]{o} += #1342#   t2[O][OVV]{o} 1@1 F_MI[O][O]{o}");
  blas->solve("t2_eqns[OO][VV]{o} += #2341# - t2[O][OVV]{o} 1@1 F_MI[O][O]{o}");
  END_TIMER(1);
}

}} /* End Namespaces */
