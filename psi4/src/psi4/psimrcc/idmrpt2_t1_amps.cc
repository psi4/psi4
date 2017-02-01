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

namespace psi{

    namespace psimrcc{

/**
 * \brief Computes the contribution
 * \f[
 * \mathcal{F}_{ai}(\mu)
 * + \sum_{e} t_i^e(\mu) \mathcal{F}_{ae}(\mu)
 * - \sum_{m} t_m^a(\mu) \mathcal{F}_{mi}(\mu)
 * + \sum_{uv} t_{iu}^{av}(\mu) \mathcal{F}_{uv}(\mu)
 * + \sum_{UV} t_{iU}^{aV}(\mu) \mathcal{F}_{UV}(\mu)
 * \f]
 */
void IDMRPT2::build_t1_ia_amplitudes()
{
  START_TIMER(1,"Building the T1_ia Amplitudes");
  blas->solve("t1_eqns[o][v]{u}  =   fock[o][v]{u}");
  blas->solve("t1_eqns[o][v]{u} +=     t1[o][v]{u} 2@2 F_ae[v][v]{u}");
  blas->solve("t1_eqns[o][v]{u} += - F_mi[o][o]{u} 1@1 t1[o][v]{u}");
  blas->solve("t1_eqns[o][v]{u} += #12# t2_ovov[aa][ov]{u} 1@1 fock[aa]{u}");
  blas->solve("t1_eqns[o][v]{u} += #12# t2_ovOV[ov][AA]{u} 2@1 fock[AA]{u}");
  END_TIMER(1);
}

/**
 * \brief Computes the contribution
 * \f[
 * \mathcal{F}_{AI}(\mu)
 * + \sum_{E} t_I^E(\mu) \mathcal{F}_{AE}(\mu)
 * - \sum_{M} t_M^A(\mu) \mathcal{F}_{MI}(\mu)
 * + \sum_{uv} t_{uI}^{vA}(\mu) \mathcal{F}_{uv}(\mu)
 * + \sum_{UV} t_{IU}^{AV}(\mu) \mathcal{F}_{UV}(\mu)
 * \f]
 */
void IDMRPT2::build_t1_IA_amplitudes()
{
  START_TIMER(1,"Building the T1_IA Amplitudes");
  // Closed-shell
  blas->solve("t1_eqns[O][V]{c} = t1_eqns[o][v]{c}");
  // Open-shell
  blas->solve("t1_eqns[O][V]{o}  =   fock[O][V]{o}");
  blas->solve("t1_eqns[O][V]{o} +=     t1[O][V]{o} 2@2 F_AE[V][V]{o}");
  blas->solve("t1_eqns[O][V]{o} += - F_MI[O][O]{o} 1@1 t1[O][V]{o}");
  blas->solve("t1_eqns[O][V]{o} += #12# t2_ovOV[aa][OV]{o} 1@1 fock[aa]{o}");
  blas->solve("t1_eqns[O][V]{o} += #12# t2_OVOV[AA][OV]{o} 1@1 fock[AA]{o}");
  END_TIMER(1);
}

}} /* End Namespaces */
