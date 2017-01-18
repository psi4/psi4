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

#include "psi4/libmoinfo/libmoinfo.h"

#include "blas.h"
#include "debugging.h"
#include "idmrpt2.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{

/**
 * \brief Computes the contribution (to be spin-factored)
 * \f[
 * f_{uv}(\mu)
 * + \sum_{e} t_u^e(\mu) f_{ve}(\mu)
 * - \sum_{m} t_m^v(\mu) f_{mu}(\mu)
 * + \sum_{me} t_{um}^{ve}(\mu) f_{me}(\mu)
 * - \sum_{nf} t_{n}^{f}(\mu) <nv||uf>
 * - \frac{1}{2} \sum_{mef} t_{um}^{ef}(\mu) <mv||ef>
 * - \frac{1}{2} \sum_{mne} t_{mn}^{ve}(\mu) <nm||eu>
 * \f]
 */
void IDMRPT2::build_Heff_uv()
{
  START_TIMER(1,"Building the Heff_uv Matrix Elements");

  // Closed-shell
  blas->solve("Hia[a][a]{c}  = fock[a][a]{c}");
  blas->solve("Hia[a][a]{c} += t1_ov[a][v]{c} 2@2 fock[a][v]{c}");
  blas->solve("Hia[a][a]{c} += - fock[o][a]{c} 1@1 t1_ov[o][a]{c}");
  blas->solve("Hia[a][a]{c} += #12# t2_ovov[aa][ov]{c} 2@1 fock[ov]{c}");
  blas->solve("Hia[a][a]{c} += #12# t2_ovOV[aa][OV]{c} 2@1 fock[ov]{c}");
  blas->solve("Hia[a][a]{c} += #12# - <[aa]|[ov]> 2@1 t1[ov]{c}");
  blas->solve("Hia[a][a]{c} += #21# 2 ([ov]|[aa]) 1@1 t1[ov]{c}");
  blas->solve("Hia[a][a]{c} += 1/2 t2_oovv[a][ovv]{c} 2@2 <[a]:[ovv]>");
  blas->solve("Hia[a][a]{c} +=     t2_oOvV[a][OvV]{c} 2@2 <[a]|[ovv]>");
  blas->solve("Hia[a][a]{c} += -1/2 <[a]:[voo]> 2@2 t2_vvoo[a][voo]{c}");
  blas->solve("Hia[a][a]{c} += - <[a]|[voo]> 2@2 t2_vVoO[a][VoO]{c}");


  // Open-shell
  blas->solve("Hia[a][a]{o}  = fock[a][a]{o}");
  blas->solve("Hia[a][a]{o} += t1_ov[a][v]{o} 2@2 fock[a][v]{o}");
  blas->solve("Hia[a][a]{o} += - fock[o][a]{o} 1@1 t1_ov[o][a]{o}");
  blas->solve("Hia[a][a]{o} += #12# t2_ovov[aa][ov]{o} 2@1 fock[ov]{o}");
  blas->solve("Hia[a][a]{o} += #12# t2_ovOV[aa][OV]{o} 2@1 fock[OV]{o}");

  blas->solve("Hia[a][a]{o} += #12# - <[aa]|[ov]> 2@1 t1[ov]{o}");
  blas->solve("Hia[a][a]{o} += #21#   ([ov]|[aa]) 1@1 t1[ov]{o}");
  blas->solve("Hia[a][a]{o} += #21#   ([ov]|[aa]) 1@1 t1[OV]{o}");
  blas->solve("Hia[a][a]{o} += 1/2 t2_oovv[a][ovv]{o} 2@2 <[a]:[ovv]>");
  blas->solve("Hia[a][a]{o} +=     t2_oOvV[a][OvV]{o} 2@2 <[a]|[ovv]>");
  blas->solve("Hia[a][a]{o} += -1/2 <[a]:[voo]> 2@2 t2_vvoo[a][voo]{o}");
  blas->solve("Hia[a][a]{o} += - <[a]|[voo]> 2@2 t2_vVoO[a][VoO]{o}");

  DEBUGGING(3,blas->print("Hia[a][a]{u}"););

  END_TIMER(1);
}

/**
 * \brief Computes the contribution (to be spin-factored)
 * \f[
 * f_{UV}(\mu)
 * + \sum_{E} t_U^E(\mu) f_{VE}(\mu)
 * - \sum_{M} t_M^V(\mu) f_{MU}(\mu)
 * + \sum_{me} t_{Um}^{Ve}(\mu) f_{me}(\mu)
 * - \sum_{nf} t_{n}^{f}(\mu) <nV||Uf>
 * - \frac{1}{2} \sum_{mef} t_{Um}^{ef}(\mu) <mV||ef>
 * - \frac{1}{2} \sum_{mne} t_{mn}^{Ve}(\mu) <nm||eU>
 * \f]
 */
void IDMRPT2::build_Heff_UV()
{
  START_TIMER(1,"Building the Heff_UV Matrix Elements");

  // Closed-shell
  blas->solve("HIA[A][A]{c} = Hia[a][a]{c}");

  // Open-shell
  blas->solve("HIA[A][A]{o} = fock[A][A]{o}");
  blas->solve("HIA[A][A]{o} += t1_OV[A][V]{o} 2@2 fock[A][V]{o}");
  blas->solve("HIA[A][A]{o} += - fock[O][A]{o} 1@1 t1_OV[O][A]{o}");
  blas->solve("HIA[A][A]{o} += #12# t2_OVOV[AA][OV]{o} 2@1 fock[OV]{o}");
  blas->solve("HIA[A][A]{o} += #12# t2_ovOV[ov][AA]{o} 1@1 fock[ov]{o}");
  blas->solve("HIA[A][A]{o} += #12# - <[aa]|[ov]> 2@1 t1[OV]{o}");
  blas->solve("HIA[A][A]{o} += #21#   ([ov]|[aa]) 1@1 t1[OV]{o}");
  blas->solve("HIA[A][A]{o} += #21#   ([ov]|[aa]) 1@1 t1[ov]{o}");
  blas->solve("HIA[A][A]{o} += 1/2 t2_OOVV[A][OVV]{o} 2@2 <[a]:[ovv]>");
  blas->solve("HIA[A][A]{o} +=     t2_OoVv[A][oVv]{o} 2@2 <[a]|[ovv]>");
  blas->solve("HIA[A][A]{o} += -1/2 <[a]:[voo]> 2@2 t2_VVOO[A][VOO]{o}");
  blas->solve("HIA[A][A]{o} += - <[a]|[voo]> 2@2 t2_VvOo[A][vOo]{o}");

  DEBUGGING(3,blas->print("HIA[A][A]{u}"););

  END_TIMER(1);
}

}} /* End Namespaces */
