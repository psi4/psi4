/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <libutil/libutil.h>

#include "idmrpt2.h"
#include "blas.h"
#include "debugging.h"

extern FILE *outfile;

namespace psi{ namespace psimrcc{

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
