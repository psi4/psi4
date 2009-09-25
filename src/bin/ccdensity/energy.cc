/*! \file
    \ingroup CCDENSITY
    \brief Calculates the one- and two-electron CC energies using the
    coresponding one- and two-particle density matrices.
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    /* ENERGY(): Compute the CC energy using the one- and two-particle
    ** density matrices.
    **
    ** E = sum_pq Dpq fpq + 1/4 sum_pqrs Gpqrs <pq||rs>
    **
    ** The individual two-electron components are:
    **
    ** E(ijkl) = 1/4 sum_ijkl Gijkl <ij||kl>
    **
    ** E(ijka) = 1/4 sum_ijka [Gijka <ij||ka> + Gijak <ij||ak> +
    **                         Giajk <ia||jk> + Gaijk <ai||jk>]
    **         = sum_ijka Gijka <ij||ka>
    **
    ** E(ijab) = 1/4 sum_ijab [Gijab <ij||ab> + Gabij <ab||ij>]
    **         = 1/2 sum_ijab Gijab <ij||ab>
    **
    ** E(iajb) = 1/4 sum_iajb [Giajb <ia||jb> + Giabj <ia||bj> +
    **                         Gaijb <ai||jb> + Gaibj <ai||bj>]
    **         = sum_iajb Giajb <ia||jb>
    **
    ** E(abci) = 1/4 sum_abci [Gabci <ab||ci> + Gabic <ab||ic> +
    **                         Gciab <ic||ab> + Gicab <ic||ab>]
    **         = sum_abci Gabci <ab||ci>
    ** E(abcd) = 1/4 sum_abcd Gabcd <ab||cd>
    **
    ** Individual spin cases are handled below.
    */

    void energy_RHF(struct RHO_Params rho_params);
    void energy_ROHF(struct RHO_Params rho_params);
    void energy_UHF(struct RHO_Params rho_params);

    void energy(struct RHO_Params rho_params)
    {
      if(params.ref == 0) energy_RHF(rho_params);
      else if(params.ref == 1) energy_ROHF(rho_params);
      else if(params.ref == 2) energy_UHF(rho_params);
    }

  }} // namespace psi::ccdensity
