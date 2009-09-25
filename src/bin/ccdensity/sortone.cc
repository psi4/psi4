/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void sortone_RHF(struct RHO_Params rho_params);
void sortone_ROHF(struct RHO_Params rho_params);
void sortone_UHF(struct RHO_Params rho_params);

void sortone(struct RHO_Params rho_params)
{
  if(params.ref == 0) sortone_RHF(rho_params);
  else if(params.ref == 1) sortone_ROHF(rho_params);
  else if(params.ref == 2) sortone_UHF(rho_params);
}

}} // namespace psi::ccdensity
