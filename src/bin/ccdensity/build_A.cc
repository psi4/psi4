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

/* BUILD_A(): Construct the molecular orbital Hessian, A.
** */

void build_A_RHF(void);
void build_A_ROHF(void);
void build_A_UHF(void);

void build_A(void)
{
  if(params.ref == 0) build_A_RHF();
  else if(params.ref == 1) build_A_ROHF();
  else if(params.ref == 2) build_A_UHF();
}


}} // namespace psi::ccdensity
