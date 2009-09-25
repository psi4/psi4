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

void Gijab_RHF(void);
void Gijab_ROHF(void);
void Gijab_UHF(void);

void Gijab(void)
{
  if(params.ref == 0) Gijab_RHF();
  else if(params.ref == 1) Gijab_ROHF();
  else if(params.ref == 2) Gijab_UHF();
}

}} // namespace psi::ccdensity
