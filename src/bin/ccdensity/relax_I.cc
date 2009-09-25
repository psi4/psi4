/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* relax_I(): Add the orbital-response contributions from the
** one-electron density matrix to the I(I,J) and I(I,A) blocks of the
** Lagrangian.  These terms arise from the first-order CPHF
** equations.  I *think* the following code is general enough to deal
** with both RHF and ROHF cases.
** */

void relax_I_RHF(void);
void relax_I_ROHF(void);
void relax_I_UHF(void);

void relax_I(void)
{
  if(params.ref == 0) relax_I_RHF();
  else if(params.ref == 1) relax_I_ROHF();
  else if(params.ref == 2) relax_I_UHF();
}
 

}} // namespace psi::ccdensity
