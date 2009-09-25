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

/* SORTI(): Place all the components of the Lagrangian into a large
** matrix, I (moinfo.I), which we also symmetrize by computing Ipq =
** 1/2 (Ipq + Iqp).  This matrix is later written to disk in dump()
** for subsequent backtransformation.  Note that some of the
** components of the Lagrangian computed into the IIJ, Iij, IIA, and
** Iia matrices remain non-symmetric (e.g., IIJ neq IJI).  I re-used
** my sortone.c code here, so don't let some of the variable names
** confuse you. */

void sortI_RHF(void);
void sortI_ROHF(void);
void sortI_UHF(void);

void sortI(void)
{
  if(params.ref == 0) sortI_RHF();
  else if(params.ref == 1) sortI_ROHF();
  else if(params.ref == 2) sortI_UHF();
}

}} // namespace psi::ccdensity
