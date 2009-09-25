/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/*
** DIIS: Direct inversion in the iterative subspace routine to
** accelerate convergence of the CCSD amplitude equations.
**
** Substantially improved efficiency of this routine:
** (1) Keeping at most two error vectors in core at once.
** (2) Limiting direct product (overlap) calculation to unique pairs.
** (3) Using LAPACK's linear equation solver DGESV instead of flin.
**
** -TDC  12/22/01
** -Modifications for ROHF and UHF, TDC, 6/03
*/

void diis_RHF(int);
void diis_ROHF(int);
void diis_UHF(int);

void diis(int iter)
{
  if(params.ref == 0) diis_RHF(iter);
  else if(params.ref == 1) diis_ROHF(iter);
  else if(params.ref == 2) diis_UHF(iter);

  return;
}
}} // namespace psi::ccenergy
