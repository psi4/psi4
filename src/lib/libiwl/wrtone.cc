/*!
  \file wrtone.c
  \ingroup (IWL)
*/
#include <stdio.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include "iwl.hpp"
#include "iwl.h"

  using namespace psi;
  
void IWL::wrtone(PSIO* psio, int itap, char *label, int ntri, double *onel_ints)
{
  psio->open(itap, PSIO_OPEN_OLD);
  psio->write_entry(itap, label, (char *) onel_ints, ntri*sizeof(double));
  psio->close(itap,1);
}

extern "C" {
/*!
** IWL_WRTONE()
**
** This function writes one-electron integrals.
**
**   itap       = tape to read ints from
**   label      = the PSIO label
**   ntri       = the size of the array (lower triangle)
**   onel_ints  = array to hold the one-electron integrals.
**
** David Sherrill, March 1995
** Revised by TDC, June 2001
** \ingroup (IWL)
*/
void iwl_wrtone(int itap, char *label, int ntri, double *onel_ints)
{
  psio_open(itap, PSIO_OPEN_OLD);
  psio_write_entry(itap, label, (char *) onel_ints, ntri*sizeof(double));
  psio_close(itap,1);
}

} /* extern "C" */