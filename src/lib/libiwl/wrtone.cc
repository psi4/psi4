/*!
  \file wrtone.c
  \ingroup (IWL)
*/
#include <stdio.h>
#include <libpsio/psio.h>
#include "iwl.h"

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