/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cmath>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include "iwl.h"
#include "iwl.hpp"

namespace psi {
  
void IWL::read_one(PSIO *psio, int itap, char *label, double *ints, 
    int ntri, int erase, int printflg, FILE *out)  
{
    int nmo;

    psio->open(itap, PSIO_OPEN_OLD);
    psio->read_entry(itap, label, (char *) ints, ntri*sizeof(double));
    psio->close(itap, !erase);

    if (printflg) {
        nmo = (int) (sqrt((double) (1 + 8 * ntri)) - 1)/2;
        print_array(ints, nmo, out);
    }
}

/*!
** IWL_RDONE()
**
** This function reads the one-electron integrals in the MO basis
**  from disk and stores them in core.  Substantially revised on
**  29 April 1998 to filter out frozen orbitals if requested.
**  This change requires a very different argument list from the 
**  previous version of this code.
**
** David Sherrill, January 1994
** Revised by David Sherrill, April 1998
** Revised by TDC, June 2001
**
**   \param itap       = tape to read ints from
**   \param label      = the PSIO label
**   \param ints       = buffer (already allocated) to store the integrals
**   \param ntri       = number of packed integrals
**   \param erase      = erase itap (1=yes, 0=no)
**   \param printflg   = printing flag.  Set to 1 to print ints; 
**                       otherwise, set to 0
**   \param out    = file pointer for output of ints or error messages
** \ingroup IWL
*/
int iwl_rdone(int itap, char *label, double *ints, int ntri, int erase, 
              int printflg, FILE *out)
{

  int nmo;

  psio_open(itap, PSIO_OPEN_OLD);
  psio_read_entry(itap, label, (char *) ints, ntri*sizeof(double));
  psio_close(itap, !erase);

  if (printflg) {
    nmo = (int) (sqrt((double) (1 + 8 * ntri)) - 1)/2;
    print_array(ints, nmo, out);
  }

  return(1);
}

}

