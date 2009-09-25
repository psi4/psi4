/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
/*
** READ_LAG.C
**
** Read the lagrangian
**
** C. David Sherrill
** University of California, Berkeley
**
** April 1998
** Updated to new libpsio libraries 8/03
*/

#include <cstdlib>
#include <cstdio>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include "globaldefs.h"
#include "globals.h"


namespace psi { namespace detcas {

void read_lagrangian(void)
{
  int nmo;

  nmo = CalcInfo.nmo;
  
  CalcInfo.lag = block_matrix(nmo, nmo);

  psio_open(Params.lag_file, PSIO_OPEN_OLD);  
  psio_read_entry(Params.lag_file, "MO-basis Lagrangian", 
    (char *) CalcInfo.lag[0], nmo*nmo*sizeof(double));

  if (Params.print_lvl > 3) {
    fprintf(outfile, "Lagrangian matrix\n");
    print_mat(CalcInfo.lag, nmo, nmo, outfile);
  }

  psio_close(Params.lag_file, Params.lag_erase ? 0 : 1);

} 

}} // end namespace psi::detcas

