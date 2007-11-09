/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

extern "C" {


void dpd_error(char *caller, FILE *outfile)
{
  fprintf(outfile, "Error in: %s\n", caller);
  dpd_close(dpd_default);
  exit(PSI_RETURN_FAILURE);
}

} /* extern "C" */
