/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

namespace psi {


void dpd_error(const char *caller, FILE *out)
{
  fprintf(out, "Error in: %s\n", caller);
  dpd_close(dpd_default);
  exit(PSI_RETURN_FAILURE);
}

} // namespace psi
