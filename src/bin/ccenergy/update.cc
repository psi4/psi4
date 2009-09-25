/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void update(void)
{
  fprintf(outfile,"  %4d      %20.15f    %4.3e    %7.6f    %7.6f    %7.6f\n",
          moinfo.iter,moinfo.ecc,moinfo.conv,moinfo.t1diag,moinfo.d1diag,moinfo.new_d1diag);
  fflush(outfile);
}
}} // namespace psi::ccenergy
