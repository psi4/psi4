/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void f_sort(void)
{
  dpdbuf4 F;

  if(params.ref == 2) {  /*** UHF ***/
    dpd_buf4_init(&F, CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    dpd_buf4_sort_ooc(&F, CC_FINTS, spqr, 27, 29, "F <iA|bC>");
    dpd_buf4_close(&F);
  }
}

}} // namespace psi::ccsort
