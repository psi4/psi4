/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void f_spinad(void)
{
  dpdbuf4 F, F1;

  if(params.ref == 0) { /*** RHF ***/
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_scmcopy(&F, CC_FINTS, "F 2<ia|bc> - <ia|cb>", 2);
    dpd_buf4_sort_ooc(&F, CC_TMP0, pqsr, 10, 5, "F <ia|cb>");
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F 2<ia|bc> - <ia|cb>");
    dpd_buf4_init(&F1, CC_TMP0, 0, 10, 5, 10, 5, 0, "F <ia|cb>");
    dpd_buf4_axpy(&F1, &F, -1);
    dpd_buf4_close(&F1);
    dpd_buf4_close(&F);

  }
}

}} // namespace psi::ccsort
