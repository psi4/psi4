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

void d_spinad(void)
{
  dpdbuf4 D, D1;

  if(params.ref == 0) { /*** RHF ***/
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_scmcopy(&D, CC_DINTS, "D 2<ij|ab> - <ij|ba>", 2);
    dpd_buf4_sort_axpy(&D, CC_DINTS, pqsr, 0, 5, "D 2<ij|ab> - <ij|ba>", -1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_sort(&D, CC_DINTS, prqs, 10, 10, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    dpd_buf4_close(&D);

  }
}

}} // namespace psi::ccsort
