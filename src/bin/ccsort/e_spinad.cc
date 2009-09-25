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

void e_spinad(void)
{
  dpdbuf4 E, E1;

  if(params.ref == 0) { /*** RHF ***/
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_scmcopy(&E, CC_EINTS, "E 2<ai|jk> - <ai|kj>", 2);
    dpd_buf4_sort(&E, CC_TMP0, pqsr, 11, 0, "E <ai|kj>");
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    dpd_buf4_init(&E1, CC_TMP0, 0, 11, 0, 11, 0, 0, "E <ai|kj>");
    dpd_buf4_axpy(&E1, &E, -1);
    dpd_buf4_close(&E1);
    dpd_buf4_close(&E);

  }
}

}} // namespace psi::ccsort
