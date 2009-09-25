/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void e_sort(void)
{
  dpdbuf4 E;

  if(params.ref == 2) {  /** UHF **/
    /*** AA ***/
    /* <ij|ka> */
    dpd_buf4_init(&E, CC_EINTS, 0, 21, 0, 21, 0, 0, "E <AI|JK>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 0, 20, "E <IJ|KA>");
    dpd_buf4_close(&E);

    /* <ij||ka> (i>j,ka) */
    dpd_buf4_init(&E, CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 2, 20, "E <IJ||KA> (I>J,KA)");
    dpd_buf4_close(&E);

    /* <ij||ka> (i>j,ak) */
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_buf4_sort(&E, CC_EINTS, pqsr, 2, 21, "E <IJ||KA> (I>J,AK)");
    dpd_buf4_close(&E);

    /*** BB ***/
    /* <ij|ka> */
    dpd_buf4_init(&E, CC_EINTS, 0, 31, 10, 31, 10, 0, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 10, 30, "E <ij|ka>");
    dpd_buf4_close(&E);

    /* <ij||ka> (i>j,ka) */
    dpd_buf4_init(&E, CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 12, 30, "E <ij||ka> (i>j,ka)");
    dpd_buf4_close(&E);

    /* <ij||ka> (i>j,ak) */
    dpd_buf4_init(&E, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&E, CC_EINTS, pqsr, 12, 31, "E <ij||ka> (i>j,ak)");
    dpd_buf4_close(&E);

    /*** AB ***/
    /* <iJ|kA> */
    dpd_buf4_init(&E, CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
    dpd_buf4_sort(&E, CC_EINTS, qrsp, 23, 27, "E <iJ|kA>");
    dpd_buf4_close(&E);

    /* <Ij|Ak> */
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_buf4_sort(&E, CC_EINTS, qpsr, 22, 26, "E <Ij|Ak>");
    dpd_buf4_close(&E);

    /* <iJ|aK> */
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_sort(&E, CC_EINTS, qpsr, 23, 25, "E <iJ|aK>");
    dpd_buf4_close(&E);

    /* <Ia|Jk> */
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_sort(&E, CC_EINTS, rspq, 24, 22, "E <Ia|Jk>");
    dpd_buf4_close(&E);

  }
  else {  /** RHF/ROHF **/
    /* <ij|ka> */
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 0, 10, "E <ij|ka>");
    dpd_buf4_close(&E);

    /* <ij||ka> (i>j,ka) */
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, srqp, 2, 10, "E <ij||ka> (i>j,ka)");
    dpd_buf4_close(&E);

    /* <ij|ka> (ij,ak) */
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_sort(&E, CC_EINTS, pqsr, 0, 11, "E <ij|ka> (ij,ak)");
    dpd_buf4_close(&E);

    /* <ij||ka> (i>j,ak) */
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&E, CC_EINTS, pqsr, 2, 11, "E <ij||ka> (i>j,ak)");
    dpd_buf4_close(&E);

    /* <ia|jk> */
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, qpsr, 10, 0, "E <ia|jk>");
    dpd_buf4_close(&E);

    /* <ij|ak> */
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_sort(&E, CC_EINTS, rspq, 0, 11, "E <ij|ak>");
    dpd_buf4_close(&E);
  }
}

}} // namespace psi::ccsort
