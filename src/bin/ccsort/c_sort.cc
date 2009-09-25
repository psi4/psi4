/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void c_sort(void)
{
  dpdbuf4 C, D;

  if(params.ref == 2) { /** UHF **/

    /* <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - <ij|ba> */

    /*** AA ***/
    dpd_buf4_init(&C, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA|JB>");
    dpd_buf4_copy(&C, CC_CINTS, "C <IA||JB>");
    dpd_buf4_close(&C);
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ|AB>");
    dpd_buf4_sort(&D, CC_TMP0, psqr, 20, 20, "D <IJ|AB> (IB,JA)");
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_TMP0, 0, 20, 20, 20, 20, 0, "D <IJ|AB> (IB,JA)");
    dpd_buf4_init(&C, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    dpd_buf4_axpy(&D, &C, -1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&C);

    /* <IA||JB> (IA,BJ) (Wmbej.c) */
    dpd_buf4_init(&C, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    dpd_buf4_sort(&C, CC_CINTS, pqsr, 20, 21, "C <IA||JB> (IA,BJ)");
    dpd_buf4_close(&C);

    /*** BB ***/
    dpd_buf4_init(&C, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia|jb>");
    dpd_buf4_copy(&C, CC_CINTS, "C <ia||jb>");
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij|ab>");
    dpd_buf4_sort(&D, CC_TMP0, psqr, 30, 30, "D <ij|ab> (ib,ja)");
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_TMP0, 0, 30, 30, 30, 30, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_init(&C, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    dpd_buf4_axpy(&D, &C, -1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&C);

    /* <ia||jb> (ia,bj) (Wmbej.c) */
    dpd_buf4_init(&C, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    dpd_buf4_sort(&C, CC_CINTS, pqsr, 30, 31, "C <ia||jb> (ia,bj)");
    dpd_buf4_close(&C);

    /*** AB ***/

    /* <Ai|Bj> (iA,Bj) (Wmbej.c) */
    dpd_buf4_init(&C, CC_CINTS, 0, 26, 26, 26, 26, 0, "C <Ai|Bj>");
    dpd_buf4_sort(&C, CC_CINTS, qpsr, 27, 27, "C <iA|jB>");
    dpd_buf4_sort(&C, CC_CINTS, qprs, 27, 26, "C <Ai|Bj> (iA,Bj)");
    dpd_buf4_close(&C);

    /* <Ia|Jb> (Ia,bJ) (Wmbej.c) */
    dpd_buf4_init(&C, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    dpd_buf4_sort(&C, CC_CINTS, pqsr, 24, 25, "C <Ia|Jb> (Ia,bJ)");
    dpd_buf4_close(&C);

  }
  else { /** RHF/ROHF **/
    /* <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - <ij|ba> */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_copy(&C, CC_CINTS, "C <ia||jb>");
    dpd_buf4_close(&C);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_sort(&D, CC_TMP0, psqr, 10, 10, "D <ij|ab> (ib,ja)");
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_TMP0, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_buf4_axpy(&D, &C, -1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&C);

    /* <ia|jb> (bi,ja) */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_sort(&C, CC_CINTS, sprq, 11, 10, "C <ia|jb> (bi,ja)");
    dpd_buf4_close(&C);

    /* <ia||jb> (bi,ja) */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_buf4_sort(&C, CC_CINTS, sprq, 11, 10, "C <ia||jb> (bi,ja)");
    dpd_buf4_close(&C);

    /* <ia|jb> (ia,bj) */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_sort(&C, CC_CINTS, pqsr, 10, 11, "C <ia|jb> (ia,bj)");
    dpd_buf4_close(&C);

    /* <ia||jb> (ia,bj) (Wmbej.c) */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_buf4_sort(&C, CC_CINTS, pqsr, 10, 11, "C <ia||jb> (ia,bj)");
    dpd_buf4_close(&C);

    /* <ai|bj> (cchbar/Wabei_RHF.c) */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_sort(&C, CC_CINTS, qpsr, 11, 11, "C <ai|bj>");
    dpd_buf4_close(&C);

  }

}

}} // namespace psi::ccsort
