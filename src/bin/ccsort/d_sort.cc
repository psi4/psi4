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

void d_sort(void)
{
  dpdbuf4 D;

  if(params.ref == 2) { /*** UHF ***/
    /*** AA ***/
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 0, 5, 1, "D <IJ|AB>");
    dpd_buf4_copy(&D, CC_DINTS, "D <IJ||AB> (I>J,A>B)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 0, 5, 1, "D <IJ|AB>");
    dpd_buf4_copy(&D, CC_DINTS, "D <IJ||AB> (I>J,AB)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 5, 1, "D <IJ|AB>");
    dpd_buf4_copy(&D, CC_DINTS, "D <IJ||AB> (IJ,A>B)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 1, "D <IJ|AB>");
    dpd_buf4_copy(&D, CC_DINTS, "D <IJ||AB>");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    dpd_buf4_sort(&D, CC_DINTS, prqs, 20, 20, "D <IJ||AB> (IA,JB)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_buf4_sort(&D, CC_DINTS, pqsr, 20, 21, "D <IJ||AB> (IA,BJ)");
    dpd_buf4_close(&D);

    /*** BB ***/
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 10, 15, 1, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 10, 15, 1, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab> (i>j,ab)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 17, 10, 15, 1, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab> (ij,a>b)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 10, 15, 1, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab>");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    dpd_buf4_sort(&D, CC_DINTS, prqs, 30, 30, "D <ij||ab> (ia,jb)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_sort(&D, CC_DINTS, pqsr, 30, 31, "D <ij||ab> (ia,bj)");
    dpd_buf4_close(&D);

    /*** AB ***/
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_sort(&D, CC_DINTS, qpsr, 23, 29, "D <iJ|aB>");
    dpd_buf4_sort(&D, CC_DINTS, psrq, 24, 26, "D <Ij|Ab> (Ib,Aj)");
    dpd_buf4_sort(&D, CC_DINTS, prqs, 20, 30, "D <Ij|Ab> (IA,jb)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_buf4_sort(&D, CC_DINTS, rspq, 30, 20, "D <Ij|Ab> (ia,JB)");
    dpd_buf4_sort(&D, CC_DINTS, pqsr, 20, 31, "D <Ij|Ab> (IA,bj)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_buf4_sort(&D, CC_DINTS, pqsr, 30, 21, "D <Ij|Ab> (ia,BJ)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_buf4_sort(&D, CC_DINTS, psrq, 27, 25, "D <iJ|aB> (iB,aJ)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_buf4_sort(&D, CC_DINTS, pqsr, 24, 27, "D <Ij|Ab> (Ib,jA)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_buf4_sort(&D, CC_DINTS, pqsr, 27, 24, "D <iJ|aB> (iB,Ja)");
    dpd_buf4_close(&D);
  }
  else {  /*** RHF/ROHF ***/
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 0, 5, 1, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 0, 5, 1, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab> (i>j,ab)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 5, 1, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab> (ij,a>b)");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 1, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_DINTS, "D <ij||ab>");
    dpd_buf4_close(&D);

    /* <ij|ab> (ia,jb) */
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_sort(&D, CC_DINTS, prqs, 10, 10, "D <ij|ab> (ia,jb)");
    dpd_buf4_close(&D);

    /* <ij|ab> (ai,jb) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_sort(&D, CC_DINTS, qprs, 11, 10, "D <ij|ab> (ai,jb)");
    dpd_buf4_close(&D);

    /* <ij|ab> (aj,ib) */
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_sort(&D, CC_DINTS, rqps, 11, 10, "D <ij|ab> (aj,ib)");
    dpd_buf4_close(&D);

    /* <ij|ab> (bi,ja) */
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_sort(&D, CC_DINTS, spqr, 11, 10, "D <ij|ab> (bi,ja)");
    dpd_buf4_close(&D);

    /* <ij||ab> (ia,jb) */
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    dpd_buf4_sort(&D, CC_DINTS, prqs, 10, 10, "D <ij||ab> (ia,jb)");
    dpd_buf4_close(&D);
  
    /* <ij|ab> (ib,ja) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_sort(&D, CC_DINTS, psrq, 10, 10, "D <ij|ab> (ib,ja)");
    dpd_buf4_close(&D);

    /* <ij|ab> (ib,aj) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_sort(&D, CC_DINTS, pqsr, 10, 11, "D <ij|ab> (ib,aj)");
    dpd_buf4_close(&D);

    /* <ij|ab> (ia,bj) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_sort(&D, CC_DINTS, pqsr, 10, 11, "D <ij|ab> (ia,bj)");
    dpd_buf4_close(&D);

    /* <ij||ab> (ia,bj) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_sort(&D, CC_DINTS, pqsr, 10, 11, "D <ij||ab> (ia,bj)");
    dpd_buf4_close(&D);

    /* <ib|aj> (ib,aj) */
    /* just use <ij|ab> (ib,aj), dummy
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_sort(&D, CC_DINTS, psrq, 10, 11, "D <ib|aj> (ib,aj)");
    dpd_buf4_close(&D);
    */
  }
}

}} // namespace psi::ccsort
