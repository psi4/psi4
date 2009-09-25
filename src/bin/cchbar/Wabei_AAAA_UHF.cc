/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* WABEI_UHF(): Computes all contributions to the ABEI spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (EI,AB) ordering and is referred to on disk as "WEIAB".
**
** The spin-orbital expression for the Wabei elements is:
**
** Wabei = <ab||ei> - Fme t_mi^ab + t_i^f <ab||ef>
**         - P(ab) t_m^b <am||ef>t_i^f + 1/2 tau_mn^ab <mn||ef> t_i^f
**         + 1/2 <mn||ei> tau_mn^ab - P(ab) <mb||ef> t_mi^af
**         - P(ab) t_m^a { <mb||ei> - t_ni^bf <mn||ef> }
**
** (cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).)
**
** For the ABEI spin case, we evaluate these contractions with two
** target orderings, (AB,EI) and (EI,AB), depending on the term.
** After all terms have been evaluated, the (AB,EI) terms are sorted
** into (EI,AB) ordering and both groups arer added together.
**
** TDC, June 2002
*/

void WABEI_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /**** Term I ****/

  /** W(EI,AB) <--- <EI||AB> **/
  dpd_buf4_init(&F, CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
  dpd_buf4_copy(&F, CC_HBAR, "WEIAB");
  dpd_buf4_close(&F);

  /**** Term II ****/

  /** W(EI,AB) <--- - F_ME t_MI^AB **/
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_buf4_init(&W, CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  dpd_contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_file2_close(&Fme);
  dpd_buf4_close(&T2);

  /**** Term III ****/

  /** W'(AB,EI) <--- <AB||EF> t_I^F **/
  dpd_buf4_init(&W, CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  dpd_buf4_init(&B, CC_BINTS, 0, 7, 5, 5, 5, 1, "B <AB|CD>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&B, &T1, &W, 3, 1, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&B);
  dpd_buf4_close(&W);

  /**** Term IV ****/

  /** WABEI <-- t_M^B <MA||EF> t_I^F - t_M^A <MB||EF> t_I^F
      Evaluate in two steps: 
          (1) Z_MBEI = <MB||EF> t_I^F 
          (2) WABEI <-- t_M^B Z_MAEI - t_M^A Z_MBEI
  **/

  /** Z(MB,EI) <-- - <MB||EF> t_I^F **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&F, &T1, &Z, 3, 1, 0, -1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** t_M^A Z(MB,EI) --> Z1(AB,EI) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(AB,EI)");
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort(&Z1, CC_TMP0, qprs, 5, 21, "Z2(BA,EI)");
  dpd_buf4_close(&Z1);

  /** Z1(AB,EI) - Z2(BA,EI) --> W'(AB,EI) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(AB,EI)");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z2(BA,EI)");
  dpd_buf4_axpy(&Z2, &Z1, -1.0);
  dpd_buf4_close(&Z2);

  dpd_buf4_init(&W, CC_TMP0, 0, 5, 21, 7, 21, 0, "W'(AB,EI)");
  dpd_buf4_axpy(&Z1, &W, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z1);

  /**** Term V ****/
  
  /** WABEI <-- 1/2 tau_MN^AB <MN||EF> t_I^F
      Evaluate in two steps:
         (1) Z_MNEI = <MN||EF> t_I^F
         (2) WABEI <-- 1/2 tau_MN^AB Z_MNEI
      Store target in W'(AB,EI)
  **/

  /** Z(MN,EI) <-- <MN||EF> t_I^F **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 2, 21, 2, 21, 0, "Z(MN,EI)");
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&D, &T1, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  /** tau_MN^AB Z(MN,EI) --> W'(AB,EI) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 2, 21, 2, 21, 0, "Z(MN,EI)");
  dpd_buf4_init(&W, CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_contract444(&T2, &Z, &W, 1, 1, 1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z);

  /**** Term VI ****/

  /** tau_MN^AB <MN||EI> --> W'(AB,EI) **/
  dpd_buf4_init(&W, CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_buf4_init(&E, CC_EINTS, 0, 2, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
  dpd_contract444(&T2, &E, &W, 1, 1, -1, 1);
  dpd_buf4_close(&E);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&W);

  /**** Term VII ****/

  /** WABEI <-- <BM||EF> t_IM^AF + <Bm|Ef> t_Im^Af - <AM||EF> t_IM^BF - <Am|Ef> t_Im^Bf
      Evaluate in six steps:
        (1) Sort <BM||EF> and <Bm|Ef> to F(BE,MF) and F(BE,mf) ordering.
        (2) Z(BE,IA) = F(BE,MF) T(IA,MF) + F(BE,mf) T(IA,mf)
        (3) Sort Z(BE,IA) --> Z'(EI,AB)
	(4) Sort Z'(EI,AB) --> Z''(EI,BA)
        (5) AXPY: Z'(EI,AB) = Z'(EI,AB) - Z''(EI,BA)
        (6) AXPY: W(EI,AB) <-- Z'(EI,AB)
      NB: The storage for the sorts is expensive and will eventually require out-of-core 
          codes.
  **/

  /** <BM||EF> --> F(BE,MF) **/
  dpd_buf4_init(&F, CC_FINTS, 0, 21, 5, 21, 5, 1, "F <AI|BC>");
  dpd_buf4_sort(&F, CC_FINTS, prqs, 5, 20, "F <AI||BC> (AB,IC)");
  dpd_buf4_close(&F);

  /** <Bm|Ef> --> (BE,mf) **/
  dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  dpd_buf4_sort(&F, CC_FINTS, prqs, 5, 30, "F <Ai|Bc> (AB,ic)");
  dpd_buf4_close(&F);

  /** <BM||EF> t_IM^AF --> Z(BE,IA) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(BE,IA)");
  dpd_buf4_init(&F, CC_FINTS, 0, 5, 20, 5, 20, 0, "F <AI||BC> (AB,IC)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  dpd_contract444(&F, &T2, &Z, 0, 0, -1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** <Bm|Ef> t_Im^Af --> Z(BE,IA) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(BE,IA)");
  dpd_buf4_init(&F, CC_FINTS, 0, 5, 30, 5, 30, 0, "F <Ai|Bc> (AB,ic)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  dpd_contract444(&F, &T2, &Z, 0, 0, -1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** Z(BE,IA) --> Z'(EI,AB) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(BE,IA)");
  dpd_buf4_sort(&Z, CC_TMP0, qrsp, 21, 5, "Z'(EI,AB)");
  dpd_buf4_close(&Z);

  /** Z'(EI,AB) --> Z''(EI,BA) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 21, 5, 21, 5, 0, "Z'(EI,AB)");
  dpd_buf4_sort(&Z, CC_TMP0, pqsr, 21, 5, "Z''(EI,BA)");
  dpd_buf4_close(&Z);

  /** Z'(EI,AB) = Z'(EI,AB) - Z''(EI,BA) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 21, 5, 21, 5, 0, "Z'(EI,AB)");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 21, 5, 21, 5, 0, "Z''(EI,BA)");
  dpd_buf4_axpy(&Z2, &Z1, -1);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&Z1);

  /** W(EI,AB) <-- Z'(EI,AB) **/
  dpd_buf4_init(&W, CC_HBAR, 0, 21, 5, 21, 7, 0, "WEIAB");
  dpd_buf4_init(&Z, CC_TMP0, 0, 21, 5, 21, 5, 0, "Z'(EI,AB)");
  dpd_buf4_axpy(&Z, &W, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /**** Terms VIII and IX ****/

  /** WABEI <-- -P(AB) t_M^A { <MB||EI> + t_IN^BF <MN||EF> + t_In^Bf <Mn|Ef> }
      Evaluate in two steps: 
         (1) Z_MBEI = <MB||EI> + t_IN^BF <MN||EF> + tIn^Bf <Mn|Ef>
         (2) WABEI <-- - t_M^A Z_MBEI + t_M^B Z_MAEI
      Store target in W'(AB,EI)
  **/

  /** Z(MB,EI) <-- <MB||EI> **/
  dpd_buf4_init(&C, CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
  dpd_buf4_copy(&C, CC_TMP0, "Z(MB,EI)");
  dpd_buf4_close(&C);
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  dpd_buf4_scm(&Z, -1);
  dpd_buf4_close(&Z);

  /** <MN||EF> t_IN^BF --> Z(ME,IB) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  /** <Mn|Ef> t_In^Bf --> Z(ME,IB) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  /** Z(ME,IB) --> Z(MB,EI) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  dpd_buf4_sort_axpy(&Z, CC_TMP0, psqr, 20, 21, "Z(MB,EI)", 1);
  dpd_buf4_close(&Z);

  /** Z(AB,EI) <-- -t_M^A Z(MB,EI) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  dpd_buf4_init(&Z1, CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z1, &Z, 0, 0, 0, -1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&Z);

  /** Z(AB,EI) --> Z'(BA,EI) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  dpd_buf4_sort(&Z, CC_TMP0, qprs, 5, 21, "Z'(BA,EI)");
  dpd_buf4_close(&Z);

  /** Z(AB,EI) = Z(AB,EI) - Z'(BA,EI) **/
  dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z'(BA,EI)");
  dpd_buf4_axpy(&Z2, &Z1, -1);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&Z1);

  /** Z(AB,EI) --> W'(AB,EI) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  dpd_buf4_init(&W, CC_TMP0, 0, 5, 21, 7, 21, 0, "W'(AB,EI)");
  dpd_buf4_axpy(&Z, &W, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z);

  /**** Combine accumulated W'(AB,EI) and W(EI,AB) terms into WEIAB ****/
  dpd_buf4_init(&W, CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  dpd_buf4_sort_axpy(&W, CC_HBAR, rspq, 21, 7, "WEIAB", 1);
  dpd_buf4_close(&W);
}

}} // namespace psi::cchbar
