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

/* WaBeI_UHF(): Computes all contributions to the aBeI spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (eI,aB) ordering and is referred to on disk as "WaBeI".
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
** For the aBeI spin case, we evaluate these contractions with two
** target orderings, (aB,eI) and (eI,aB), depending on the term.
** After all terms have been evaluated, the (aB,eI) terms are sorted
** into (eI,aB) ordering and both groups arer added together.
**
** TDC, June 2002
*/

void WaBeI_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /**** Term I ****/

  /** W(eI,aB) <--- <eI|aB> **/
  dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  dpd_buf4_copy(&F, CC_HBAR, "WeIaB");
  dpd_buf4_close(&F);

  /**** Term II ****/

  /** W(eI,aB) <--- - F_me t_mI^aB **/
  dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
  dpd_buf4_init(&W, CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  dpd_contract244(&Fme, &T2, &W, 0, 0, 0, -1, 1);
  dpd_buf4_close(&W);
  dpd_file2_close(&Fme);
  dpd_buf4_close(&T2);

  /**** Term III ****/  /** This will require special out-of-core code **/

  /** Z(Ie,Ba) <--- t_I^F <Fe|Ba> **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 24, 28, 24, 28, 0, "Z(Ie,Ba)");
  dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &B, &Z, 1, 0, 0, 1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&B);
  /** Z(Ie,Ba) --> W'(aB,eI) **/
  dpd_buf4_sort(&Z, CC_TMP0, srqp, 29, 25, "W'(aB,eI)");
  dpd_buf4_close(&Z);

  /**** Term IV ****/

  /** WaBeI <-- - t_M^B <Ma|Fe> t_I^F - t_m^a <mB|eF> t_I^F
      Evaluate in three steps: 
          (1) Z_MaeI = - <aM|eF> t_I^F [stored (aM,eI)]
	  (2) Z_mBeI = <mB|eF> t_I^F   [stored (mB,eI)]
          (3) WaBeI <-- t_M^B Z_MaeI - t_m^a Z_mBeI
       Store targets in:  W(eI,aB)  and   W'(aB,eI)  
  **/

  /** Z(aM,eI) <-- - <aM|eF> t_I^F **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&F, &T1, &Z, 3, 1, 0, -1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** t_M^B Z(aM,eI) --> W(eI,aB) **/
  dpd_buf4_init(&W, CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  dpd_buf4_init(&Z, CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /** Z(mB,eI) <-- <mB|eF> t_I^F **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
  dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&F, &T1, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** - t_m^a Z(mB,eI) --> W'(aB,eI) **/
  dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  dpd_buf4_init(&Z, CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /**** Term V ****/

  /** WaBeI <-- tau_mN^aB <mN|eF> t_I^F
      Evaluate in two steps:
         (1) Z_mNeI = <mN|eF> t_I^F
         (2) WaBeI <-- tau_mN^aB Z_mNeI
      Store target in W'(aB,eI)
  **/

  /** Z(mN,eI) <-- <mN|eF> t_I^F **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 23, 25, 23, 25, 0, "Z(mN,eI)");
  dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&D, &T1, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&D);

  /** tau_mN^aB Z(mN,eI) --> W'(aB,eI) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 23, 25, 23, 25, 0, "Z(mN,eI)");
  dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
  dpd_contract444(&T2, &Z, &W, 1, 1, 1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z);

  /**** Term VI ****/

  /** tau_mN^aB <mN|eI> --> W'(aB,eI) **/
  dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
  dpd_buf4_init(&E, CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
  dpd_contract444(&T2, &E, &W, 1, 1, 1, 1);
  dpd_buf4_close(&E);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&W);

  /**** Term VII ****/

  /** WaBeI <-- <Bm|Fe> t_Im^Fa - <am||ef> t_Im^Bf - <aM|eF> t_IM^BF
      Evaluate in five steps:
        (1) Sort <Bm|Fe> to F(Be,mF) ordering.  (** Note that we assume a sort has already 
            been done for <am||ef> and <aM|eF> in the WABEI and Wabei terms. **)
        (2) Z'(Be,Ia) = F(Be,mF) T(Ia,mF)
        (3) Sort Z'(Be,Ia) --> W(eI,aB)
        (4) Z''(ae,IB) = - F(ae,mf) T(IB,mf) - F(ae,MF) T(IB,MF)
	(5) Sort Z''(ae,IB) --> W(eI,aB)

      NB: The storage for the sorts is expensive and will eventually require out-of-core 
          codes.
  **/

  /** <Bm|Fe> --> F(Be,mF) **/
  dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  dpd_buf4_sort(&F, CC_FINTS, psqr, 28, 27, "F <Ai|Bc> (Ac,iB)");
  dpd_buf4_close(&F);

  /** <Bm|Fe> t_Im^Fa --> Z(Be,Ia) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 28, 24, 28, 24, 0, "Z(Be,Ia)");
  dpd_buf4_init(&F, CC_FINTS, 0, 28, 27, 28, 27, 0, "F <Ai|Bc> (Ac,iB)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
  dpd_contract444(&F, &T2, &Z, 0, 0, -1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** Z(Be,Ia) --> W(eI,aB) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 28, 24, 28, 24, 0, "Z(Be,Ia)");
  dpd_buf4_sort_axpy(&Z, CC_HBAR, qrsp, 25, 29, "WeIaB", 1);
  dpd_buf4_close(&Z);

  /** Z''(ae,IB) <-- - <am||ef> t_Im^Bf **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 15, 20, 15, 20, 0, "Z(ae,IB)");
  dpd_buf4_init(&F, CC_FINTS, 0, 15, 30, 15, 30, 0, "F <ai||bc> (ab,ic)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  dpd_contract444(&F, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** Z''(ae,IB) <-- -<aM|eF> t_IM^BF **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 15, 20, 15, 20, 0, "Z(ae,IB)");
  dpd_buf4_init(&F, CC_FINTS, 0, 15, 20, 15, 20, 0, "F <aI|bC> (ab,IC)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  dpd_contract444(&F, &T2, &Z, 0, 0, 1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** Z''(ai,IB) --> W(eI,aB) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 15, 20, 15, 20, 0, "Z(ae,IB)");
  dpd_buf4_sort_axpy(&Z, CC_HBAR, qrps, 25, 29, "WeIaB", 1);
  dpd_buf4_close(&Z);


  /**** Terms VIII and IX ****/

  /** WaBeI <-- - t_m^a { <mB|eI> + t_In^Bf <mn||ef> + t_IN^BF <mN|eF> }
                + t_M^B {-<Ma|Ie> + t_In^Fa <Mn|Fe> } 
      Evaluate in three steps: 
         (1) Z_mBeI =  <mB|eI> + t_In^Bf <mn||ef> + tIN^BF <mN|eF>  [stored (mB,eI)]
         (2) Z_MaeI = -<Ma|Ie> + t_In^Fa <Mn|Fe>                    [stored (aM,eI)]
         (3) WaBeI <-- - t_m^a Z_mBeI + t_M^B Z_MaeI
      Store targets in     W'(aB,eI) and  W(eI,aB)
  **/

  /** Z(mB,eI) <-- <mB|eI> **/
  dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
  dpd_buf4_copy(&D, CC_TMP0, "Z(mB,eI)");
  dpd_buf4_close(&D);

  /** <mn||ef> t_In^Bf --> Z(me,IB) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  /** <mN|eF> t_IN^BF --> Z(me,IB) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  /** Z(me,IB) --> Z(mB,eI) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  dpd_buf4_sort_axpy(&Z, CC_TMP0, psqr, 27, 25, "Z(mB,eI)", 1);
  dpd_buf4_close(&Z);

  /** W'(aB,eI) <-- - t_m^a Z(mB,eI) **/
  dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  dpd_buf4_init(&Z, CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /** Z(aM,eI) <-- - <Ma|Ie> **/
  dpd_buf4_init(&C, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
  dpd_buf4_sort(&C, CC_TMP0, qpsr, 25, 25, "Z(aM,eI)");
  dpd_buf4_close(&C);
  dpd_buf4_init(&Z, CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  dpd_buf4_scm(&Z, -1);
  dpd_buf4_close(&Z);

  /** Z(Me,Ia) <-- t_In^Fa <Mn|Fe> **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ia)");
  dpd_buf4_init(&D, CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  /** Z(Me,Ia) --> Z(aM,eI) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ia)");
  dpd_buf4_sort_axpy(&Z, CC_TMP0, spqr, 25, 25, "Z(aM,eI)", 1);
  dpd_buf4_close(&Z);

  /** W(eI,aB) <-- t_M^B Z_MaeI **/
  dpd_buf4_init(&W, CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  dpd_buf4_init(&Z, CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /**** Combine accumulated W'(aB,eI) and W(eI,aB) terms into WeIaB ****/
  dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  dpd_buf4_sort_axpy(&W, CC_HBAR, rspq, 25, 29, "WeIaB", 1);
  dpd_buf4_close(&W);
}

}} // namespace psi::cchbar
