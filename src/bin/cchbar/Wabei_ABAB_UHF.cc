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

/* WAbEi_UHF(): Computes all contributions to the AbEi spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (Ei,Ab) ordering and is referred to on disk as "WAbEi".
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
** For the AbEi spin case, we evaluate these contractions with two
** target orderings, (Ab,Ei) and (Ei,Ab), depending on the term.
** After all terms have been evaluated, the (Ab,Ei) terms are sorted
** into (Ei,Ab) ordering and both groups arer added together.
**
** TDC, June 2002
*/

void WAbEi_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /**** Term I ****/

  /** W(Ei,Ab) <--- <Ei|Ab> **/
  dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  dpd_buf4_copy(&F, CC_HBAR, "WEiAb");
  dpd_buf4_close(&F);

  /**** Term II ****/

  /** W(Ei,Ab) <--- - F_ME t_Mi^Ab **/
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  dpd_contract244(&Fme, &T2, &W, 0, 0, 0, -1, 1);
  dpd_buf4_close(&W);
  dpd_file2_close(&Fme);
  dpd_buf4_close(&T2);

  /**** Term III ****/

  /** <Ab|Ef> t_i^f **/
  dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&B, &T1, &W, 3, 1, 0, 1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&B);
  dpd_buf4_close(&W);

  /**** Term IV ****/

  /** WAbEi <-- - t_m^b <mA|fE> t_i^f - t_M^A <Mb|Ef> t_i^f
      Evaluate in three steps: 
          (1) Z_mAEi = - <Am|Ef> t_i^f [stored (Am,Ei)]
	  (2) Z_MbEi = <Mb|Ef> t_i^f   [stored (Mb,Ei)]
          (3) WAbEi <-- t_m^b Z_mAEi - t_M^A Z_MbEi
       Store targets in:  W(Ei,Ab)  and   W'(Ab,Ei)  
  **/

  /** Z(Am,Ei) <-- - <Am|Ef> t_i^f **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&F, &T1, &Z, 3, 1, 0, -1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** t_m^b Z(Am,Ei) --> W(Ei,Ab) **/
  dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  dpd_buf4_init(&Z, CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /** Z(Mb,Ei) <-- <Mb|Ef> t_i^f **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
  dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&F, &T1, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** - t_M^A Z(Mb,Ei) --> W'(Ab,Ei) **/
  dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  dpd_buf4_init(&Z, CC_TMP0, 0, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /**** Term V ****/

  /** WAbEi <-- tau_Mn^Ab <Mn|Ef> t_i^f
      Evaluate in two steps:
         (1) Z_MnEi = <Mn|Ef> t_i^f
         (2) WAbEi <-- tau_Mn^Ab Z_MnEi
      Store target in W'(Ab,Ei)
  **/

  /** Z(Mn,Ei) <-- <Mn|Ef> t_i^f **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 22, 26, 22, 26, 0, "Z(Mn,Ei)");
  dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&D, &T1, &Z, 3, 1, 0, 1, 0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&D);

  /** tau_Mn^Ab Z1(Mn,Ei) --> W'(Ab,Ei) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 22, 26, 22, 26, 0, "Z(Mn,Ei)");
  dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  dpd_contract444(&T2, &Z, &W, 1, 1, 1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z);

  /**** Term VI ****/

  /** tau_Mn^Ab <Mn|Ei> --> Z(Ab,Ei) **/
  dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
  dpd_contract444(&T2, &E, &W, 1, 1, 1, 1);
  dpd_buf4_close(&E);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&W);

  /**** Term VII ****/

  /** WAbEi <-- <bM|fE> t_iM^fA - <AM||EF> t_iM^bF - <Am|Ef> t_im^bf
      Evaluate in five steps:
        (1) Sort <bM|fE> to F(bE,Mf) ordering.  (** Note that we assume a sort has already 
            been done for <AM||EF> and <Am|Ef> in the WABEI and Wabei terms. **)
        (2) Z'(bE,iA) = F(bE,Mf) T(iA,Mf)
        (3) Sort Z'(bE,iA) --> W(Ei,Ab)
        (4) Z''(AE,ib) = - F(AE,MF) T(ib,MF) - F(AE,mf) T(ib,mf)
	(5) Sort Z''(AE,ib) --> W(Ei,Ab)

      NB: The storage for the sorts is expensive and will eventually require out-of-core 
          codes.
  **/

  /** <bM|fE> --> F(bE,Mf) **/
  dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  dpd_buf4_sort(&F, CC_FINTS, psqr, 29, 24, "F <aI|bC> (aC,Ib)");
  dpd_buf4_close(&F);

  /** <bM|fE> t_iM^fA --> Z(bE,iA) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 29, 27, 29, 27, 0, "Z(bE,iA)");
  dpd_buf4_init(&F, CC_FINTS, 0, 29, 24, 29, 24, 0, "F <aI|bC> (aC,Ib)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 27, 24, 27, 24, 0, "tiBJa");
  dpd_contract444(&F, &T2, &Z, 0, 0, -1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** Z(bE,iA) --> W(Ei,Ab) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 29, 27, 29, 27, 0, "Z(bE,iA)");
  dpd_buf4_sort_axpy(&Z, CC_HBAR, qrsp, 26, 28, "WEiAb", 1);
  dpd_buf4_close(&Z);

  /** Z''(AE,ib) <-- - <AM||EF> t_iM^bF **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 30, 5, 30, 0, "Z(AE,ib)");
  dpd_buf4_init(&F, CC_FINTS, 0, 5, 20, 5, 20, 0, "F <AI||BC> (AB,IC)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
  dpd_contract444(&F, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** Z''(AE,ib) <-- -<Am|Ef> t_im^bf **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 30, 5, 30, 0, "Z(AE,ib)");
  dpd_buf4_init(&F, CC_FINTS, 0, 5, 30, 5, 30, 0, "F <Ai|Bc> (AB,ic)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  dpd_contract444(&F, &T2, &Z, 0, 0, 1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&Z);

  /** Z''(AE,ib) --> W(Ei,Ab) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 30, 5, 30, 0, "Z(AE,ib)");
  dpd_buf4_sort_axpy(&Z, CC_HBAR, qrps, 26, 28, "WEiAb", 1);
  dpd_buf4_close(&Z);

  /**** Terms VIII and IX ****/

  /** WAbEi <-- - t_M^A { <Mb|Ei> + t_iN^bF <MN||EF> + t_in^bf <Mn|Ef> }
                + t_m^b {-<mA|iE> + t_iN^fA <mN|fE> } 
      Evaluate in three steps: 
         (1) Z_MbEi =  <Mb|Ei> + t_iN^bF <MN||EF> + tin^bf <Mn|Ef>  [stored (Mb,Ei)]
         (2) Z_mAEi = -<mA|iE> + t_iN^fA <mN|fE>                    [stored (Am,Ei)]
         (3) WAbEi <-- - t_M^A Z_MbEi + t_m^b Z_mAEi
      Store targets in     W'(Ab,Ei) and  W(Ei,AB)
  **/

  /** Z(Mb,Ei) <-- <Mb|Ei> **/
  dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
  dpd_buf4_copy(&D, CC_TMP0, "Z(Mb,Ei)");
  dpd_buf4_close(&D);

  /** <MN||EF> t_iN^bF --> Z(ME,ib) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,ib)");
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  /** <Mn|Ef> t_in^bf --> Z(ME,ib) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,ib)");
  dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  /** Z(ME,ib) --> Z(Mb,Ei) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,ib)");
  dpd_buf4_sort_axpy(&Z, CC_TMP0, psqr, 24, 26, "Z(Mb,Ei)", 1);
  dpd_buf4_close(&Z);

  /** W'(Ab,Ei) <-- - t_M^A Z(Mb,Ei) **/
  dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  dpd_buf4_init(&Z, CC_TMP0, 0, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /** Z(Am,Ei) <-- - <mA|iE> **/
  dpd_buf4_init(&C, CC_CINTS, 0, 26, 26, 26, 26, 0, "C <Ai|Bj>");
  dpd_buf4_copy(&C, CC_TMP0, "Z(Am,Ei)");
  dpd_buf4_close(&C);
  dpd_buf4_init(&Z, CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  dpd_buf4_scm(&Z, -1);
  dpd_buf4_close(&Z);

  /** Z(mE,iA) <-- t_iN^fA <mN|fE> **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(mE,iA)");
  dpd_buf4_init(&D, CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 27, 24, 27, 24, 0, "tiBJa");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_close(&Z);

  /** Z(mE,iA) --> Z(Am,Ei) **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(mE,iA)");
  dpd_buf4_sort_axpy(&Z, CC_TMP0, spqr, 26, 26, "Z(Am,Ei)", 1);
  dpd_buf4_close(&Z);

  /** W(Ei,AB) <-- t_m^b Z_mAEi **/
  dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
  dpd_buf4_init(&Z, CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /**** Combine accumulated W'(Ab,Ei) and W(Ei,Ab) terms into WEiAb ****/
  dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "W'(Ab,Ei)");
  dpd_buf4_sort_axpy(&W, CC_HBAR, rspq, 26, 28, "WEiAb", 1);
  dpd_buf4_close(&W);
}

}} // namespace psi::cchbar
