/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

/* WmbejL2(): Computes the contributions of the Wmbej HBAR matrix
** elements to the Lambda double deexcitation amplitude equations.
** These contributions are written in spin orbitals as:
**
** L_ij^ab <-- P(ij) P(ab) L_im^ae Wjebm =
**     L_im^ae Wjebm - L_jm^ae Wiebm - L_im^be Wjeam + L_jm^be Wieam
**
** The matrix elements for all six spin cases are stored in (me,jb)
** ordering.  This leads to the following contractions for the three
** L2 spin cases:
** 
** L(IA,JB) <-- L(IA,ME) W(JB,ME) + L(IA,me) W(JB,me) 
**            - L(JA,ME) W(IB,ME) - L(JA,me) W(IB,me)
**            - L(IB,ME) W(JA,ME) - L(IB,me) W(JA,me)
**            + L(JB,ME) W(IA,ME) + L(JB,me) W(IA,me)
**    (only two unique contractions)
** L(ia,jb) <-- L(ia,me) W(jb,me) + L(ia,ME) W(jb,ME) 
**            - L(ja,me) W(ib,me) - L(ja,ME) W(ib,ME)
**            - L(ib,me) W(ja,me) - L(ib,ME) W(ja,ME)
**            + L(jb,me) W(ia,me) + L(jb,ME) W(ia,ME)
**    (only two unique contractions)
** L(IA,jb) <-- L(IA,ME) W(jb,ME) + L(IA,me) W(jb,me)
**            - L(jA,Me) W(Ib,Me)
**            - L(Ib,Me) W(jA,Me)
**            + L(jb,ME) W(IA,ME) + L(jb,me) W(IA,me)
**    (all six contractions are unique for ROHF and UHF orbitals)
**
** TDC, July 2002
**
** For RHF contractions, we evaluate the contractions as follows:
**
** + L(jb,ME) W(IA,ME) + L(jb,me) W(IA,me)                I
** + L(jA,Me) W(Ib,Me) + L(Ib,mE) W(jA,mE)           II   +   III
** + L(IA,ME) W(jb,ME) + L(IA,me) W(jb,me)               IV
**
** Similar to what we did in ccenergy/WmbejT2.c, the AB L2 terms labelled I and IV
** above may be written as (apart from index swapping):
**
** 1/2 [2 L(jb,ME) - L(jB,Me)] [2 W(me,IA) + W(Me,Ia)] + 1/2 L(jB,Me) W(Me,Ia)
**
** Terms II and III are exactly the same as the last term above, apart from the
** factor of 1/2 and index swapping.  So, for RHF orbitals, we need to evaluate two
** contractions.
**
** TDC, March 2004
** 
*/

void WmbejL2(int L_irr)
{
  dpdbuf4 newL2, L2, W, Z, Z2;

  /* RHS += P(ij)P(ab)Limae * Wjebm */
  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 10, 10, 10, 10, 0, "Z(Ib,jA)");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIbjA");
    dpd_contract444(&W, &L2, &Z, 0, 1, 1, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);
    dpd_buf4_sort(&Z, CC_TMP0, psrq, 10, 10, "Z(IA,jb) III");
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 10, 10, 10, 10, 0, "Z(IA,jb) I");

    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "2 LIAjb - LIbjA");
    dpd_contract444(&W, &L2, &Z, 0, 1, 0.5, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Z2, CC_TMP0, L_irr, 10, 10, 10, 10, 0, "Z(Ib,jA)");
    dpd_buf4_axpy(&Z2, &Z, 0.5);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z2, CC_TMP0, L_irr, 10, 10, 10, 10, 0, "Z(IA,jb) III");
    dpd_buf4_axpy(&Z2, &Z, 1);
    dpd_buf4_close(&Z2);

    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, prqs, 0, 5, "New LIjAb", 1);
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, rpsq, 0, 5, "New LIjAb", 1);
    dpd_buf4_close(&Z);

  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 10, 10, 10, 10, 0, "Z(IA,JB)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_sort(&Z, CC_TMP1, rqps, 10, 10, "Z(JA,IB)");
    dpd_buf4_sort(&Z, CC_TMP2, psrq, 10, 10, "Z(IB,JA)");
    dpd_buf4_sort(&Z, CC_TMP3, rspq, 10, 10, "Z(JB,IA)");
    dpd_buf4_init(&Z2, CC_TMP1, L_irr, 10, 10, 10, 10, 0, "Z(JA,IB)");
    dpd_buf4_axpy(&Z2, &Z, -1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z2, CC_TMP2, L_irr, 10, 10, 10, 10, 0, "Z(IB,JA)");
    dpd_buf4_axpy(&Z2, &Z, -1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z2, CC_TMP3, L_irr, 10, 10, 10, 10, 0, "Z(JB,IA)");
    dpd_buf4_axpy(&Z2, &Z, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_sort(&Z, CC_TMP1, prqs, 0, 5, "Z(IJ,AB)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 0, 5, 0, 5, 0, "Z(IJ,AB)");
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&Z, &newL2, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&newL2);

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 10, 10, 10, 10, 0, "Z(ia,jb)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "Liajb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LiaJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_sort(&Z, CC_TMP1, rqps, 10, 10, "Z(ja,ib)");
    dpd_buf4_sort(&Z, CC_TMP2, psrq, 10, 10, "Z(ib,ja)");
    dpd_buf4_sort(&Z, CC_TMP3, rspq, 10, 10, "Z(jb,ia)");
    dpd_buf4_init(&Z2, CC_TMP1, L_irr, 10, 10, 10, 10, 0, "Z(ja,ib)");
    dpd_buf4_axpy(&Z2, &Z, -1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z2, CC_TMP2, L_irr, 10, 10, 10, 10, 0, "Z(ib,ja)");
    dpd_buf4_axpy(&Z2, &Z, -1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z2, CC_TMP3, L_irr, 10, 10, 10, 10, 0, "Z(jb,ia)");
    dpd_buf4_axpy(&Z2, &Z, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_sort(&Z, CC_TMP1, prqs, 0, 5, "Z(ij,ab)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 0, 5, 0, 5, 0, "Z(ij,ab)");
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "New Lijab");
    dpd_buf4_axpy(&Z, &newL2, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&newL2);


    dpd_buf4_init(&Z, CC_TMP0, L_irr, 10, 10, 10, 10, 0, "Z(IA,jb)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "Liajb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
    dpd_contract444(&W, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LiaJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
    dpd_contract444(&W, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_sort(&Z, CC_TMP1, prqs, 0, 5, "Z(Ij,Ab)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_axpy(&Z, &newL2, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&newL2);

    dpd_buf4_init(&Z, CC_TMP0, L_irr, 10, 10, 10, 10, 0, "Z(Ib,jA)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LIbjA");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    dpd_contract444(&W, &L2, &Z, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 10, 10, 10, 10, 0, "LjAIb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBEj");
    dpd_contract444(&L2, &W, &Z, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);
    dpd_buf4_sort(&Z, CC_TMP1, prqs, 0, 5, "Z(Ij,bA)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP1, L_irr, 0, 5, 0, 5, 0, "Z(Ij,bA)");
    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 0, 5, "Z(Ij,Ab)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, L_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_axpy(&Z, &newL2, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&newL2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&Z, CC_TMP2, L_irr, 20, 20, 20, 20, 0, "Z(IA,JB)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_sort(&Z, CC_TMP2, rqps, 20, 20, "Z(JA,IB)");
    dpd_buf4_sort(&Z, CC_TMP2, psrq, 20, 20, "Z(IB,JA)");
    dpd_buf4_sort(&Z, CC_TMP2, rspq, 20, 20, "Z(JB,IA)");
    dpd_buf4_init(&Z2, CC_TMP2, L_irr, 20, 20, 20, 20, 0, "Z(JA,IB)");
    dpd_buf4_axpy(&Z2, &Z, -1);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z2, CC_TMP2, L_irr, 20, 20, 20, 20, 0, "Z(IB,JA)");
    dpd_buf4_axpy(&Z2, &Z, -1);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z2, CC_TMP2, L_irr, 20, 20, 20, 20, 0, "Z(JB,IA)");
    dpd_buf4_axpy(&Z2, &Z, 1);
    dpd_buf4_close(&Z2);
    dpd_buf4_sort(&Z, CC_TMP2, prqs, 0, 5, "Z(IJ,AB)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP2, L_irr, 0, 5, 0, 5, 0, "Z(IJ,AB)");
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "New LIJAB");
    dpd_buf4_axpy(&Z, &newL2, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&newL2);

    dpd_buf4_init(&Z, CC_TMP2, L_irr, 30, 30, 30, 30, 0, "Z(ia,jb)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 30, 30, 30, 30, 0, "Liajb");
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 30, 20, 30, 20, 0, "LiaJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_sort(&Z, CC_TMP2, rqps, 30, 30, "Z(ja,ib)");
    dpd_buf4_sort(&Z, CC_TMP2, psrq, 30, 30, "Z(ib,ja)");
    dpd_buf4_sort(&Z, CC_TMP2, rspq, 30, 30, "Z(jb,ia)");
    dpd_buf4_init(&Z2, CC_TMP2, L_irr, 30, 30, 30, 30, 0, "Z(ja,ib)");
    dpd_buf4_axpy(&Z2, &Z, -1);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z2, CC_TMP2, L_irr, 30, 30, 30, 30, 0, "Z(ib,ja)");
    dpd_buf4_axpy(&Z2, &Z, -1);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z2, CC_TMP2, L_irr, 30, 30, 30, 30, 0, "Z(jb,ia)");
    dpd_buf4_axpy(&Z2, &Z, 1);
    dpd_buf4_close(&Z2);
    dpd_buf4_sort(&Z, CC_TMP2, prqs, 10, 15, "Z(ij,ab)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP2, L_irr, 10, 15, 10, 15, 0, "Z(ij,ab)");
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 10, 15, 12, 17, 0, "New Lijab");
    dpd_buf4_axpy(&Z, &newL2, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&newL2);

    dpd_buf4_init(&Z, CC_TMP2, L_irr, 20, 30, 20, 30, 0, "Z(IA,jb)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    dpd_contract444(&L2, &W, &Z, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 30, 30, 30, 30, 0, "Liajb");
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    dpd_contract444(&W, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 30, 20, 30, 20, 0, "LiaJB");
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    dpd_contract444(&W, &L2, &Z, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, prqs, 22, 28, "New LIjAb", 1);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP2, L_irr, 24, 27, 24, 27, 0, "Z(Ib,jA)");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 24, 27, 24, 27, 0, "LIbjA");
    dpd_buf4_init(&W, CC_HBAR, 0, 24, 24, 24, 24, 0, "WMbeJ");
    dpd_contract444(&W, &L2, &Z, 0, 1, 1, 0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 27, 24, 27, 24, 0, "LjAIb");
    dpd_buf4_init(&W, CC_HBAR, 0, 27, 27, 27, 27, 0, "WmBEj");
    dpd_contract444(&L2, &W, &Z, 1, 0, 1, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC_LAMBDA, prsq, 22, 28, "New LIjAb", 1);
    dpd_buf4_close(&Z);

  }
}

}} // namespace psi::cclambda
