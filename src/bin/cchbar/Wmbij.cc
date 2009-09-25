/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <math.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/** Wmbij_build(): Constructs the Wmbij HBAR intermediate, defined in 
 ** spin orbitals as:
 **
 ** Wmbij = <mb||ij> - Fme t_ij^be - t_n^b Wmnij + 1/2 <mb||ef> tau_ij^ef
 **    + P(ij) <mn||ie> t_jn^be + P(ij) t_i^e { <mb||ej> - t_nj^bf <mn||ef> }
 **
 ** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995)]
 **
 ** For RHF orbitals, there is only one unique spin case: WMbIj.
 **
 ** TDC, March 2004
 */

void Wmbij_build(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 W, E, T2, Wmnij, I, Tau, Z, Z1, Z2, C, D;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_sort(&E, CC_HBAR, rspq, 10, 0, "WMbIj");
    dpd_buf4_close(&E);

  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    /** <MB||IJ> **/
    dpd_buf4_sort(&E, CC_HBAR, rspq, 10, 2, "WMBIJ");
    /** <mb||ij> **/
    dpd_buf4_sort(&E, CC_HBAR, rspq, 10, 2, "Wmbij");
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    /** <Mb|Ij> **/
    dpd_buf4_sort(&E, CC_HBAR, rspq, 10, 0, "WMbIj");
    /** <mB|iJ> **/
    dpd_buf4_sort(&E, CC_HBAR, rspq, 10, 0, "WmBiJ");
    dpd_buf4_close(&E);

  }
  else if(params.ref == 2) { /** UHF **/

    /** <MB||IJ> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_buf4_sort(&E, CC_HBAR, rspq, 20, 2, "WMBIJ");
    dpd_buf4_close(&E);

    /** <mb||ij> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&E, CC_HBAR, rspq, 30, 12, "Wmbij");
    dpd_buf4_close(&E);

    /** <Mb|Ij> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_sort(&E, CC_HBAR, rspq, 24, 22, "WMbIj");
    dpd_buf4_close(&E);

    /** <mB|iJ> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_buf4_sort(&E, CC_HBAR, rspq, 27, 23, "WmBiJ");
    dpd_buf4_close(&E);

  }

  if(params.ref == 0) { /** RHF **/

    /** F_ME t_Ij^Eb --> W(Mb,Ij) **/
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    dpd_file2_close(&Fme);

  }
  else if(params.ref ==1) { /** ROHF **/

    /** F_ME t_IJ^EB --> W(MB,IJ) **/
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
    dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    dpd_file2_close(&Fme);

    /** F_me t_ij^eb --> W(mb,ij) **/
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
    dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    dpd_file2_close(&Fme);

    /** F_ME t_Ij^Eb --> W(Mb,Ij) **/
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    dpd_file2_close(&Fme);

    /** F_me t_iJ^eB --> W(mB,iJ) **/
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    dpd_file2_close(&Fme);

  }
  else if(params.ref == 2) { /** UHF **/

    /** F_ME t_IJ^EB --> W(MB,IJ) **/
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 2, 20, 2, 0, "WMBIJ");
    dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    dpd_file2_close(&Fme);

    /** F_me t_ij^eb --> W(mb,ij) **/
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 12, 30, 12, 0, "Wmbij");
    dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    dpd_file2_close(&Fme);

    /** F_ME t_Ij^Eb --> W(Mb,Ij) **/
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&W, CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    dpd_file2_close(&Fme);

    /** F_me t_iJ^eB --> W(mB,iJ) **/
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&W, CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    dpd_contract244(&Fme, &T2, &W, 1, 2, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&T2);
    dpd_file2_close(&Fme);

  }

  if(params.ref == 0) { /** RHF **/

    /** - t_n^b W_MnIj **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);
    dpd_file2_close(&T1);

  }
  else if(params.ref == 1) { /** ROHF **/

    /** - t_N^B W_MNIJ --> W(MB,IJ) **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 0, 2, 2, 2, 0, "WMNIJ");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
    dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);
    dpd_file2_close(&T1);

    /** - t_n^b W_mnij **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 0, 2, 2, 2, 0, "Wmnij");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
    dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);
    dpd_file2_close(&T1);

    /** - t_n^b W_MnIj **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);
    dpd_file2_close(&T1);

    /** - t_N^B W_mNiJ **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_buf4_sort(&Wmnij, CC_TMP0, qprs, 0, 0, "WnMIj");
    dpd_buf4_close(&Wmnij);
    dpd_buf4_init(&Wmnij, CC_TMP0, 0, 0, 0, 0, 0, 0, "WnMIj");
    dpd_buf4_sort(&Wmnij, CC_TMP1, pqsr, 0, 0, "WnMjI");
    dpd_buf4_close(&Wmnij);
    dpd_buf4_init(&Wmnij, CC_TMP1, 0, 0, 0, 0, 0, 0, "WnMjI");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);
    dpd_file2_close(&T1);
  }
  else if(params.ref == 2) { /** UHF **/

    /** - t_N^B W_MNIJ --> W(MB,IJ) **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 0, 2, 2, 2, 0, "WMNIJ");
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 2, 20, 2, 0, "WMBIJ");
    dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);
    dpd_file2_close(&T1);

    /** - t_n^b W_mnij **/
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 10, 12, 12, 12, 0, "Wmnij");
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 12, 30, 12, 0, "Wmbij");
    dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);
    dpd_file2_close(&T1);

    /** - t_n^b W_MnIj **/
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    dpd_buf4_init(&W, CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);
    dpd_file2_close(&T1);

    /** - t_N^B W_mNiJ **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Wmnij, CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");
    dpd_buf4_sort(&Wmnij, CC_TMP0, qpsr, 23, 23, "WmNiJ");
    dpd_buf4_close(&Wmnij);
    dpd_buf4_init(&Wmnij, CC_TMP0, 0, 23, 23, 23, 23, 0, "WmNiJ");
    dpd_buf4_init(&W, CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    dpd_contract424(&Wmnij, &T1, &W, 1, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Wmnij);
    dpd_file2_close(&T1);

  }

  if(params.ref == 0) { /** RHF **/

    /** <Mb|Ef> tau_Ij^Ef **/
    dpd_buf4_init(&I, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0); /* should run OCC, if needed */
    dpd_buf4_close(&W);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&I);

  }
  else if(params.ref == 1) { /** ROHF **/

    /** <MB||EF> tau_IJ^EF **/
    dpd_buf4_init(&I, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
    dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&I);

    /* <mb||ef> tau_ij^ef **/
    dpd_buf4_init(&I, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
    dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&I);

    /** <Mb|Ef> tau_Ij^Ef **/
    dpd_buf4_init(&I, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&I);

    /** <mB|eF> tau_iJ^eF **/
    dpd_buf4_init(&I, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&I);
  }
  else if(params.ref == 2) { /** UHF **/

    /** <MB||EF> tau_IJ^EF **/
    dpd_buf4_init(&I, CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 2, 20, 2, 0, "WMBIJ");
    dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&I);
    dpd_buf4_close(&W);

    /* <mb||ef> tau_ij^ef **/
    dpd_buf4_init(&I, CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 12, 30, 12, 0, "Wmbij");
    dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&I);

    /** <Mb|Ef> tau_Ij^Ef **/
    dpd_buf4_init(&I, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_buf4_init(&W, CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&I);

    /** <mB|eF> tau_iJ^eF **/
    dpd_buf4_init(&I, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
    dpd_buf4_init(&W, CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    dpd_contract444(&I, &Tau, &W, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&I);

  }

  /* Sort <ij||ka> integrals for the E*T2 contributions */
  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&I, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_sort(&I, CC_TMP0, prqs, 0, 10, "<Mn|Ie> (MI,ne)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    dpd_buf4_sort(&I, CC_TMP0, sqrp, 0, 10, "2 <Mn|Ie> - <Nm|Ie> (MI,ne)");
    dpd_buf4_close(&I);

  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_buf4_init(&I, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&I, CC_TMP0, prqs, 0, 10, "I(MI,NE)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_sort(&I, CC_TMP1, prqs, 0, 10, "I(MI,NE)");
    dpd_buf4_close(&I);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&I, CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_buf4_sort(&I, CC_TMP0, prqs, 0, 20, "<MN||IE> (MI,NE)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&I, CC_TMP0, prqs, 10, 30, "<mn||ie> (mi,ne)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_sort(&I, CC_TMP0, prqs, 0, 30, "<Mn|Ie> (MI,ne)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_buf4_sort(&I, CC_TMP0, prqs, 10, 20, "<mN|iE> (mi,NE)");
    dpd_buf4_close(&I);

  }

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(MI,jb)");

    dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "2 <Mn|Ie> - <Nm|Ie> (MI,ne)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "<Mn|Ie> (MI,ne)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_contract444(&I, &T2, &Z, 0, 1, -1, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    dpd_buf4_sort_axpy(&Z, CC_HBAR, psqr, 10, 0, "WMbIj", 1);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(Mj,Ib)");
    dpd_buf4_init(&I, CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 11, 10, 11, 0, "tIbAj");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    dpd_buf4_sort_axpy(&Z, CC_HBAR, psrq, 10, 0, "WMbIj", -1);
    dpd_buf4_close(&Z);

  }
  else if(params.ref == 1) { /** ROHF **/

    /** <MN||IE> t_JN^BE **/
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,JB)");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** <Mn||Ie> t_Jn^Be **/
    dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** <MN||JE> t_IN^BE **/
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_buf4_init(&Z, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
    dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** <Mn||Je> t_In^Be **/
    dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    dpd_buf4_init(&Z1, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,JB)");
    dpd_buf4_sort(&Z1, CC_TMP4, prqs, 0, 10, "Z(MJ,IB)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP4, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
    dpd_buf4_init(&Z2, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
    dpd_buf4_axpy(&Z1, &Z2, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z2, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(MJ,IB)");
    dpd_buf4_sort(&Z2, CC_TMP4, psrq, 10, 0, "Z(MB,IJ)");
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(MB,IJ)");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 2, 0, "WMBIJ");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /** <mn||ie> t_jn^be **/
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,jb)");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** <mN||iE> t_jN^bE **/
    dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** <mn||je> t_in^be **/
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_buf4_init(&Z, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
    dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** <mN||jE> t_iN^bE **/
    dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    dpd_buf4_init(&Z1, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,jb)");
    dpd_buf4_sort(&Z1, CC_TMP4, prqs, 0, 10, "Z(mj,ib)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z1, CC_TMP4, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
    dpd_buf4_init(&Z2, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
    dpd_buf4_axpy(&Z1, &Z2, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z2, CC_TMP3, 0, 0, 10, 0, 10, 0, "Z(mj,ib)");
    dpd_buf4_sort(&Z2, CC_TMP4, psrq, 10, 0, "Z(mb,ij)");
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(mb,ij)");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 2, 0, "Wmbij");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /** <MN||IE> t_jN^bE **/
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,jb)");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_sort(&Z, CC_TMP3, psrq, 10, 0, "Z(Mb,jI)");
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);
    dpd_buf4_init(&Z, CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(Mb,jI)");
    dpd_buf4_sort(&Z, CC_TMP4, pqsr, 10, 0, "Z(Mb,Ij)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /** <Mn|Ie> t_jn^be **/
    dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(MI,jb)");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_sort(&Z, CC_TMP3, psrq, 10, 0, "Z(Mb,jI)");
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);
    dpd_buf4_init(&Z, CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(Mb,jI)");
    dpd_buf4_sort(&Z, CC_TMP4, pqsr, 10, 0, "Z(Mb,Ij)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /** <mn||ie> t_Jn^Be **/
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,JB)");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_sort(&Z, CC_TMP3, psrq, 10, 0, "Z(mB,Ji)");
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);
    dpd_buf4_init(&Z, CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(mB,Ji)");
    dpd_buf4_sort(&Z, CC_TMP4, pqsr, 10, 0, "Z(mB,iJ)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(mB,iJ)");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /** <mN|iE> t_JN^BE **/
    dpd_buf4_init(&I, CC_TMP1, 0, 0, 10, 0, 10, 0, "I(MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mi,JB)");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_sort(&Z, CC_TMP3, psrq, 10, 0, "Z(mB,Ji)");
    dpd_buf4_close(&Z);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);
    dpd_buf4_init(&Z, CC_TMP3, 0, 10, 0, 10, 0, 0, "Z(mB,Ji)");
    dpd_buf4_sort(&Z, CC_TMP4, pqsr, 10, 0, "Z(mB,iJ)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP4, 0, 10, 0, 10, 0, 0, "Z(mB,iJ)");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);
  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 20, 0, 20, 0, "Z(MI,JB)");

    /** <MN||IE> t_JN^BE **/
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 20, 0, 20, 0, "<MN||IE> (MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** <Mn|Ie> t_Jn^Be **/
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 30, 0, 30, 0, "<Mn|Ie> (MI,ne)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** Z(MI,JB) --> Z(MB,IJ) **/
    dpd_buf4_sort(&Z, CC_TMP0, psqr, 20, 0, "Z(MB,IJ)");
    dpd_buf4_close(&Z);

    /** Z(MB,IJ) = Z(MB,IJ) - Z(MB,JI) **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 20, 0, 20, 0, 0, "Z(MB,IJ)");
    dpd_buf4_sort(&Z1, CC_TMP0, pqsr, 20, 0, "Z(MB,JI)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 20, 0, 20, 0, 0, "Z(MB,JI)");
    dpd_buf4_axpy(&Z2, &Z1, -1.0);
    dpd_buf4_close(&Z2);

    /** Z(MB,IJ) --> W(MB,IJ) **/
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 0, 20, 2, 0, "WMBIJ");
    dpd_buf4_axpy(&Z1, &W, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&W);


    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 30, 10, 30, 0, "Z(mi,jb)");

    /** <mn||ie> t_jn^be **/
    dpd_buf4_init(&I, CC_TMP0, 0, 10, 30, 10, 30, 0, "<mn||ie> (mi,ne)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** <mN||iE> t_jN^bE **/
    dpd_buf4_init(&I, CC_TMP0, 0, 10, 20, 10, 20, 0, "<mN|iE> (mi,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** Z(mi,jb) --> Z(mb,ij) **/
    dpd_buf4_sort(&Z, CC_TMP0, psqr, 30, 10, "Z(mb,ij)");
    dpd_buf4_close(&Z);

    /** Z(mb,ij) = Z(mb,ij) - Z(mb,ji) **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 30, 10, 30, 10, 0, "Z(mb,ij)");
    dpd_buf4_sort(&Z1, CC_TMP0, pqsr, 30, 10, "Z(mb,ji)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 30, 10, 30, 10, 0, "Z(mb,ji)");
    dpd_buf4_axpy(&Z2, &Z1, -1.0);
    dpd_buf4_close(&Z2);

    /** Z(mb,ij) --> W(mb,ij) **/
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 10, 30, 12, 0, "Wmbij");
    dpd_buf4_axpy(&Z1, &W, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&W);


    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 30, 0, 30, 0, "Z(MI,jb)");

    /** <MN||IE> t_jN^bE **/
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 20, 0, 20, 0, "<MN||IE> (MI,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** <Mn|Ie> t_jn^be **/
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 30, 0, 30, 0, "<Mn|Ie> (MI,ne)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** Z(MI,jb) --> Z(Mb,Ij) **/
    dpd_buf4_sort(&Z, CC_TMP0, psqr, 24, 22, "Z(Mb,Ij)");
    dpd_buf4_close(&Z);


    /** -<Mn|Ej> t_In^Eb **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 22, 24, 22, 24, 0, "Z(Mj,Ib)");
    dpd_buf4_init(&I, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 24, 26, 24, 26, 0, "tIbAj");
    dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    dpd_buf4_sort(&Z, CC_TMP0, psrq, 24, 22, "Z1(Mb,Ij)");
    dpd_buf4_close(&Z);

    /** Z(Mb,Ij) --> W(Mb,Ij) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 24, 22, 24, 22, 0, "Z(Mb,Ij)");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 24, 22, 24, 22, 0, "Z1(Mb,Ij)");
    dpd_buf4_axpy(&Z1, &Z, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&W, CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);


    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 20, 10, 20, 0, "Z(mi,JB)");

    /** <mn||ie> t_Jn^Be **/
    dpd_buf4_init(&I, CC_TMP0, 0, 10, 30, 10, 30, 0, "<mn||ie> (mi,ne)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** <mN|iE> t_JN^BE **/
    dpd_buf4_init(&I, CC_TMP0, 0, 10, 20, 10, 20, 0, "<mN|iE> (mi,NE)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_contract444(&I, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    /** Z(mi,JB) --> Z(mB,iJ) **/
    dpd_buf4_sort(&Z, CC_TMP0, psqr, 27, 23, "Z(mB,iJ)");
    dpd_buf4_close(&Z);

    /** Z(mB,iJ) --> W(mB,iJ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 27, 23, 27, 23, 0, "Z(mB,iJ)");
    dpd_buf4_init(&W, CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&W);

    dpd_buf4_close(&Z);

    /** -<mN|eJ> t_iN^eB **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 23, 27, 23, 27, 0, "Z(mJ,iB)");
    dpd_buf4_init(&I, CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 27, 25, 27, 25, 0, "tiBaJ");
    dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);

    dpd_buf4_sort(&Z, CC_TMP0, psrq, 27, 23, "Z(mB,iJ)");
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, 0, 27, 23, 27, 23, 0, "Z(mB,iJ)");
    dpd_buf4_init(&W, CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

  } /** UHF **/


  if(params.ref == 1) { /** RHF or ROHF **/

    /* Sort <ai||jk> integrals for remaining E*T2 contributions */
    /** THIS IS NOT ACTUALLY NECESSARY!  FIX THESE CONTRACTIONS! (11/14/01) **/
    dpd_buf4_init(&I, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_buf4_sort(&I, CC_TMP0, rspq, 0, 11, "I(Mn,Ej)");
    dpd_buf4_close(&I);
    dpd_buf4_init(&I, CC_TMP0, 0, 0, 11, 0, 11, 0, "I(Mn,Ej)");
    dpd_buf4_sort(&I, CC_TMP1, psrq, 0, 11, "I(Mj,En)");
    dpd_buf4_close(&I);


    /** -<Mn|Ej> t_In^EB **/ /** Can be written as: - <Mj|En> t_In^Eb !! **/
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&T2, CC_TMP0, psrq, 10, 11, "T2(Ib,En)");
    dpd_buf4_close(&T2);
    dpd_buf4_init(&I, CC_TMP1, 0, 0, 11, 0, 11, 0, "I(Mj,En)");
    dpd_buf4_init(&T2, CC_TMP0, 0, 10, 11, 10, 11, 0, "T2(Ib,En)");
    dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(Mj,Ib)");
    dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    dpd_buf4_sort(&Z, CC_TMP0, psrq, 10, 0, "Z(Mb,Ij)");
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);
  

    /** -<mN|eJ> t_iN^eB **/ /** Can be written as: -<mJ|eN> t_iN^Eb !! **/
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_sort(&T2, CC_TMP0, psrq, 10, 11, "T2(iB,eN)");
    dpd_buf4_close(&T2);
    dpd_buf4_init(&I, CC_TMP1, 0, 0, 11, 0, 11, 0, "I(Mj,En)");
    dpd_buf4_init(&T2, CC_TMP0, 0, 10, 11, 10, 11, 0, "T2(iB,eN)");
    dpd_buf4_init(&Z, CC_TMP2, 0, 0, 10, 0, 10, 0, "Z(mJ,iB)");
    dpd_contract444(&I, &T2, &Z, 0, 0, -1.0, 0.0);
    dpd_buf4_sort(&Z, CC_TMP0, psrq, 10, 0, "Z(mB,iJ)");
    dpd_buf4_close(&T2);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(mB,iJ)");
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    dpd_buf4_axpy(&Z, &W, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);
  }

  /** Prepare intermediates for final term of Wmbij **/

  if(params.ref == 0) { /** RHF **/

    /* Z(ME,jb) = { <Mb|Ej> + t_jN^bF [2 <Mn|Ef> - <Mn|Fe>] - t_jN^Fb <Mn|Ef> } */
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ME,jb)");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_contract444(&D, &T2, &Z, 0, 0, -1, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_axpy(&D, &Z, 1);
    dpd_buf4_close(&D);

    /* W(Mb,Ij) <-- Z(ME,jb) t_I^E */
    dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(MI,jb)");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract424(&Z, &T1, &Z1, 1, 1, 1, 1, 0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z);

    dpd_buf4_sort_axpy(&Z1, CC_HBAR, psqr, 10, 0, "WMbIj", 1);
    dpd_buf4_close(&Z1);


    /* Z(Me,Ib) = { - <Mb|Ie> + t_In^Fb <Mn|Fe> } */
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Me,Ib)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&C, CC_CINTS, 0, 11, 10, 11, 10, 0, "C <ia|jb> (bi,ja)");
    dpd_buf4_sort_axpy(&C, CC_TMP0, qprs, 10, 10, "Z(Me,Ib)", -1);
    dpd_buf4_close(&C);

    /* W(Mb,Ij) <-- Z(Me,Ib) t_j^e */
    dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 10, 0, 10, 0, "Z(Mj,Ib)");
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Me,Ib)");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract424(&Z, &T1, &Z1, 1, 1, 1, 1, 0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z);

    dpd_buf4_sort_axpy(&Z1, CC_HBAR, psrq, 10, 0, "WMbIj", -1);
    dpd_buf4_close(&Z1);
  }
  else if(params.ref == 1) { /** ROHF **/

    /** t_JN^BF <MN||EF> --> Z_MBJE **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ME,JB)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z, CC_TMP1, psrq, 10, 10, "Z(MB,JE)");
    dpd_buf4_close(&Z);

    /** t_I^E ( <MB||JE> + Z1_MBJE ) --> Z2_MBIJ **/
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(MB,JE)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_buf4_axpy(&C, &Z1, -1.0);
    dpd_buf4_close(&C);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(MB,JI)");
    dpd_contract424(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z1);
    dpd_buf4_sort(&Z2, CC_TMP1, pqsr, 10, 0, "Z2(MB,IJ)");
    dpd_buf4_close(&Z2);

    /** Z1_MBJI(TMP0) - Z2_MBIJ(TMP1) --> W_MBIJ **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(MB,JI)");
    dpd_buf4_init(&Z2, CC_TMP1, 0, 10, 0, 10, 0, 0, "Z2(MB,IJ)");
    dpd_buf4_axpy(&Z1, &Z2, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 2, 0, "WMBIJ");
    dpd_buf4_axpy(&Z2, &W, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&W);

    /** t_jn^bf <mn||ef> --> Z_mbje **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(me,jb)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z, CC_TMP1, psrq, 10, 10, "Z(mb,je)");
    dpd_buf4_close(&Z);

    /** t_i^e ( <mb||je> + Z1_mbje ) --> Z2_mbij **/
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(mb,je)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_buf4_axpy(&C, &Z1, -1.0);
    dpd_buf4_close(&C);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(mb,ji)");
    dpd_contract424(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z1);
    dpd_buf4_sort(&Z2, CC_TMP1, pqsr, 10, 0, "Z2(mb,ij)");
    dpd_buf4_close(&Z2);

    /** Z1_mbji(TMP0) - Z2_mbij(TMP1) --> W_mbij **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z1(mb,ji)");
    dpd_buf4_init(&Z2, CC_TMP1, 0, 10, 0, 10, 0, 0, "Z2(mb,ij)");
    dpd_buf4_axpy(&Z1, &Z2, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 2, 0, "Wmbij");
    dpd_buf4_axpy(&Z2, &W, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&W);


    /** <Mn|Ef> t_jn^bf + <MN||EF> t_jN^bF --> Z1_MEjb(TMP0) **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ME,jb)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiajb");
    dpd_contract444(&D, &T2, &Z1, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tiaJB");
    dpd_contract444(&D, &T2, &Z1, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_close(&Z1);

    /** <Mn|Fe> t_In^Fb --> Z2_MeIb (TMP3) **/
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_sort(&D, CC_TMP1, psrq, 10, 11, "D(Me,Fn)");
    dpd_buf4_close(&D);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&T2, CC_TMP2, psrq, 10, 11, "T2(Ib,Fn)");
    dpd_buf4_close(&T2);
    dpd_buf4_init(&Z2, CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(Me,Ib)");
    dpd_buf4_init(&T2, CC_TMP2, 0, 10, 11, 10, 11, 0, "T2(Ib,Fn)");
    dpd_buf4_init(&D, CC_TMP1, 0, 10, 11, 10, 11, 0, "D(Me,Fn)");
    dpd_contract444(&D, &T2, &Z2, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&Z2);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ME,jb)");
    dpd_buf4_sort(&Z1, CC_TMP1, psrq, 10, 10, "Z1(Mb,jE)");
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z2, CC_TMP3, 0, 10, 10, 10, 10, 0, "Z(Me,Ib)");
    dpd_buf4_sort(&Z2, CC_TMP0, psrq, 10, 10, "Z(Mb,Ie)");
    dpd_buf4_close(&Z2);

    /** t_I^E ( <Mj|Eb> + Z(Mb,jE)(TMP1) ) --> Z(Mb,jI)(TMP1) **/
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z1(Mb,jE)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_axpy(&D, &Z1, 1.0);
    dpd_buf4_close(&D);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(Mb,jI)");
    dpd_contract424(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z1);
    dpd_buf4_sort(&Z, CC_TMP1, pqsr, 10, 0, "Z(Mb,Ij)");
    dpd_buf4_close(&Z);

    /** t_j^e ( <Mb|Ie> - Z(Mb,Ie)(TMP0) ) --> Z(Mb,Ij)(TMP2) **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Ie)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_axpy(&C, &Z1, -1.0);
    dpd_buf4_close(&C);
    dpd_buf4_scm(&Z1, -1.0);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_buf4_init(&Z, CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    dpd_contract424(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&Z);

    /** Z(Mb,Ij) (TMP1) + Z(Mb,Ij) (TMP2) --> W_MbIj **/
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    dpd_buf4_init(&Z2, CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(Mb,Ij)");
    dpd_buf4_axpy(&Z2, &Z1, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_buf4_axpy(&Z1, &W, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&W);

    /** t_JN^BF <mN|eF> + t_Jn^Bf <mn||ef> --> Z(mB,Je) (TMP1) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(me,JB)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAJB");
    dpd_contract444(&D, &T2, &Z, 0 , 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z, CC_TMP1, psrq, 10, 10, "Z(mB,Je)");
    dpd_buf4_close(&Z);

    /** -t_Ni^Bf <mN|fE> --> Z(mB,iE) (TMP0) **/
    dpd_buf4_init(&Z, CC_TMP2, 0, 10, 10, 10, 10, 0, "Z(mE,iB)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb");
    dpd_contract444(&D, &T2, &Z, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z, CC_TMP0, psrq, 10, 10, "Z(mB,iE)");
    dpd_buf4_close(&Z);

    /** t_i^e ( <mJ|eB> + Z(mB,Je) ) --> Z1(mB,iJ) (TMP1) **/
    dpd_buf4_init(&Z, CC_TMP1, 0, 10, 10, 10, 10, 0, "Z(mB,Je)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_axpy(&D, &Z, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z1, CC_TMP2, 0, 10, 0, 10, 0, 0, "Z(mB,Ji)");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract424(&Z, &T1, &Z1, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z);
    dpd_buf4_sort(&Z1, CC_TMP1, pqsr, 10, 0, "Z1(mB,iJ)");
    dpd_buf4_close(&Z1);

    /** t_J^E ( <mB|iE> + Z(mB,iE) ) + Z1(mB,Ij) (TMP1) --> Z2(mB,iJ) (TMP2) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(mB,iE)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_axpy(&C, &Z, 1.0);
    dpd_buf4_close(&C);
    dpd_buf4_init(&Z2, CC_TMP2, 0, 10, 0, 10, 0, 0, "Z2(mB,iJ)");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract424(&Z, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z1, CC_TMP1, 0, 10, 0, 10, 0, 0, "Z1(mB,iJ)");
    dpd_buf4_axpy(&Z1, &Z2, 1.0);
    dpd_buf4_close(&Z1);

    /** Z2(mB,iJ) --> W_mBiJ **/
    dpd_buf4_init(&W, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    dpd_buf4_axpy(&Z2, &W, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z2);
  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /** t_JN^BF <MN||EF> --> Z_MBJE **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,JB)");
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z, CC_TMP0, psrq, 20, 20, "Z(MB,JE)");
    dpd_buf4_close(&Z);

    /** t_I^E ( <MB||JE> + Z1_MBJE ) --> Z2_MBIJ **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(MB,JE)");
    dpd_buf4_init(&C, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    dpd_buf4_axpy(&C, &Z1, -1.0);
    dpd_buf4_close(&C);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 20, 0, 20, 0, 0, "Z1(MB,JI)");
    dpd_contract424(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z1);
    dpd_buf4_sort(&Z2, CC_TMP0, pqsr, 20, 0, "Z2(MB,IJ)");
    dpd_buf4_close(&Z2);

    /** Z1_MBJI - Z2_MBIJ --> W_MBIJ **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 20, 0, 20, 0, 0, "Z1(MB,JI)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 20, 0, 20, 0, 0, "Z2(MB,IJ)");
    dpd_buf4_axpy(&Z1, &Z2, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&W, CC_HBAR, 0, 20, 0, 20, 2, 0, "WMBIJ");
    dpd_buf4_axpy(&Z2, &W, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&W);

    /** t_jn^bf <mn||ef> --> Z_mbje **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 30, 30, 30, 30, 0, "Z(me,jb)");
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z, CC_TMP0, psrq, 30, 30, "Z(mb,je)");
    dpd_buf4_close(&Z);

    /** t_i^e ( <mb||je> + Z1_mbje ) --> Z2_mbij **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 30, 30, 30, 30, 0, "Z(mb,je)");
    dpd_buf4_init(&C, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    dpd_buf4_axpy(&C, &Z1, -1.0);
    dpd_buf4_close(&C);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 30, 10, 30, 10, 0, "Z1(mb,ji)");
    dpd_contract424(&Z1, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z1);
    dpd_buf4_sort(&Z2, CC_TMP0, pqsr, 30, 10, "Z2(mb,ij)");
    dpd_buf4_close(&Z2);

    /** Z1_mbji - Z2_mbij --> W_mbij **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 30, 10, 30, 10, 0, "Z1(mb,ji)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 30, 10, 30, 10, 0, "Z2(mb,ij)");
    dpd_buf4_axpy(&Z1, &Z2, -1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&W, CC_HBAR, 0, 30, 10, 30, 12, 0, "Wmbij");
    dpd_buf4_axpy(&Z2, &W, 1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&W);


    /** <Mn|Ef> t_jn^bf + <MN||EF> t_jN^bF --> Z1_MEjb **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 20, 30, 20, 30, 0, "Z(ME,jb)");
    dpd_buf4_init(&D, CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z, CC_TMP0, psrq, 24, 27, "Z(Mb,jE)");
    dpd_buf4_close(&Z);


    /** <Mn|Fe> t_In^Fb --> Z2_MeIb **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ib)");
    dpd_buf4_init(&D, CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&T2);
    dpd_buf4_sort(&Z, CC_TMP0, psrq, 24, 24, "Z(Mb,Ie)");
    dpd_buf4_close(&Z);


    /** t_I^E ( <Mj|Eb> + Z(Mb,jE) ) --> Z(Mb,jI) **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 24, 27, 24, 27, 0, "Z(Mb,jE)");
    dpd_buf4_init(&D, CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    dpd_buf4_axpy(&D, &Z1, 1.0);
    dpd_buf4_close(&D);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Z, CC_TMP0, 0, 24, 23, 24, 23, 0, "Z(Mb,jI)");
    dpd_contract424(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z1);
    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 24, 22, "Z1(Mb,Ij)");
    dpd_buf4_close(&Z);

    /** t_j^e ( <Mb|Ie> - Z(Mb,Ie) ) --> Z(Mb,Ij) **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Mb,Ie)");
    dpd_buf4_init(&C, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    dpd_buf4_axpy(&C, &Z1, -1.0);
    dpd_buf4_close(&C);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_buf4_init(&Z, CC_TMP0, 0, 24, 22, 24, 22, 0, "Z2(Mb,Ij)");
    dpd_contract424(&Z1, &T1, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&Z);

    /** Z(Mb,Ij) (TMP1) + Z(Mb,Ij) (TMP2) --> W_MbIj **/
    dpd_buf4_init(&Z1, CC_TMP0, 0, 24, 22, 24, 22, 0, "Z1(Mb,Ij)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 24, 22, 24, 22, 0, "Z2(Mb,Ij)");
    dpd_buf4_axpy(&Z2, &Z1, -1.0);
    dpd_buf4_close(&Z2);
    dpd_buf4_init(&W, CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    dpd_buf4_axpy(&Z1, &W, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&W);



    /** t_JN^BF <mN|eF> + t_Jn^Bf <mn||ef> --> Z(mB,Je) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,JB)");
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
    dpd_contract444(&D, &T2, &Z, 0 , 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z, CC_TMP0, psrq, 27, 24, "Z(mB,Je)");
    dpd_buf4_close(&Z);

    /** t_Ni^Bf <mN|fE> --> Z(mB,iE) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(mE,iB)");
    dpd_buf4_init(&D, CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb");
    dpd_contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z, CC_TMP0, psrq, 27, 27, "Z(mB,iE)");
    dpd_buf4_close(&Z);

    /** t_i^e ( <mJ|eB> + Z(mB,Je) ) --> Z1(mB,iJ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 27, 24, 27, 24, 0, "Z(mB,Je)");
    dpd_buf4_init(&D, CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    dpd_buf4_axpy(&D, &Z, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 27, 22, 27, 22, 0, "Z(mB,Ji)");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract424(&Z, &T1, &Z1, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z);
    dpd_buf4_sort(&Z1, CC_TMP0, pqsr, 27, 23, "Z1(mB,iJ)");
    dpd_buf4_close(&Z1);

    /** -t_J^E ( -<mB|iE> + Z(mB,iE) ) + Z1(mB,Ij) --> Z1(mB,iJ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 27, 27, 27, 27, 0, "Z(mB,iE)");
    dpd_buf4_init(&C, CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    dpd_buf4_axpy(&C, &Z, -1.0);
    dpd_buf4_close(&C);
    dpd_buf4_init(&Z2, CC_TMP0, 0, 27, 23, 27, 23, 0, "Z2(mB,iJ)");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract424(&Z, &T1, &Z2, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&T1);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z1, CC_TMP0, 0, 27, 23, 27, 23, 0, "Z1(mB,iJ)");
    dpd_buf4_axpy(&Z2, &Z1, -1.0);
    dpd_buf4_close(&Z2);

    /** Z2(mB,iJ) --> W_mBiJ **/
    dpd_buf4_init(&W, CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    dpd_buf4_axpy(&Z1, &W, 1.0);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&W);

  }
}

}} // namespace psi::cchbar
