/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/* cc3_Wmbij(): Compute the Wmbij components of the
** T1-similarity-transformed Hamiltonian matrix, which is given in
** spin-orbitals as:
**
** Wmbij = <mb||ij> + P(ij) t_i^e <mb||ej> - t_n^b Wmnij + t_i^e t_j^f <mb||ef>
**
** where the Wmnij intermediate is described in cc3_Wmnij.c.
**
** TDC, Feb 2004
*/

void purge_Wmbij(void);

void cc3_Wmbij(void)
{
  dpdbuf4 C, D, E, F, W, W1, Z, X, Z1;
  dpdfile2 t1, tia, tIA;

  if(params.ref == 0) { /** RHF **/

    /* W(Mb,Ij) <-- <Mb|Ij> */
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0,"E <ij|ka>");
    dpd_buf4_sort(&E, CC3_HET1, rspq, 10, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_close(&E);

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

    /* W(Mb,Ij) <-- - t(n,b) W(Mn,Ij) */
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&W1, CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract424(&W1, &t1, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&W1);
    dpd_buf4_close(&W);


    /* W(Mb,Ij) <-- + t(j,e) <Mb|Ie) */
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_contract424(&C, &t1, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&C);
    dpd_buf4_close(&W);

    /* Z(Mb,Ej) = <Mb|Ej> + t(j,f) * <Mb|Ef> */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "CC3 ZMbEj (Mb,Ej)");
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC3 ZMbEj (Mb,Ej)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &t1, &Z, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);

    /* W(Mb,Ij) <-- t(I,E) * Z(Mb,Ej) */
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC3 ZMbEj (Mb,Ej)");
    dpd_contract244(&t1, &Z, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_sort(&W, CC3_HET1, rspq, 0, 10, "CC3 WMbIj (Ij,Mb)");
    dpd_buf4_close(&W);

    dpd_file2_close(&t1);

  }

  else if (params.ref == 1) {
    /** W(MB,I>J) <--- <MB||IJ> **/
     /** W(mb,i>j) <--- <mb||ij> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&E, CC3_HET1, rspq, 10, 2, "CC3 WMBIJ (MB,I>J)");
    dpd_buf4_sort(&E, CC3_HET1, rspq, 10, 2, "CC3 Wmbij (mb,i>j)");
    dpd_buf4_close(&E);

     /** W(Mb,Ij) <--- <Mb|Ij> **/
    /** W(mB,iJ) <--- <mB|iJ> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_sort(&E, CC3_HET1, rspq, 10, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_sort(&E, CC3_HET1, rspq, 10, 0, "CC3 WmBiJ (mB,iJ)");
    dpd_buf4_close(&E);

    /**** Term II ****/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /**** W(MB,I>J) <-- -ZMBJI <-- P(I/J)( -<JE||MB> * t1[I][E] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z (MB,JI)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_contract424(&C, &tIA, &Z, 1, 1, 0, -1, 0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 10, 0, "X (MB,IJ)");
    dpd_buf4_init(&X, CC_TMP0, 0, 10, 0, 10, 0, 0, "X (MB,IJ)");
    dpd_buf4_axpy(&Z, &X, -1);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 2, 0, "CC3 WMBIJ (MB,I>J)");
    dpd_buf4_axpy(&X, &W, 1);
    dpd_buf4_close(&X);
    dpd_buf4_close(&W);

    /**** W(mb,i>j) <-- -Zmbji <-- P(i/j)( -<je||mb> * t1[i][e] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z (mb,ji)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_contract424(&C, &tia, &Z, 1, 1, 0, -1, 0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 10, 0, "X (mb,ij)");
    dpd_buf4_init(&X, CC_TMP0, 0, 10, 0, 10, 0, 0, "X (mb,ij)");
    dpd_buf4_axpy(&Z, &X, -1);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 2, 0, "CC3 Wmbij (mb,i>j)");
    dpd_buf4_axpy(&X, &W, 1);
    dpd_buf4_close(&X);
    dpd_buf4_close(&W);

    /**** W(Mb,Ij) <-- ( <Mb|Ej> * t1[I][E] ) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tIA, &D, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /**** W(Mb,Ij) <-- ( <Mb|Ie> * t1[j][e] ) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_contract424(&C, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&C);
    dpd_buf4_close(&W);

    /**** W(mB,iJ) <-- ( <mB|eJ> * t1[i][e] ) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WmBiJ (mB,iJ)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tia, &D, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /**** W(mB,iJ) <-- ( <mB|iE> * t1[J][E] ) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WmBiJ (mB,iJ)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_contract424(&C, &tIA, &W, 3, 1, 0, 1.0, 1);
    dpd_buf4_close(&C);
    dpd_buf4_close(&W);

    /**** Term III ****/

    /**** W(MB,I>J) <-- ( t1[N][B] * W(MN,I>J) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 2, 10, 2, 0, "CC3 WMBIJ (MB,I>J)");
    dpd_buf4_init(&Z, CC3_HET1, 0, 0, 2, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    dpd_contract424(&Z, &tIA, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(mb,i>j) <-- ( t1[n][b] * W(mn,i>j) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 2, 10, 2, 0, "CC3 Wmbij (mb,i>j)");
    dpd_buf4_init(&Z, CC3_HET1, 0, 0, 2, 2, 2, 0, "CC3 Wmnij (m>n,i>j)");
    dpd_contract424(&Z, &tia, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(Mb,Ij) <-- ( t1[n][b] * W(Mn,Ij) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&Z, CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract424(&Z, &tia, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(mB,iJ) <-- ( t1[N][B] * W(mN,iJ) ****/

    dpd_buf4_init(&Z1, CC_TMP0, 0, 11, 0, 11, 0, 0, "Z (Bm,Ji)");
    dpd_buf4_init(&Z, CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract244(&tIA, &Z, &Z1, 0, 0, 0, -1, 0);
    dpd_buf4_close(&Z);
    dpd_buf4_sort_axpy(&Z1, CC3_HET1, qpsr, 10, 0, "CC3 WmBiJ (mB,iJ)", 1.0);
    dpd_buf4_close(&Z1);

    /**** Term IV ****/

    /**** W(MB,I>J) <-- 0.5*P(I/J)XMBIJ <--- ( t1[I][E] * ZMBEJ ) <--  <MB||EF> * t1[J][F] ****/

    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z (MB,EJ)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 2, 0, "CC3 WMBIJ (MB,I>J)");
    dpd_contract244(&tIA, &Z, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(mb,i>j) <-- P(i/j) (Zmbif * t1[j][f]) <-- 0.5*( t1[i][e] * <mb||ef> ) ****/

    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z (mb,ej)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_contract424(&F, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 2, 0, "CC3 Wmbij (mb,i>j)");
    dpd_contract244(&tia, &Z, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(Mb,Ij) <-- (ZIfMb * t1[j][f]) <--  t1[I][E] * <Mb|Ef> ****/

    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (If,Mb)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract244(&tIA, &F, &Z, 1, 2, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_contract424(&Z, &tia, &W, 1, 1, 0, 1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(mB,iJ) <-- ZiFmB * t1[J][F] <-- t1[i][e] * <mB|eF> ****/

    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (iF,mB)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract244(&tia, &F, &Z, 1, 2, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WmBiJ (mB,iJ)");
    dpd_contract424(&Z, &tIA, &W, 1, 1, 0, 1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /* do purge before sort */
    purge_Wmbij();

    /* do final sort to get (Ij,Mb) */
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 2, 10, 2, 0, "CC3 WMBIJ (MB,I>J)");
    dpd_buf4_sort(&W, CC3_HET1, rspq, 2, 10, "CC3 WMBIJ (I>J,MB)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 2, 10, 2, 0, "CC3 Wmbij (mb,i>j)");
    dpd_buf4_sort(&W, CC3_HET1, rspq, 2, 10, "CC3 Wmbij (i>j,mb)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_sort(&W, CC3_HET1, rspq, 0, 10, "CC3 WMbIj (Ij,Mb)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 0, 10, 0, 0, "CC3 WmBiJ (mB,iJ)");
    dpd_buf4_sort(&W, CC3_HET1, rspq, 0, 10, "CC3 WmBiJ (iJ,mB)");
    dpd_buf4_close(&W); 

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  }

  else if (params.ref == 2) {

    /** W(MB,I>J) <--- <MB||IJ> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_buf4_sort(&E, CC3_HET1, rspq, 20, 2, "CC3 WMBIJ (MB,I>J)");
    dpd_buf4_close(&E);

     /** W(mb,i>j) <--- <mb||ij> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_buf4_sort(&E, CC3_HET1, rspq, 30, 12, "CC3 Wmbij (mb,i>j)");
    dpd_buf4_close(&E);

     /** W(Mb,Ij) <--- <Mb|Ij> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_buf4_sort(&E, CC3_HET1, rspq, 24, 22, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_close(&E);

    /** W(mB,iJ) <--- <mB|iJ> **/
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    dpd_buf4_sort(&E, CC3_HET1, rspq, 27, 23, "CC3 WmBiJ (mB,iJ)");
    dpd_buf4_close(&E);

    /**** Term II ****/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /**** W(MB,I>J) <-- -ZMBJI <-- P(I/J)( -<JE||MB> * t1[I][E] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 20, 0, 20, 0, 0, "Z (MB,JI)");
    dpd_buf4_init(&C, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    dpd_contract424(&C, &tIA, &Z, 1, 1, 0, -1, 0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 20, 0, "X (MB,IJ)");
    dpd_buf4_init(&X, CC_TMP0, 0, 20, 0, 20, 0, 0, "X (MB,IJ)");
    dpd_buf4_axpy(&Z, &X, -1);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&W, CC3_HET1, 0, 20, 0, 20, 2, 0, "CC3 WMBIJ (MB,I>J)");
    dpd_buf4_axpy(&X, &W, 1);
    dpd_buf4_close(&X);
    dpd_buf4_close(&W);

    /**** W(mb,i>j) <-- -Zmbji <-- P(i/j)( -<je||mb> * t1[i][e] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 30, 10, 30, 10, 0, "Z (mb,ji)");
    dpd_buf4_init(&C, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    dpd_contract424(&C, &tia, &Z, 1, 1, 0, -1, 0);
    dpd_buf4_close(&C);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 30, 10, "X (mb,ij)");
    dpd_buf4_init(&X, CC_TMP0, 0, 30, 10, 30, 10, 0, "X (mb,ij)");
    dpd_buf4_axpy(&Z, &X, -1);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 10, 30, 12, 0, "CC3 Wmbij (mb,i>j)");
    dpd_buf4_axpy(&X, &W, 1);
    dpd_buf4_close(&X);
    dpd_buf4_close(&W);

    /**** W(Mb,Ij) <-- ( <Mb|Ej> * t1[I][E] ) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 24, 22, 24, 22, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_contract244(&tIA, &D, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /**** W(Mb,Ij) <-- ( <Mb|Ie> * t1[j][e] ) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 24, 22, 24, 22, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&C, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    dpd_contract424(&C, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&C);
    dpd_buf4_close(&W);

    /**** W(mB,iJ) <-- ( <mB|eJ> * t1[i][e] ) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 27, 23, 27, 23, 0, "CC3 WmBiJ (mB,iJ)");
    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_contract244(&tia, &D, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /**** W(mB,iJ) <-- ( <mB|iE> * t1[J][E] ) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 27, 23, 27, 23, 0, "CC3 WmBiJ (mB,iJ)");
    dpd_buf4_init(&C, CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    dpd_contract424(&C, &tIA, &W, 3, 1, 0, 1.0, 1);
    dpd_buf4_close(&C);
    dpd_buf4_close(&W);

    /**** Term III ****/

    /**** W(MB,I>J) <-- ( t1[N][B] * W(MN,I>J) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 20, 2, 20, 2, 0, "CC3 WMBIJ (MB,I>J)");
    dpd_buf4_init(&Z, CC3_HET1, 0, 0, 2, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    dpd_contract424(&Z, &tIA, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(mb,i>j) <-- ( t1[n][b] * W(mn,i>j) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 30, 12, 30, 12, 0, "CC3 Wmbij (mb,i>j)");
    dpd_buf4_init(&Z, CC3_HET1, 0, 10, 12, 12, 12, 0, "CC3 Wmnij (m>n,i>j)");
    dpd_contract424(&Z, &tia, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(Mb,Ij) <-- ( t1[n][b] * W(Mn,Ij) ****/

    dpd_buf4_init(&W, CC3_HET1, 0, 24, 22, 24, 22, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&Z, CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract424(&Z, &tia, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(mB,iJ) <-- ( t1[N][B] * W(mN,iJ) ****/

    dpd_buf4_init(&Z1, CC_TMP0, 0, 26, 22, 26, 22, 0, "Z (Bm,Ji)");
    dpd_buf4_init(&Z, CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract244(&tIA, &Z, &Z1, 0, 0, 0, -1, 0);
    dpd_buf4_close(&Z);
    dpd_buf4_sort_axpy(&Z1, CC3_HET1, qpsr, 27, 23, "CC3 WmBiJ (mB,iJ)", 1.0);
    dpd_buf4_close(&Z1);

    /**** Term IV ****/

    /**** W(MB,I>J) <-- 0.5*P(I/J)XMBIJ <--- ( t1[I][E] * ZMBEJ ) <--  <MB||EF> * t1[J][F] ****/

    dpd_buf4_init(&Z, CC_TMP0, 0, 20, 21, 20, 21, 0, "Z (MB,EJ)");
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC3_HET1, 0, 20, 0, 20, 2, 0, "CC3 WMBIJ (MB,I>J)");
    dpd_contract244(&tIA, &Z, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(mb,i>j) <-- P(i/j) (Zmbif * t1[j][f]) <-- 0.5*( t1[i][e] * <mb||ef> ) ****/

    dpd_buf4_init(&Z, CC_TMP0, 0, 30, 31, 30, 31, 0, "Z (mb,ej)");
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_contract424(&F, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC3_HET1, 0, 30, 10, 30, 12, 0, "CC3 Wmbij (mb,i>j)");
    dpd_contract244(&tia, &Z, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(Mb,Ij) <-- (ZIfMb * t1[j][f]) <--  t1[I][E] * <Mb|Ef> ****/

    dpd_buf4_init(&Z, CC_TMP0, 0, 24, 24, 24, 24, 0, "Z (If,Mb)");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract244(&tIA, &F, &Z, 1, 2, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC3_HET1, 0, 24, 22, 24, 22, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_contract424(&Z, &tia, &W, 1, 1, 0, 1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(mB,iJ) <-- ZiFmB * t1[J][F] <-- t1[i][e] * <mB|eF> ****/

    dpd_buf4_init(&Z, CC_TMP0, 0, 27, 27, 27, 27, 0, "Z (iF,mB)");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract244(&tia, &F, &Z, 1, 2, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC3_HET1, 0, 27, 23, 27, 23, 0, "CC3 WmBiJ (mB,iJ)");
    dpd_contract424(&Z, &tIA, &W, 1, 1, 0, 1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    dpd_buf4_init(&W, CC3_HET1, 0, 20, 2, 20, 2, 0, "CC3 WMBIJ (MB,I>J)");
    dpd_buf4_sort(&W, CC3_HET1, rspq, 2, 20, "CC3 WMBIJ (I>J,MB)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 12, 30, 12, 0, "CC3 Wmbij (mb,i>j)");
    dpd_buf4_sort(&W, CC3_HET1, rspq, 12, 30, "CC3 Wmbij (i>j,mb)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 24, 22, 24, 22, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_sort(&W, CC3_HET1, rspq, 22, 24, "CC3 WMbIj (Ij,Mb)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 27, 23, 27, 23, 0, "CC3 WmBiJ (mB,iJ)");
    dpd_buf4_sort(&W, CC3_HET1, rspq, 23, 27, "CC3 WmBiJ (iJ,mB)");
    dpd_buf4_close(&W);
  }
}

void purge_Wmbij(void) {
  dpdfile4 W;
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n;
  int    A, B, E, F, I, J, M, N;
  int mn, ei, ma, ef, me, jb, mb, ij, ab;
  int asym, bsym, esym, fsym, isym, jsym, msym, nsym;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi, nirreps;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  dpd_file4_init(&W, CC3_HET1, 0, 10, 2,"CC3 WMBIJ (MB,I>J)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      b = W.params->roworb[h][mb][1];
      bsym = W.params->qsym[b];
      B = b - vir_off[bsym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
        if (B >= (virtpi[bsym] - openpi[bsym]))
          W.matrix[h][mb][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HET1, 0, 10, 2,"CC3 Wmbij (mb,i>j)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      m = W.params->roworb[h][mb][0];
      msym = W.params->psym[m];
      M = m - occ_off[msym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
        i = W.params->colorb[h][ij][0];
        j = W.params->colorb[h][ij][1];
        isym = W.params->rsym[i];
        jsym = W.params->ssym[j];
        I = i - occ_off[isym];
        J = j - occ_off[jsym];
        if ((M >= (occpi[msym] - openpi[msym])) ||
            (I >= (occpi[isym] - openpi[isym])) ||
            (J >= (occpi[jsym] - openpi[jsym])) )
          W.matrix[h][mb][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HET1, 0, 10, 0,"CC3 WMbIj (Mb,Ij)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      for(ij=0; ij<W.params->coltot[h]; ij++) {
        j = W.params->colorb[h][ij][1];
        jsym = W.params->ssym[j];
        J = j - occ_off[jsym];
        if (J >= (occpi[jsym] - openpi[jsym]))
          W.matrix[h][mb][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HET1, 0, 10, 0,"CC3 WmBiJ (mB,iJ)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(mb=0; mb<W.params->rowtot[h]; mb++) {
      m = W.params->roworb[h][mb][0];
      b = W.params->roworb[h][mb][1];
      msym = W.params->psym[m];
      bsym = W.params->qsym[b];
      M = m - occ_off[msym];
      B = b - vir_off[bsym];
      for(ij=0; ij<W.params->coltot[h]; ij++) {
        i = W.params->colorb[h][ij][0];
        isym = W.params->rsym[i];
        I = i - occ_off[isym];
        if ((M >= (occpi[msym] - openpi[msym])) ||
            (B >= (virtpi[bsym] - openpi[bsym])) ||
            (I >= (occpi[isym] - openpi[isym])) )
          W.matrix[h][mb][ij] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);
}

}} // namespace psi::ccenergy
