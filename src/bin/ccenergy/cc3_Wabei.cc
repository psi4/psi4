/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/* cc3_Wabei(): Compute the Wabei matrix from CC3 theory, which is
** given in spin orbitals as:
**
** Wabei = <ab||ei> - P(ab) t_m^a <mb||ei> + t_i^f <ab||ef> 
**         - P(ab) t_i^f t_m^b <am||ef> + t_m^a t_n^b <mn||ei>
**         + t_m^a t_i^f t_n^b <mn||ef>
**
** The basic strategy for this code is to generate two intermediate
** quantities, Z1(Ab,EI) and Z2(Ei,Ab), which are summed in the final
** step to give the complete W(Ei,Ab) intermediate.  This is sorted
** to W(iE,bA) storage for use in the triples equations.
**
** TDC, Feb 2004
*/

void purge_Wabei(void);

void cc3_Wabei(void)
{
  int omit = 0;
  dpdfile2 T1, t1, tIA, tia;
  dpdbuf4 Z, Z1, Z2, Z3;
  dpdbuf4 B, C, D, E, F, W;
  int Gef, Gei, Gab, Ge, Gf, Gi;
  int EE, e;
  int nrows, ncols, nlinks;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_buf4_copy(&F, CC_TMP0, "CC3 Z(Ei,Ab)");
    dpd_buf4_close(&F);

    /* Z1(Ab,Ei) <-- <Ab|Ef> * t(i,f) */
/*
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "CC3 Z(Ab,Ei)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract424(&B, &t1, &Z1, 3, 1, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&Z1);
    dpd_buf4_close(&Z2);
    dpd_file2_close(&t1);
*/
    // Added new B(+)/B(-) code from cchbar 12/29/09, -TDC
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);  
    dpd_buf4_init(&Z1, CC_TMP0, 0, 11, 8, 11, 8, 0, "Z1(ei,a>=b)");
    dpd_buf4_scm(&Z1, 0.0); /* this scm is necessary for cases with empty occpi or virtpi irreps */
    for(Gef=0; Gef < moinfo.nirreps; Gef++) {    
      Gei = Gab = Gef; /* W and B are totally symmetric */
      for(Ge=0; Ge < moinfo.nirreps; Ge++) {      
        Gf = Ge ^ Gef; Gi = Gf;  /* T1 is totally symmetric */
        B.matrix[Gef] = dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gef]);
        Z1.matrix[Gef] = dpd_block_matrix(moinfo.occpi[Gi],Z1.params->coltot[Gei])
;
        nrows = moinfo.occpi[Gi]; 
        ncols = Z1.params->coltot[Gef]; 
        nlinks = moinfo.virtpi[Gf];
        if(nrows && ncols && nlinks) {        
          for(EE=0; EE < moinfo.virtpi[Ge]; EE++) {
            e = moinfo.vir_off[Ge] + EE;
            dpd_buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
            C_DGEMM('n','n',nrows,ncols,nlinks,0.5,T1.matrix[Gi][0],nlinks,
                    B.matrix[Gef][0],ncols,0.0,Z1.matrix[Gei][0],ncols);
            dpd_buf4_mat_irrep_wrt_block(&Z1,Gei,Z1.row_offset[Gei][e],moinfo.occpi[Gi]);
          }
        }
        dpd_free_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gef]);
        dpd_free_block(Z1.matrix[Gef], moinfo.occpi[Gi], Z1.params->coltot[Gei]);
      }
    }
    dpd_buf4_close(&Z1);
    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
    dpd_buf4_close(&B);

    dpd_buf4_init(&B, CC_BINTS, 0, 5, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&T1);
    dpd_file2_mat_rd(&T1);
    dpd_buf4_init(&Z2, CC_TMP0, 0, 11, 9, 11, 9, 0, "Z2(ei,a>=b)");
    dpd_buf4_scm(&Z2, 0.0); /* this scm is necessary for cases with empty occpi or virtpi irreps */
    for(Gef=0; Gef < moinfo.nirreps; Gef++) {
      Gei = Gab = Gef; /* W and B are totally symmetric */
      for(Ge=0; Ge < moinfo.nirreps; Ge++) {
        Gf = Ge ^ Gef; Gi = Gf;  /* T1 is totally symmetric */
        B.matrix[Gef] = dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gef]);
        Z2.matrix[Gef] = dpd_block_matrix(moinfo.occpi[Gi],Z2.params->coltot[Gei]);
        nrows = moinfo.occpi[Gi]; 
        ncols = Z2.params->coltot[Gef]; 
        nlinks = moinfo.virtpi[Gf];
        if(nrows && ncols && nlinks) {
          for(EE=0; EE < moinfo.virtpi[Ge]; EE++) {
            e = moinfo.vir_off[Ge] + EE;
            dpd_buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
            C_DGEMM('n','n',nrows,ncols,nlinks,0.5,T1.matrix[Gi][0],nlinks,
                    B.matrix[Gef][0],ncols,0.0,Z2.matrix[Gei][0],ncols);
            dpd_buf4_mat_irrep_wrt_block(&Z2, Gei, Z2.row_offset[Gei][e], moinfo.occpi[Gi]);
          }
        }
        dpd_free_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gef]);
        dpd_free_block(Z2.matrix[Gef], moinfo.occpi[Gi], Z2.params->coltot[Gei]);
      }
    }
    dpd_buf4_close(&Z2);
    dpd_file2_mat_close(&T1);
    dpd_file2_close(&T1);
    dpd_buf4_close(&B);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 11, 5, 11, 8, 0, "Z1(ei,a>=b)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 11, 5, 11, 9, 0, "Z2(ei,a>=b)");
    dpd_buf4_init(&W, CC_TMP0, 0, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");
    dpd_buf4_axpy(&Z1, &W, 1);
    dpd_buf4_axpy(&Z2, &W, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&Z1);

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "CC3 Z(Ab,Ei)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");

    /* Z(Mb,Ei) <-- <Mb|Ei> + <Mb|Ef> t(i,f) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "CC3 Z(Mb,Ei)");
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC3 Z(Mb,Ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &t1, &Z, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    /* Z1(Ab,Ei) <-- - t(M,A) * Z(Mb,Ei) */
    dpd_contract244(&t1, &Z, &Z1, 0, 0, 0, -1, 0);
    dpd_buf4_close(&Z);

    /* Z(Ei,Am) <-- <Ei|Am> + <Am|Ef> t(i,f) */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_sort(&C, CC_TMP0, qpsr, 11, 11, "CC3 Z(Ei,Am)");
    dpd_buf4_close(&C);
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "CC3 Z(Am,Ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_contract424(&F, &t1, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);
    dpd_buf4_sort_axpy(&Z, CC_TMP0, rspq, 11, 11, "CC3 Z(Ei,Am)", 1);
    dpd_buf4_close(&Z);
    /* Z2(Ei,Ab) <-- - Z(Ei,Am) t(m,b) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "CC3 Z(Ei,Am)");
    dpd_contract424(&Z, &t1, &Z2, 3, 0, 0, -1, 1);
    dpd_buf4_close(&Z);

    /* Z(Mn,Ei) = <Mn|Ei> + <Mn|Ef> t_i^f */
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    dpd_buf4_copy(&E, CC_TMP0, "CC3 Z(Mn,Ei)");
    dpd_buf4_close(&E);
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "CC3 Z(Mn,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &t1, &Z, 3, 1, 0, 1, 1);
    dpd_buf4_close(&D);
    /* Z'(An,Ei) = t_M^A Z(Mn,Ei) */
    dpd_buf4_init(&Z3, CC_TMP0, 0, 11, 11, 11, 11, 0, "CC3 Z(An,Ei)");
    dpd_contract244(&t1, &Z, &Z3, 0, 0, 0, 1, 0);
    dpd_buf4_close(&Z);
    /* Z2(Ei,Ab) <-- Z'(An,Ei) t_n^b */
    dpd_contract424(&Z3, &t1, &Z2, 1, 0, 0, 1, 1);
    dpd_buf4_close(&Z3);

    dpd_buf4_close(&Z2);
    dpd_buf4_close(&Z1);

    /* W(Ab,Ei) = Z1(Ab,Ei) + Z2(Ei,Ab) */
    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "CC3 Z(Ab,Ei)");
    dpd_buf4_sort_axpy(&Z1, CC_TMP0, rspq, 11, 5, "CC3 Z(Ei,Ab)", 1);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z2, CC_TMP0, 0, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");
    dpd_buf4_sort(&Z2, CC3_HET1, qpsr, 10, 5, "CC3 WAbEi (iE,bA)");
    dpd_buf4_close(&Z2);

    dpd_file2_close(&t1);
  }

  else if (params.ref == 1) { /* ROHF */
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /** term 1 **/

    /* W(A>B,EI) <--- <AB||EI> */
    /* W(a>b,ei) <--- <ab||ei> */
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 7, 11, 5, 1, "F <ai|bc>");
    dpd_buf4_sort(&F, CC_TMP2, rspq, 7, 11, "CC3 WABEI (A>B,EI)");
    dpd_buf4_sort(&F, CC_TMP2, rspq, 7, 11, "CC3 Wabei (a>b,ei)");
    dpd_buf4_close(&F);

    /* W(Ab,Ei) <--- <Ab|Ei> */
    /* W(aB,eI) <--- <aB|eI> */
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_buf4_sort(&F, CC_TMP2, rspq, 5, 11, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_sort(&F, CC_TMP2, rspq, 5, 11, "CC3 WaBeI (aB,eI)");
    dpd_buf4_close(&F);

    /** term 2 **/

    /** W(A>B,EI) <--- <AB||EF> * t1[I][F] **/
    /** W(a>b,ei) <--- <ab||ef> * t1[i][f] **/
    dpd_buf4_init(&W, CC_TMP2, 0, 7, 11, 7, 11, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 5, 5, 5, 1, "B <ab|cd>");
    dpd_contract424(&B, &tIA, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP2, 0, 7, 11, 7, 11, 0, "CC3 Wabei (a>b,ei)");
    dpd_contract424(&B, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&B);
    dpd_buf4_close(&W);

    /** W(Ab,Ei) <--- <Ab|Ef> * t1[i][f] **/
    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract424(&B, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&B);
    dpd_buf4_close(&W);

    /** W(aB,eI) <--- t1[I][F] * <aB|eF>  **/
    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WaBeI (aB,eI)");
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract424(&B, &tIA, &W, 3, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&B);
    dpd_buf4_close(&W);

    /** term 3 **/

    /** W(A>B,EI) <--- P(A/B)( -<BM||IE> * t1[M][A] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 5, 11, 5, 11, 0, "CC3 ZABEI (AB,EI)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    dpd_contract244(&tIA, &C, &Z, 0, 0, 0, 1, 0);
    dpd_buf4_close(&C);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 7, 11, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z, CC_TMP2, qprs, 7, 11, "CC3 WABEI (A>B,EI)", -1);
    dpd_buf4_close(&Z);

    /** W(a>b,ei) <--- P(a/b)( -<bm||ie> * t1[m][a] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 5, 11, 5, 11, 0, "CC3 Zabei (ab,ei)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia||jb> (ia,bj)");
    dpd_contract244(&tia, &C, &Z, 0, 0, 0, 1, 0);
    dpd_buf4_close(&C);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 7, 11, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z, CC_TMP2, qprs, 7, 11, "CC3 Wabei (a>b,ei)", -1);
    dpd_buf4_close(&Z);

    /** W(Ab,Ei) <---  -t1[M][A] * <Mb|Ei>  - <mA|iE> * t1[m][b] **/
    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
 
    /* sort the C integrals */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 11, 10, 11, 0, "C <ia|jb> (ia,bj)");
    dpd_buf4_sort(&C, CC_CINTS, qprs, 11, 11, "C <ia|jb> (ai,bj)");
    dpd_buf4_close(&C);

    dpd_buf4_init(&C, CC_CINTS, 0, 11, 11, 11, 11, 0, "C <ia|jb> (ai,bj)");
    dpd_contract424(&C, &tia, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&C);

    dpd_buf4_close(&W);

    /** W(aB,eI) <---  - t1[m][a] * <mB|eI>  - <aM|eI> * t1[M][B] **/
    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WaBeI (aB,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj)");
    dpd_contract244(&tia, &D, &W, 0, 0, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 5, 11, 5, 0, "CC3 ZaBeI (eI,aB)");
    dpd_buf4_init(&C, CC_CINTS, 0, 11, 11, 11, 11, 0, "C <ia|jb> (ai,bj)");
    dpd_contract424(&C, &tIA, &Z, 1, 0, 0, -1, 0);
    dpd_buf4_close(&C);
    dpd_buf4_sort_axpy(&Z, CC_TMP2, rspq, 5, 11, "CC3 WaBeI (aB,eI)",1);
    dpd_buf4_close(&Z);

    /** term 4 **/

    /** W(A>B,EI) <--- -P(A/B)( <BM||FE>*t1[I][F]*t1[M][A] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(MB,EI)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z1(BA,EI)");
    dpd_contract244(&tIA, &Z, &Z1, 0, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 7, 11, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z1, &W, -1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z1, CC_TMP2, qprs, 7, 11, "CC3 WABEI (A>B,EI)", 1);
    dpd_buf4_close(&Z1);

    /** W(a>b,ei) <--- -P(a/b)( <bm||fe>*t1[i][f]*t1[m][a] ) **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(mb,ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_contract424(&F, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z1(ba,ei)");
    dpd_contract244(&tia, &Z, &Z1, 0, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 7, 11, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z1, &W, -1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z1, CC_TMP2, qprs, 7, 11, "CC3 Wabei (a>b,ei)", 1);
    dpd_buf4_close(&Z1);

    /** W(Ab,Ei) <--- -<bM|fE>*t1[i][f]*t1[M][A] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(Mb,Ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z1(Ab,Ei)");
    dpd_contract244(&tIA, &Z, &Z1, 0, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_axpy(&Z1, &W, -1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z1);

    /** W(Ab,Ei) <--- -<Am|Ef>*t1[i][f]*t1[m][b] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 5, 5, 5, 5, 0, "Z(Ef,Ab)");
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_contract424(&F, &tia, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z1(Ab,Ei)");
    dpd_contract424(&Z, &tia, &Z1, 1, 1, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_axpy(&Z1, &W, -1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z1);

    /** W(aB,eI) <--- -<Bm|Fe>*t1[I][F]*t1[m][a] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "Z(mB,eI)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z1(aB,eI)");
    dpd_contract244(&tia, &Z, &Z1, 0, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WaBeI (aB,eI)");
    dpd_buf4_axpy(&Z1, &W, -1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z1);

    /** W(aB,eI) <--- -<aM|eF>*t1[I][F]*t1[M][B] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(aM,eI)");
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 11, 5, 11, 5, 0, "Z1(eI,aB)");
    dpd_contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_sort_axpy(&Z1, CC_TMP2, rspq, 5, 11, "CC3 WaBeI (aB,eI)", -1);
    dpd_buf4_close(&Z1);

    /** term 5 **/

    /** W(A>B,EI) <---  0.5 * P(A/B)( <NM||EI>*t1[M][B]*t1[N][A] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z(EI,NB)");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    dpd_contract424(&E, &tIA, &Z, 1, 0, 0, -1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z1(AB,EI)");
    dpd_contract244(&tIA, &Z, &Z1, 0, 2, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 7, 11, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z1, &W, 0.5);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z1, CC_TMP2, qprs, 7, 11, "CC3 WABEI (A>B,EI)", -0.5);
    dpd_buf4_close(&Z1);

    /** W(a>b,ei) <--- 0.5 * P(a/b)( <nm||ei>*t1[m][b]*t1[n][a] ) **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z(ei,nb)");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 2, 11, 0, "E <ij||ka> (i>j,ak)");
    dpd_contract424(&E, &tia, &Z, 1, 0, 0, -1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z1(ab,ei)");
    dpd_contract244(&tia, &Z, &Z1, 0, 2, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 7, 11, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z1, &W, 0.5);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z1, CC_TMP2, qprs, 7, 11, "CC3 Wabei (a>b,ei)", -0.5);
    dpd_buf4_close(&Z1);

    /** W(Ab,Ei) <--- 0.5 * ( <Nm|Ei>*t1[m][b]*t1[N][A] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z(Ei,Nb)");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    dpd_contract424(&E, &tia, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_contract244(&tIA, &Z, &W, 0, 2, 0, 0.5, 1);
    dpd_buf4_close(&W);

    /** W(Ab,Ei) <--- 0.5 * ( <Mn|Ei>*t1[M][A]*t1[n][b] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z(Ei,Mb)");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    dpd_contract424(&E, &tia, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_contract244(&tIA, &Z, &W, 0, 2, 0, 0.5, 1);
    dpd_buf4_close(&W);

    /** W(aB,eI) <---  0.5 * ( <nM|eI>*t1[M][B]*t1[n][a] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z(eI,nB)");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    dpd_contract424(&E, &tIA, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WaBeI (aB,eI)");
    dpd_contract244(&tia, &Z, &W, 0, 2, 0, 0.5, 1);
    dpd_buf4_close(&W);

    /** W(aB,eI) <--- 0.5 * ( <mN|eI>*t1[m][a]*t1[N][B] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z(eI,mB)");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    dpd_contract424(&E, &tIA, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WaBeI (aB,eI)");
    dpd_contract244(&tia, &Z, &W, 0, 2, 0, 0.5, 1);
    dpd_buf4_close(&W);

    /** term 6 **/

    /** W(A>B,EI) <---  0.5 * P(A/B) <NM||EF> * t1[M][B]*t1[I][F]*t1[N][A] **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z (NM,EI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z1 (EI,NB)");
    dpd_contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z2, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z2 (AB,EI)");
    dpd_contract244(&tIA, &Z1, &Z2, 0, 2, 0, 0.5, 0);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 7, 11, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z2, &W, 1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z2, CC_TMP2, qprs, 7, 11, "CC3 WABEI (A>B,EI)", -1);
    dpd_buf4_close(&Z2);

    /** W(a>b,ei) <---  0.5 * P(a/b) <nm||ef> * t1[m][b]*t1[i][f]*t1[n][a] **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z (nm,ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z1 (ei,nb)");
    dpd_contract424(&Z, &tia, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z2, CC_TMP0, 0, 5, 11, 5, 11, 0, "Z2 (ab,ei)");
    dpd_contract244(&tia, &Z1, &Z2, 0, 2, 0, 0.5, 0);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 7, 11, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z2, &W, 1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z2, CC_TMP2, qprs, 7, 11, "CC3 Wabei (a>b,ei)", -1);
    dpd_buf4_close(&Z2);

    /** W(Ab,Ei) <---  <Nm|Ef> * t1[m][b]*t1[i][f]*t1[N][A] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z (Nm,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z1 (Ei,Nb)");
    dpd_contract424(&Z, &tia, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_contract244(&tIA, &Z1, &W, 0, 2, 0, 1, 1);
    dpd_buf4_close(&W);

    /** W(aB,eI) <---  <nM|eF> * t1[M][B]*t1[I][F]*t1[n][a] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z (nM,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z1 (eI,nB)");
    dpd_contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WaBeI (aB,eI)");
    dpd_contract244(&tia, &Z1, &W, 0, 2, 0, 1, 1);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    /* sort to Wabei (ei,ab) */
    dpd_buf4_init(&W, CC_TMP2, 0, 7, 11, 7, 11, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_sort(&W, CC_TMP2, rspq, 11, 7, "CC3 WABEI (EI,A>B)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP2, 0, 7, 11, 7, 11, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_sort(&W, CC_TMP2, rspq, 11, 7, "CC3 Wabei (ei,a>b)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_sort(&W, CC_TMP2, rspq, 11, 5, "CC3 WAbEi (Ei,Ab)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP2, 0, 5, 11, 5, 11, 0, "CC3 WaBeI (aB,eI)");
    dpd_buf4_sort(&W, CC_TMP2, rspq, 11, 5, "CC3 WaBeI (eI,aB)");
    dpd_buf4_close(&W);

    /* purge before final sort */

    purge_Wabei();

    /* sort to Wabei (ie,ba) */
    dpd_buf4_init(&W, CC_TMP2, 0, 11, 7, 11, 7, 0, "CC3 WABEI (EI,A>B)");
    dpd_buf4_sort(&W, CC3_HET1, qprs, 10, 7, "CC3 WABEI (IE,B>A)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 7, 10, 7, 0, "CC3 WABEI (IE,B>A)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP2, 0, 11, 7, 11, 7, 0, "CC3 Wabei (ei,a>b)");
    dpd_buf4_sort(&W, CC3_HET1, qprs, 10, 7, "CC3 Wabei (ie,b>a)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 7, 10, 7, 0, "CC3 Wabei (ie,b>a)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP2, 0, 11, 5, 11, 5, 0, "CC3 WAbEi (Ei,Ab)");
    dpd_buf4_sort(&W, CC3_HET1, qpsr, 10, 5, "CC3 WAbEi (iE,bA)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP2, 0, 11, 5, 11, 5, 0, "CC3 WaBeI (eI,aB)");
    dpd_buf4_sort(&W, CC3_HET1, qpsr, 10, 5, "CC3 WaBeI (Ie,Ba)");
    dpd_buf4_close(&W);
  }

  else if (params.ref == 2) { /* UHF */

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /** term 1 **/

    /* W(A>B,EI) <--- <AB||EI> */
    dpd_buf4_init(&F, CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    dpd_buf4_sort(&F, CC_TMP0, rspq, 7, 21, "CC3 WABEI (A>B,EI)");
    dpd_buf4_close(&F);

    /* W(a>b,ei) <--- <ab||ei> */
    dpd_buf4_init(&F, CC_FINTS, 0, 31, 17, 31,15, 1, "F <ai|bc>");
    dpd_buf4_sort(&F, CC_TMP0, rspq, 17, 31, "CC3 Wabei (a>b,ei)");
    dpd_buf4_close(&F);

    /* W(Ab,Ei) <--- <Ab|Ei> */
    dpd_buf4_init(&F, CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    dpd_buf4_copy(&F, CC_TMP0, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_close(&F);

    /* W(aB,eI) <--- <aB|eI> */
    dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    dpd_buf4_sort(&F, CC_TMP0, psrq, 29, 25, "CC3 WaBeI (aB,eI)");
    dpd_buf4_close(&F);

    /** term 2 **/

    /** W(A>B,EI) <--- <AB||EF> * t1[I][F] **/
    dpd_buf4_init(&W, CC_TMP0, 0, 7, 21, 7, 21, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_init(&B, CC_BINTS, 0, 7, 5, 5, 5, 1, "B <AB|CD>");
    dpd_contract424(&B, &tIA, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&B);
    dpd_buf4_close(&W);

    /** W(a>b,ei) <--- <ab||ef> * t1[i][f] **/
    dpd_buf4_init(&W, CC_TMP0, 0, 17, 31, 17, 31, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_init(&B, CC_BINTS, 0, 17, 15, 15, 15, 1, "B <ab|cd>");
    dpd_contract424(&B, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&B);
    dpd_buf4_close(&W);

    /** W(Ab,Ei) <--- <Ab|Ef> * t1[i][f] **/
    dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    dpd_contract424(&B, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&B);
    dpd_buf4_close(&W);

    /** W(aB,eI) <--- t1[I][F] * <aB|eF>  **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 24, 28, 24, 28, 0, "CC3 ZIeBa (Ie,Ba)");
    dpd_buf4_init(&B, CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
    dpd_contract244(&tIA, &B, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&B);
    dpd_buf4_sort_axpy(&Z, CC_TMP0, srqp, 29, 25, "CC3 WaBeI (aB,eI)", 1);
    dpd_buf4_close(&Z);

    /** term 3 **/

    /** W(A>B,EI) <--- P(A/B)( -<BM||IE> * t1[M][A] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 5, 21, 5, 21, 0, "CC3 ZABEI (AB,EI)");
    dpd_buf4_init(&C, CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
    dpd_contract244(&tIA, &C, &Z, 0, 0, 0, 1, 0);
    dpd_buf4_close(&C);

    dpd_buf4_init(&W, CC_TMP0, 0, 5, 21, 7, 21, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z, CC_TMP0, qprs, 7, 21, "CC3 WABEI (A>B,EI)", -1);
    dpd_buf4_close(&Z);

    /** W(a>b,ei) <--- P(a/b)( -<bm||ie> * t1[m][a] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 15, 31, 15, 31, 0, "CC3 Zabei (ab,ei)");
    dpd_buf4_init(&C, CC_CINTS, 0, 30, 31, 30, 31, 0, "C <ia||jb> (ia,bj)");
    dpd_contract244(&tia, &C, &Z, 0, 0, 0, 1, 0);
    dpd_buf4_close(&C);

    dpd_buf4_init(&W, CC_TMP0, 0, 15, 31, 17, 31, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z, CC_TMP0, qprs, 17, 31, "CC3 Wabei (a>b,ei)", -1);
    dpd_buf4_close(&Z);

    /** W(Ab,Ei) <---  -t1[M][A] * <Mb|Ei>  - <mA|iE> * t1[m][b] **/
    dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);

    dpd_buf4_init(&C, CC_CINTS, 0, 26, 26, 26, 26, 0, "C <Ai|Bj>");
    dpd_contract424(&C, &tia, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&C);

    dpd_buf4_close(&W);

    /** W(aB,eI) <---  - t1[m][a] * <mB|eI>  - <aM|eI> * t1[M][B] **/
    dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "CC3 WaBeI (aB,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    dpd_contract244(&tia, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /* sort the C integrals */
    dpd_buf4_init(&C, CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    dpd_buf4_sort(&C, CC_CINTS, qpsr, 25, 25, "C <aI|bJ>");
    dpd_buf4_close(&C);

    dpd_buf4_init(&Z, CC_TMP0, 0, 25, 29, 25, 29, 0, "CC3 ZaBeI (eI,aB)");
    dpd_buf4_init(&C, CC_CINTS, 0, 25, 25, 25, 25, 0, "C <aI|bJ>");
    dpd_contract424(&C, &tIA, &Z, 1, 0, 0, -1, 0);
    dpd_buf4_close(&C);
    dpd_buf4_sort_axpy(&Z, CC_TMP0, rspq, 29, 25, "CC3 WaBeI (aB,eI)",1);
    dpd_buf4_close(&Z);

    /** term 4 **/

    /** W(A>B,EI) <--- -P(A/B)( <BM||FE>*t1[I][F]*t1[M][A] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(BA,EI)");
    dpd_contract244(&tIA, &Z, &Z1, 0, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP0, 0, 5, 21, 7, 21, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z1, &W, -1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 7, 21, "CC3 WABEI (A>B,EI)", 1);
    dpd_buf4_close(&Z1);

    /** W(a>b,ei) <--- -P(a/b)( <bm||fe>*t1[i][f]*t1[m][a] ) **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 30, 31, 30, 31, 0, "Z(mb,ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_contract424(&F, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 15, 31, 15, 31, 0, "Z1(ba,ei)");
    dpd_contract244(&tia, &Z, &Z1, 0, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP0, 0, 15, 31, 17, 31, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z1, &W, -1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 17, 31, "CC3 Wabei (a>b,ei)", 1);
    dpd_buf4_close(&Z1);

    /** W(Ab,Ei) <--- -<bM|fE>*t1[i][f]*t1[M][A] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 24, 26, 24, 26, 0, "Z(Mb,Ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract424(&F, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 28, 26, 28, 26, 0, "Z1(Ab,Ei)");
    dpd_contract244(&tIA, &Z, &Z1, 0, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_axpy(&Z1, &W, -1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z1);

    /** W(Ab,Ei) <--- -<Am|Ef>*t1[i][f]*t1[m][b] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 26, 26, 26, 26, 0, "Z(Am,Ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    dpd_contract424(&F, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 26, 28, 26, 28, 0, "Z1(Ei,Ab)");
    dpd_contract424(&Z, &tia, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);
    dpd_buf4_sort_axpy(&Z1, CC_TMP0, rspq, 28, 26, "CC3 WAbEi (Ab,Ei)", -1);
    dpd_buf4_close(&Z1);

    /** W(aB,eI) <--- -<Bm|Fe>*t1[I][F]*t1[m][a] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 29, 25, 29, 25, 0, "Z1(aB,eI)");
    dpd_contract244(&tia, &Z, &Z1, 0, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "CC3 WaBeI (aB,eI)");
    dpd_buf4_axpy(&Z1, &W, -1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z1);

    /** W(aB,eI) <--- -<aM|eF>*t1[I][F]*t1[M][B] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
    dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    dpd_contract424(&F, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 25, 29, 25, 29, 0, "Z1(eI,aB)");
    dpd_contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, rspq, 29, 25, "CC3 WaBeI (aB,eI)", -1);
    dpd_buf4_close(&Z1);

    /** term 5 **/

    /** W(A>B,EI) <---  0.5 * P(A/B)( <NM||EI>*t1[M][B]*t1[N][A] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 21, 20, 21, 20, 0, "Z(EI,NB)");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");

    dpd_contract424(&E, &tIA, &Z, 1, 0, 0, -1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(AB,EI)");
    dpd_contract244(&tIA, &Z, &Z1, 0, 2, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP0, 0, 5, 21, 7, 21, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z1, &W, 0.5);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 7, 21, "CC3 WABEI (A>B,EI)", -0.5);
    dpd_buf4_close(&Z1);

    /** W(a>b,ei) <--- 0.5 * P(a/b)( <nm||ei>*t1[m][b]*t1[n][a] ) **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 31, 30, 31, 30, 0, "Z(ei,nb)");
    dpd_buf4_init(&E, CC_EINTS, 0, 10, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    dpd_contract424(&E, &tia, &Z, 1, 0, 0, -1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 15, 31, 15, 31, 0, "Z1(ab,ei)");
    dpd_contract244(&tia, &Z, &Z1, 0, 2, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP0, 0, 15, 31, 17, 31, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z1, &W, 0.5);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 17, 31, "CC3 Wabei (a>b,ei)", -0.5);
    dpd_buf4_close(&Z1);

    /** W(Ab,Ei) <--- 0.5 * ( <Nm|Ei>*t1[m][b]*t1[N][A] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 26, 24, 26, 24, 0, "Z(Ei,Nb)");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_contract424(&E, &tia, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_contract244(&tIA, &Z, &W, 0, 2, 0, 0.5, 1);
    dpd_buf4_close(&W);

    /** W(Ab,Ei) <--- 0.5 * ( <Mn|Ei>*t1[M][A]*t1[n][b] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 26, 24, 26, 24, 0, "Z(Ei,Mb)");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_contract424(&E, &tia, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_contract244(&tIA, &Z, &W, 0, 2, 0, 0.5, 1);
    dpd_buf4_close(&W);

    /** W(aB,eI) <---  0.5 * ( <nM|eI>*t1[M][B]*t1[n][a] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 25, 27, 25, 27, 0, "Z(eI,nB)");
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    dpd_contract424(&E, &tIA, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "CC3 WaBeI (aB,eI)");
    dpd_contract244(&tia, &Z, &W, 0, 2, 0, 0.5, 1);
    dpd_buf4_close(&W);

    /** W(aB,eI) <--- 0.5 * ( <mN|eI>*t1[m][a]*t1[N][B] ) **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 25, 27, 25, 27, 0, "Z(eI,mB)");
    dpd_buf4_init(&E, CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
    dpd_contract424(&E, &tIA, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "CC3 WaBeI (aB,eI)");
    dpd_contract244(&tia, &Z, &W, 0, 2, 0, 0.5, 1);
    dpd_buf4_close(&W);

    /** term 6 **/

    /** W(A>B,EI) <---  0.5 * P(A/B) <NM||EF> * t1[M][B]*t1[I][F]*t1[N][A] **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 21, 0, 21, 0, "Z (NM,EI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 21, 20, 21, 20, 0, "Z1 (EI,NB)");
    dpd_contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z2, CC_TMP0, 0, 5, 21, 5, 21, 0, "Z2 (AB,EI)");
    dpd_contract244(&tIA, &Z1, &Z2, 0, 2, 0, 0.5, 0);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC_TMP0, 0, 5, 21, 7, 21, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_axpy(&Z2, &W, 1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z2, CC_TMP0, qprs, 7, 21, "CC3 WABEI (A>B,EI)", -1);
    dpd_buf4_close(&Z2);

    /** W(a>b,ei) <---  0.5 * P(a/b) <nm||ef> * t1[m][b]*t1[i][f]*t1[n][a] **/
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 31, 10, 31, 0, "Z (nm,ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 31, 30, 31, 30, 0, "Z1 (ei,nb)");
    dpd_contract424(&Z, &tia, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z2, CC_TMP0, 0, 15, 31, 15, 31, 0, "Z2 (ab,ei)");
    dpd_contract244(&tia, &Z1, &Z2, 0, 2, 0, 0.5, 0);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC_TMP0, 0, 15, 31, 17, 31, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_axpy(&Z2, &W, 1);
    dpd_buf4_close(&W);

    dpd_buf4_sort_axpy(&Z2, CC_TMP0, qprs, 17, 31, "CC3 Wabei (a>b,ei)", -1);
    dpd_buf4_close(&Z2);

    /** W(Ab,Ei) <---  <Nm|Ef> * t1[m][b]*t1[i][f]*t1[N][A] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 22, 26, 22, 26, 0, "Z (Nm,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 26, 24, 26, 24, 0, "Z1 (Ei,Nb)");
    dpd_contract424(&Z, &tia, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_contract244(&tIA, &Z1, &W, 0, 2, 0, 1, 1);
    dpd_buf4_close(&W);

    /** W(aB,eI) <---  <nM|eF> * t1[M][B]*t1[I][F]*t1[n][a] **/

    dpd_buf4_init(&Z, CC_TMP0, 0, 23, 25, 23, 25, 0, "Z (nM,eI)");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z1, CC_TMP0, 0, 25, 27, 25, 27, 0, "Z1 (eI,nB)");
    dpd_contract424(&Z, &tIA, &Z1, 1, 0, 0, 1, 0);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "CC3 WaBeI (aB,eI)");
    dpd_contract244(&tia, &Z1, &W, 0, 2, 0, 1, 1);
    dpd_buf4_close(&W);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    /* sort to Wabei (ei,ab) */
    dpd_buf4_init(&W, CC_TMP0, 0, 7, 21, 7, 21, 0, "CC3 WABEI (A>B,EI)");
    dpd_buf4_sort(&W, CC_TMP2, rspq, 21, 7, "CC3 WABEI (EI,A>B)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, 0, 17, 31, 17, 31, 0, "CC3 Wabei (a>b,ei)");
    dpd_buf4_sort(&W, CC_TMP2, rspq, 31, 17, "CC3 Wabei (ei,a>b)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, 0, 28, 26, 28, 26, 0, "CC3 WAbEi (Ab,Ei)");
    dpd_buf4_sort(&W, CC_TMP2, rspq, 26, 28, "CC3 WAbEi (Ei,Ab)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_TMP0, 0, 29, 25, 29, 25, 0, "CC3 WaBeI (aB,eI)");
    dpd_buf4_sort(&W, CC_TMP2, rspq, 25, 29, "CC3 WaBeI (eI,aB)");
    dpd_buf4_close(&W);

    /* sort to Wabei (ie,ab) */
    dpd_buf4_init(&W, CC_TMP2, 0, 21, 7, 21, 7, 0, "CC3 WABEI (EI,A>B)");
    dpd_buf4_sort(&W, CC3_HET1, qprs, 20, 7, "CC3 WABEI (IE,B>A)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 20, 7, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP2, 0, 31, 17, 31, 17, 0, "CC3 Wabei (ei,a>b)");
    dpd_buf4_sort(&W, CC3_HET1, qprs, 30, 17, "CC3 Wabei (ie,b>a)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 17, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP2, 0, 26, 28, 26, 28, 0, "CC3 WAbEi (Ei,Ab)");
    dpd_buf4_sort(&W, CC3_HET1, qpsr, 27, 29, "CC3 WAbEi (iE,bA)");
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_TMP2, 0, 25, 29, 25, 29, 0, "CC3 WaBeI (eI,aB)");
    dpd_buf4_sort(&W, CC3_HET1, qpsr, 24, 28, "CC3 WaBeI (Ie,Ba)");
    dpd_buf4_close(&W);
  }
}

void purge_Wabei(void) {
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

  /* Purge Wabei matrix elements */
  dpd_file4_init(&W, CC_TMP2, 0, 11, 7,"CC3 WABEI (EI,A>B)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      e = W.params->roworb[h][ei][0];
      esym = W.params->psym[e];
      E = e - vir_off[esym]; 
      for(ab=0; ab<W.params->coltot[h]; ab++) {
        a = W.params->colorb[h][ab][0];
        b = W.params->colorb[h][ab][1];
        asym = W.params->rsym[a];
        bsym = W.params->ssym[b]; 
        A = a - vir_off[asym];
        B = b - vir_off[bsym];
        if ((E >= (virtpi[esym] - openpi[esym])) ||
            (A >= (virtpi[asym] - openpi[asym])) ||
            (B >= (virtpi[bsym] - openpi[bsym])) )
          W.matrix[h][ei][ab] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_TMP2, 0, 11, 7,"CC3 Wabei (ei,a>b)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      i = W.params->roworb[h][ei][1];
      isym = W.params->qsym[i]; 
      I = i - occ_off[isym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
        if (I >= (occpi[isym] - openpi[isym]))
          W.matrix[h][ei][ab] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_TMP2, 0, 11, 5,"CC3 WAbEi (Ei,Ab)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      e = W.params->roworb[h][ei][0];
      i = W.params->roworb[h][ei][1];
      esym = W.params->psym[e];
      isym = W.params->qsym[i];
      E = e - vir_off[esym];
      I = i - occ_off[isym];
      for(ab=0; ab<W.params->coltot[h]; ab++) {
        a = W.params->colorb[h][ab][0];
        asym = W.params->rsym[a];
        bsym = W.params->ssym[b];
        A = a - vir_off[asym];
        if ((E >= (virtpi[esym] - openpi[esym])) ||
            (I >= (occpi[isym] - openpi[isym])) ||
            (A >= (virtpi[asym] - openpi[asym])) )
          W.matrix[h][ei][ab] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC_TMP2, 0, 11, 5,"CC3 WaBeI (eI,aB)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ei=0; ei<W.params->rowtot[h]; ei++) {
      for(ab=0; ab<W.params->coltot[h]; ab++) {
        b = W.params->colorb[h][ab][1];
        bsym = W.params->ssym[b];
        B = b - vir_off[bsym];
        if (B >= (virtpi[bsym] - openpi[bsym]))
          W.matrix[h][ei][ab] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

}

}} // namespace psi::ccenergy
