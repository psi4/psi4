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

/* cc3_Wamef(): Compute the Wamef matrix from CC3 theory, which is
** given in spin-orbitals as:
** 
** Wamef = <am||ef> - t_n^a <nm||ef>
**
** TDC, Feb 2004
*/

void purge_Wamef(void);

void cc3_Wamef(void)
{
  dpdbuf4 F, D, W;
  dpdfile2 t1,tia,tIA;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_sort(&F, CC3_HET1, qpsr, 11, 5, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract244(&t1, &D, &W, 0, 0, 0, -1, 1);
    dpd_file2_close(&t1);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&W, CC3_HET1, qprs, 10, 5, "CC3 WAmEf (mA,Ef)");
    dpd_buf4_close(&W);
  }

  else if (params.ref == 1) { /** ROHF **/

    /** W(AM,E>F) <--- <AM||EF> **/
    /** W(am,e>f) <--- <am||ef> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 7, 11, 5, 1, "F <ai|bc>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WAMEF (AM,E>F)");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 Wamef (am,e>f)");
    dpd_buf4_close(&F);

    /** W(Am,Ef) <--- <Am|Ef> **/
    /** W(aM,eF) <--- <aM|eF> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WaMeF (aM,eF)");
    dpd_buf4_close(&F);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /* t(N,A) <NM||EF> --> W(AM,E>F) */
    dpd_buf4_init(&W, CC3_HET1, 0, 11, 7, 11, 7, 0, "CC3 WAMEF (AM,E>F)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&W, CC3_HET1, qprs, 10, 7, "CC3 WAMEF (MA,F>E)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 7, 10, 7, 0, "CC3 WAMEF (MA,F>E)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    /* t(n,a) <nm||ef> --> W(am,e>f) */
    dpd_buf4_init(&W, CC3_HET1, 0, 11, 7, 11, 7, 0, "CC3 Wamef (am,e>f)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&tia, &D, &W, 0, 0, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&W, CC3_HET1, qprs, 10, 7, "CC3 Wamef (ma,f>e)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 10, 7, 10, 7, 0, "CC3 Wamef (ma,f>e)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    /* t(N,A) <Nm|Ef> --> W(Am,Ef) */
    dpd_buf4_init(&W, CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&W, CC3_HET1, qpsr, 10, 5, "CC3 WAmEf (mA,fE)");
    dpd_buf4_close(&W);

    /* t(n,a) <nM|eF> --> W(aM,eF) */
    dpd_buf4_init(&W, CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WaMeF (aM,eF)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&tia, &D, &W, 0, 0, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&W, CC3_HET1, qpsr, 10, 5, "CC3 WaMeF (Ma,Fe)");
    dpd_buf4_close(&W);

    dpd_file2_close(&tia);
    dpd_file2_close(&tIA);

    purge_Wamef();
  }
  
  else if (params.ref == 2) {

    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    /** W(AM,E>F) <--- <AM||EF> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WAMEF (AM,E>F)");
    dpd_buf4_close(&F);

    /** W(am,e>f) <--- <am||ef> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 Wamef (am,e>f)");
    dpd_buf4_close(&F);

    /** W(Am,Ef) <--- <Am|Ef> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_close(&F);

    /** W(aM,eF) <--- <aM|eF> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WaMeF (aM,eF)");
    dpd_buf4_close(&F);

    /** W(AM,E>F) <--- tNA * <NM||EF> **/
    dpd_buf4_init(&W, CC3_HET1, 0, 21, 7, 21, 7, 0, "CC3 WAMEF (AM,E>F)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&W, CC3_HET1, qprs, 20, 7, "CC3 WAMEF (MA,F>E)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 20, 7, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    /** W(am,e>f) <--- tna * <nm||ef> **/
    dpd_buf4_init(&W, CC3_HET1, 0, 31, 17, 31, 17, 0, "CC3 Wamef (am,e>f)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&tia, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&W, CC3_HET1, qprs, 30, 17, "CC3 Wamef (ma,f>e)");
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC3_HET1, 0, 30, 17, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    dpd_buf4_scm(&W, -1.0);
    dpd_buf4_close(&W);

    /** W(Am,Ef) <--- tNA * <Nm|Ef> **/
    dpd_buf4_init(&W, CC3_HET1, 0, 26, 28, 26, 28, 0, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&W, CC3_HET1, qpsr, 27, 29, "CC3 WAmEf (mA,fE)");
    dpd_buf4_close(&W);

    /** W(aM,eF) <--- tna * <nM|eF> **/
    dpd_buf4_init(&W, CC3_HET1, 0, 25, 29, 25, 29, 0, "CC3 WaMeF (aM,eF)");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract244(&tia, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&W, CC3_HET1, qpsr, 24, 28, "CC3 WaMeF (Ma,Fe)");
    dpd_buf4_close(&W);

    dpd_file2_close(&tia);
    dpd_file2_close(&tIA);
  }
}

void purge_Wamef(void) {
  dpdfile2 FAE, Fmi, FME, Fme;
  dpdfile4 W; 
  int *occpi, *virtpi;
  int h, a, b, e, f, i, j, m, n, omit;
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

  /* Purge Wamef matrix elements */
  dpd_file4_init(&W, CC3_HET1, 0, 11, 7,"CC3 WAMEF (AM,E>F)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      a = W.params->roworb[h][ma][0];
      asym = W.params->psym[a];
      A = a - vir_off[asym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
        e = W.params->colorb[h][ef][0];
        f = W.params->colorb[h][ef][1];
        esym = W.params->rsym[e];
        fsym = W.params->ssym[f];
        E = e - vir_off[esym];
        F = f - vir_off[fsym];
        if ((A >= (virtpi[asym] - openpi[asym])) ||
            (E >= (virtpi[esym] - openpi[esym])) ||
            (F >= (virtpi[fsym] - openpi[fsym])) )
          W.matrix[h][ma][ef] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HET1, 0, 11, 7,"CC3 Wamef (am,e>f)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      m = W.params->roworb[h][ma][1];
      msym = W.params->qsym[m];
      M = m - occ_off[msym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
        if (M >=  (occpi[msym] - openpi[msym]))
          W.matrix[h][ma][ef] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HET1, 0, 11, 5,"CC3 WAmEf (Am,Ef)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      a = W.params->roworb[h][ma][0];
      m = W.params->roworb[h][ma][1];
      asym = W.params->psym[a];
      msym = W.params->qsym[m];
      M = m - occ_off[msym];
      A = a - vir_off[asym];
      for(ef=0; ef< W.params->coltot[h]; ef++) {
        e = W.params->colorb[h][ef][0];
        esym = W.params->rsym[e];
        E = e - vir_off[esym];
        if ((A >= (virtpi[asym] - openpi[asym])) ||
            (M >=  (occpi[msym] - openpi[msym])) ||
            (E >= (virtpi[esym] - openpi[esym])) )
          W.matrix[h][ma][ef] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  dpd_file4_init(&W, CC3_HET1, 0, 11, 5,"CC3 WaMeF (aM,eF)");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&W, h);
    dpd_file4_mat_irrep_rd(&W, h);
    for(ma=0; ma < W.params->rowtot[h]; ma++) {
      for(ef=0; ef< W.params->coltot[h]; ef++) {
        f = W.params->colorb[h][ef][1];
        fsym = W.params->ssym[f];
        F = f - vir_off[fsym];
        if (F >= (virtpi[fsym] - openpi[fsym]))
          W.matrix[h][ma][ef] = 0.0;
      }
    }
    dpd_file4_mat_irrep_wrt(&W, h);
    dpd_file4_mat_irrep_close(&W, h);
  }
  dpd_file4_close(&W);

  return;
}
}} // namespace psi::ccenergy
