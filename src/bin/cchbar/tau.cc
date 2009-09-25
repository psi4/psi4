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

void tau_build(void)
{
  int h, ij, ab, i, j, a, b, I, J, A, B;
  int Isym, Jsym, Asym, Bsym;
  int nirreps;
  dpdbuf4 tauIJAB, tauijab, tauIjAb, tauiJaB, tauIjbA;
  dpdbuf4 tIJAB, tijab, tIjAb;
  dpdfile2 tIA, tia;

  nirreps = moinfo.nirreps;

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/

    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_copy(&tIJAB, CC_TAMPS, "tauIJAB");
    dpd_buf4_close(&tIJAB);

    dpd_buf4_init(&tijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
    dpd_buf4_copy(&tijab, CC_TAMPS, "tauijab");
    dpd_buf4_close(&tijab);

    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_copy(&tIjAb, CC_TAMPS, "tauIjAb");
    dpd_buf4_close(&tIjAb);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&tIA);
    dpd_file2_mat_rd(&tIA);
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_mat_init(&tia);
    dpd_file2_mat_rd(&tia);

    dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");

    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&tauIJAB, h);
      dpd_buf4_mat_irrep_rd(&tauIJAB, h);

      for(ij=0; ij < tauIJAB.params->rowtot[h]; ij++) {
	i = tauIJAB.params->roworb[h][ij][0];
	j = tauIJAB.params->roworb[h][ij][1];
	I = tIA.params->rowidx[i];
	J = tIA.params->rowidx[j];
	Isym = tIA.params->psym[i];
	Jsym = tIA.params->psym[j];
	for(ab=0; ab < tauIJAB.params->coltot[h]; ab++) {
	  a = tauIJAB.params->colorb[h][ab][0];
	  b = tauIJAB.params->colorb[h][ab][1];
	  A = tIA.params->colidx[a];
	  B = tIA.params->colidx[b];
	  Asym = tIA.params->qsym[a];
	  Bsym = tIA.params->qsym[b];

	  if((Isym==Asym) && (Jsym==Bsym))
	    tauIJAB.matrix[h][ij][ab] +=
	      (tIA.matrix[Isym][I][A] * tIA.matrix[Jsym][J][B]);
	  if((Isym==Bsym) && (Jsym==Asym))
	    tauIJAB.matrix[h][ij][ab] -=
	      (tIA.matrix[Isym][I][B] * tIA.matrix[Jsym][J][A]);

	}
      }

      dpd_buf4_mat_irrep_wrt(&tauIJAB, h);
      dpd_buf4_mat_irrep_close(&tauIJAB, h);
    }

    dpd_buf4_close(&tauIJAB);

    dpd_buf4_init(&tauijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");

    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&tauijab, h);
      dpd_buf4_mat_irrep_rd(&tauijab, h);

      for(ij=0; ij < tauijab.params->rowtot[h]; ij++) {
	i = tauijab.params->roworb[h][ij][0];
	j = tauijab.params->roworb[h][ij][1];
	I = tia.params->rowidx[i];
	J = tia.params->rowidx[j];
	Isym = tia.params->psym[i];
	Jsym = tia.params->psym[j];
	for(ab=0; ab < tauijab.params->coltot[h]; ab++) {
	  a = tauijab.params->colorb[h][ab][0];
	  b = tauijab.params->colorb[h][ab][1];
	  A = tia.params->colidx[a];
	  B = tia.params->colidx[b];
	  Asym = tia.params->qsym[a];
	  Bsym = tia.params->qsym[b];

	  if((Isym==Asym) && (Jsym==Bsym))
	    tauijab.matrix[h][ij][ab] +=
	      (tia.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);
	  if((Isym==Bsym) && (Jsym==Asym))
	    tauijab.matrix[h][ij][ab] -=
	      (tia.matrix[Isym][I][B] * tia.matrix[Jsym][J][A]);

	}
      }

      dpd_buf4_mat_irrep_wrt(&tauijab, h);
      dpd_buf4_mat_irrep_close(&tauijab, h);
    }

    dpd_buf4_close(&tauijab);

    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&tauIjAb, h);
      dpd_buf4_mat_irrep_rd(&tauIjAb, h);

      for(ij=0; ij < tauIjAb.params->rowtot[h]; ij++) {
	i = tauIjAb.params->roworb[h][ij][0];
	j = tauIjAb.params->roworb[h][ij][1];
	I = tIA.params->rowidx[i];
	J = tia.params->rowidx[j];
	Isym = tIA.params->psym[i];
	Jsym = tia.params->psym[j];
	for(ab=0; ab < tauIjAb.params->coltot[h]; ab++) {
	  a = tauIjAb.params->colorb[h][ab][0];
	  b = tauIjAb.params->colorb[h][ab][1];
	  A = tIA.params->colidx[a];
	  B = tia.params->colidx[b];
	  Asym = tIA.params->qsym[a];
	  Bsym = tia.params->qsym[b];

	  if((Isym==Asym) && (Jsym==Bsym))
	    tauIjAb.matrix[h][ij][ab] +=
	      (tIA.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);

	}
      }

      dpd_buf4_mat_irrep_wrt(&tauIjAb, h);
      dpd_buf4_mat_irrep_close(&tauIjAb, h);
    }

    /* This will generate the tauBA and tauIjbA files from tauIjAb */
    dpd_buf4_sort(&tauIjAb, CC_TAMPS, pqsr, 0, 5, "tauIjbA");
    dpd_buf4_sort(&tauIjAb, CC_TAMPS, psqr, 10, 10, "tauIjAb (Ib,jA)");
    dpd_buf4_close(&tauIjAb);
    dpd_buf4_init(&tauIjbA, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjbA");
    dpd_buf4_sort(&tauIjbA, CC_TAMPS, qprs,  0, 5, "tauiJaB");
    dpd_buf4_close(&tauIjbA);

    dpd_file2_mat_close(&tIA);
    dpd_file2_close(&tIA);
    dpd_file2_mat_close(&tia);
    dpd_file2_close(&tia);
  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_copy(&tIJAB, CC_TAMPS, "tauIJAB");
    dpd_buf4_close(&tIJAB);

    dpd_buf4_init(&tijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
    dpd_buf4_copy(&tijab, CC_TAMPS, "tauijab");
    dpd_buf4_close(&tijab);

    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_copy(&tIjAb, CC_TAMPS, "tauIjAb");
    dpd_buf4_close(&tIjAb);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_mat_init(&tIA);
    dpd_file2_mat_rd(&tIA);
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_mat_init(&tia);
    dpd_file2_mat_rd(&tia);

    dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");

    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&tauIJAB, h);
      dpd_buf4_mat_irrep_rd(&tauIJAB, h);

      for(ij=0; ij < tauIJAB.params->rowtot[h]; ij++) {
	i = tauIJAB.params->roworb[h][ij][0];
	j = tauIJAB.params->roworb[h][ij][1];
	I = tIA.params->rowidx[i];
	J = tIA.params->rowidx[j];
	Isym = tIA.params->psym[i];
	Jsym = tIA.params->psym[j];
	for(ab=0; ab < tauIJAB.params->coltot[h]; ab++) {
	  a = tauIJAB.params->colorb[h][ab][0];
	  b = tauIJAB.params->colorb[h][ab][1];
	  A = tIA.params->colidx[a];
	  B = tIA.params->colidx[b];
	  Asym = tIA.params->qsym[a];
	  Bsym = tIA.params->qsym[b];

	  if((Isym==Asym) && (Jsym==Bsym))
	    tauIJAB.matrix[h][ij][ab] +=
	      (tIA.matrix[Isym][I][A] * tIA.matrix[Jsym][J][B]);
	  if((Isym==Bsym) && (Jsym==Asym))
	    tauIJAB.matrix[h][ij][ab] -=
	      (tIA.matrix[Isym][I][B] * tIA.matrix[Jsym][J][A]);

	}
      }

      dpd_buf4_mat_irrep_wrt(&tauIJAB, h);
      dpd_buf4_mat_irrep_close(&tauIJAB, h);
    }

    dpd_buf4_close(&tauIJAB);

    dpd_buf4_init(&tauijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");

    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&tauijab, h);
      dpd_buf4_mat_irrep_rd(&tauijab, h);

      for(ij=0; ij < tauijab.params->rowtot[h]; ij++) {
	i = tauijab.params->roworb[h][ij][0];
	j = tauijab.params->roworb[h][ij][1];
	I = tia.params->rowidx[i];
	J = tia.params->rowidx[j];
	Isym = tia.params->psym[i];
	Jsym = tia.params->psym[j];
	for(ab=0; ab < tauijab.params->coltot[h]; ab++) {
	  a = tauijab.params->colorb[h][ab][0];
	  b = tauijab.params->colorb[h][ab][1];
	  A = tia.params->colidx[a];
	  B = tia.params->colidx[b];
	  Asym = tia.params->qsym[a];
	  Bsym = tia.params->qsym[b];

	  if((Isym==Asym) && (Jsym==Bsym))
	    tauijab.matrix[h][ij][ab] +=
	      (tia.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);
	  if((Isym==Bsym) && (Jsym==Asym))
	    tauijab.matrix[h][ij][ab] -=
	      (tia.matrix[Isym][I][B] * tia.matrix[Jsym][J][A]);

	}
      }

      dpd_buf4_mat_irrep_wrt(&tauijab, h);
      dpd_buf4_mat_irrep_close(&tauijab, h);
    }

    dpd_buf4_close(&tauijab);

    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");

    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&tauIjAb, h);
      dpd_buf4_mat_irrep_rd(&tauIjAb, h);

      for(ij=0; ij < tauIjAb.params->rowtot[h]; ij++) {
	i = tauIjAb.params->roworb[h][ij][0];
	j = tauIjAb.params->roworb[h][ij][1];
	I = tIA.params->rowidx[i];
	J = tia.params->rowidx[j];
	Isym = tIA.params->psym[i];
	Jsym = tia.params->psym[j];
	for(ab=0; ab < tauIjAb.params->coltot[h]; ab++) {
	  a = tauIjAb.params->colorb[h][ab][0];
	  b = tauIjAb.params->colorb[h][ab][1];
	  A = tIA.params->colidx[a];
	  B = tia.params->colidx[b];
	  Asym = tIA.params->qsym[a];
	  Bsym = tia.params->qsym[b];

	  if((Isym==Asym) && (Jsym==Bsym))
	    tauIjAb.matrix[h][ij][ab] +=
	      (tIA.matrix[Isym][I][A] * tia.matrix[Jsym][J][B]);

	}
      }

      dpd_buf4_mat_irrep_wrt(&tauIjAb, h);
      dpd_buf4_mat_irrep_close(&tauIjAb, h);
    }
    dpd_buf4_close(&tauIjAb);

    dpd_file2_mat_close(&tIA);
    dpd_file2_close(&tIA);
    dpd_file2_mat_close(&tia);
    dpd_file2_close(&tia);

    /* This will generate the tauBA and tauIjbA files from tauIjAb */
    dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_buf4_sort(&tauIjAb, CC_TAMPS, pqsr, 22, 29, "tauIjbA");
    dpd_buf4_sort(&tauIjAb, CC_TAMPS, qpsr, 23, 29, "tauiJaB");
    dpd_buf4_close(&tauIjAb);

  } /** UHF **/

}

}} // namespace psi::cchbar
