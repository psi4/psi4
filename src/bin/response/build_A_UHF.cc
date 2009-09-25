/*! \file
    \ingroup RESPONSE
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace response {

/* BUILD_A_UHF(): Construct the components of the molecular orbital
** Hessian, A, for UHF orbitals.  Building
** all symmetry blocks of A.
**
** A(AI,JB) = delta_IJ f_AB - delta_AB f_IJ + <IJ||AB> - <JA||IB> 
**
** A(ai,jb) = delta_ij f_ab - delta_ab f_ij + <ij||ab> - <ja||ib> 
**
** A(AI,jb) = 2<Ij|Ab>
**
** TDC, January 2003
*/

void build_A_UHF(void)
{
  int h, nirreps, row, col;
  int a, i, b, j;
  int A, I, B, J;
  int Asym, Isym, Bsym, Jsym;
  dpdfile2 fIJ, fij, fAB, fab, fIA, fia;
  dpdbuf4 Amat, D, C;

  nirreps = moinfo.nirreps;

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_rd(&fij);
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fAB);
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
  dpd_file2_mat_init(&fab);
  dpd_file2_mat_rd(&fab);
  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_mat_init(&fIA);
  dpd_file2_mat_rd(&fIA);
  dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
  dpd_file2_mat_init(&fia);
  dpd_file2_mat_rd(&fia);

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 1, "D <IJ|AB>");
  dpd_buf4_sort(&D, PSIF_MO_HESS, rpsq, 21, 21, "A(AI,BJ)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
  dpd_buf4_sort_axpy(&D, PSIF_MO_HESS, qrsp, 21, 21, "A(AI,BJ)", -1.0);
  dpd_buf4_close(&D);

  dpd_buf4_init(&Amat, PSIF_MO_HESS, 0, 21, 21, 21, 21, 0, "A(AI,BJ)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Amat, h);
    dpd_buf4_mat_irrep_rd(&Amat, h);

    for(row=0; row < Amat.params->rowtot[h]; row++) {
      a = Amat.params->roworb[h][row][0];
      i = Amat.params->roworb[h][row][1];
      A = fAB.params->rowidx[a];
      I = fIJ.params->rowidx[i];
      Asym = fAB.params->psym[a];
      Isym = fIJ.params->psym[i];
      for(col=0; col < Amat.params->coltot[h]; col++) {
	b = Amat.params->colorb[h][col][0];
	j = Amat.params->colorb[h][col][1];
	B = fAB.params->colidx[b];
	J = fIJ.params->colidx[j];
	Bsym = fAB.params->qsym[b];
	Jsym = fIJ.params->qsym[j];

	if((I==J) && (Asym==Bsym))
	  Amat.matrix[h][row][col] += fAB.matrix[Asym][A][B];
	if((A==B) && (Isym==Jsym))
	  Amat.matrix[h][row][col] -= fIJ.matrix[Isym][I][J];

      }
    }

    dpd_buf4_mat_irrep_wrt(&Amat, h);
    dpd_buf4_mat_irrep_close(&Amat, h);
  }
  dpd_buf4_close(&Amat);

  dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 10, 15, 1, "D <ij|ab>");
  dpd_buf4_sort(&D, PSIF_MO_HESS, rpsq, 31, 31, "A(ai,bj)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
  dpd_buf4_sort_axpy(&D, PSIF_MO_HESS, qrsp, 31, 31, "A(ai,bj)", -1.0);
  dpd_buf4_close(&D);

  dpd_buf4_init(&Amat, PSIF_MO_HESS, 0, 31, 31, 31, 31, 0, "A(ai,bj)");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Amat, h);
    dpd_buf4_mat_irrep_rd(&Amat, h);

    for(row=0; row < Amat.params->rowtot[h]; row++) {
      a = Amat.params->roworb[h][row][0];
      i = Amat.params->roworb[h][row][1];
      A = fab.params->rowidx[a];
      I = fij.params->rowidx[i];
      Asym = fab.params->psym[a];
      Isym = fij.params->psym[i];
      for(col=0; col < Amat.params->coltot[h]; col++) {
	b = Amat.params->colorb[h][col][0];
	j = Amat.params->colorb[h][col][1];
	B = fab.params->colidx[b];
	J = fij.params->colidx[j];
	Bsym = fab.params->qsym[b];
	Jsym = fij.params->qsym[j];

	if((I==J) && (Asym==Bsym))
	  Amat.matrix[h][row][col] += fab.matrix[Asym][A][B];
	if((A==B) && (Isym==Jsym))
	  Amat.matrix[h][row][col] -= fij.matrix[Isym][I][J];

      }
    }

    dpd_buf4_mat_irrep_wrt(&Amat, h);
    dpd_buf4_mat_irrep_close(&Amat, h);
  }
  dpd_buf4_close(&Amat);

  dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  dpd_buf4_sort(&D, PSIF_MO_HESS, rpsq, 21, 31, "A(AI,bj)");
  dpd_buf4_close(&D);
  dpd_buf4_init(&Amat, PSIF_MO_HESS, 0, 21, 31, 21, 31, 0, "A(AI,bj)");
  dpd_buf4_scm(&Amat, 2);
  dpd_buf4_close(&Amat);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_close(&fIJ);
  dpd_file2_mat_close(&fij);
  dpd_file2_close(&fij);
  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fAB);
  dpd_file2_mat_close(&fab);
  dpd_file2_close(&fab);
  dpd_file2_mat_close(&fIA);
  dpd_file2_close(&fIA);
  dpd_file2_mat_close(&fia);
  dpd_file2_close(&fia);
}


}} // namespace psi::response
