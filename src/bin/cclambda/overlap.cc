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

void overlap(int L_irr)
{
  int h, nirreps;
  int row, col;
  int i,j,a,b,I,J,A,B,Isym,Jsym,Asym,Bsym;
  dpdfile2 T1, L1, T1A, T1B;
  dpdbuf4 T2, L2;
  double value = 1.0;
  double ST1A, ST1B, ST2AA, ST2BB, ST2AB, ST12AA, ST12BB, ST12AB;

  nirreps = moinfo.nirreps;

  dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "LIA");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  ST1A = dpd_file2_dot(&T1, &L1);
  dpd_file2_close(&L1);
  dpd_file2_close(&T1);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
  }
  ST1B = dpd_file2_dot(&T1, &L1);
  dpd_file2_close(&L1);
  dpd_file2_close(&T1);

  dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  ST2AA = dpd_buf4_dot(&L2, &T2);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
  }
  ST2BB = dpd_buf4_dot(&L2, &T2);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  }
  ST2AB = dpd_buf4_dot(&L2, &T2);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);

  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1A);
  dpd_file2_mat_rd(&T1A);
  if(params.ref == 0 || params.ref == 1) /** RHF/ROHF **/
    dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  else if(params.ref == 2) /** UHF **/
    dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_file2_mat_init(&T1B);
  dpd_file2_mat_rd(&T1B);

  ST12AA = 0.0;
  dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&L2, h);
    dpd_buf4_mat_irrep_rd(&L2, h);
    for(row=0; row < L2.params->rowtot[h]; row++) {
      i = L2.params->roworb[h][row][0];
      j = L2.params->roworb[h][row][1];
      I = T1A.params->rowidx[i]; Isym = T1A.params->psym[i];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      for(col=0; col < L2.params->coltot[h^L_irr]; col++) {
	a = L2.params->colorb[h^L_irr][col][0];
	b = L2.params->colorb[h^L_irr][col][1];
	A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
	B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];
	if((Isym == Asym) && (Jsym == Bsym))
	  ST12AA += L2.matrix[h][row][col] * 
	    T1A.matrix[Isym][I][A] * T1A.matrix[Jsym][J][B];
	if((Isym == Bsym) && (Jsym == Asym))
	  ST12AA -= L2.matrix[h][row][col] * 
	    T1A.matrix[Isym][I][B] * T1A.matrix[Jsym][J][A];
      }
    }
    dpd_buf4_mat_irrep_close(&L2, h);
  }
  dpd_buf4_close(&L2);

  ST12BB = 0.0;

  if(params.ref == 0 || params.ref == 1)
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
  else if(params.ref == 2)
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&L2, h);
    dpd_buf4_mat_irrep_rd(&L2, h);
    for(row=0; row < L2.params->rowtot[h]; row++) {
      i = L2.params->roworb[h][row][0];
      j = L2.params->roworb[h][row][1];
      I = T1B.params->rowidx[i]; Isym = T1B.params->psym[i];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      for(col=0; col < L2.params->coltot[h^L_irr]; col++) {
	a = L2.params->colorb[h^L_irr][col][0];
	b = L2.params->colorb[h^L_irr][col][1];
	A = T1B.params->colidx[a]; Asym = T1B.params->qsym[a];
	B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
	if((Isym == Asym) && (Jsym == Bsym))
	  ST12BB += L2.matrix[h][row][col] *
	    T1B.matrix[Isym][I][A] * T1B.matrix[Jsym][J][B];
	if((Isym == Bsym) && (Jsym == Asym))
	  ST12BB -= L2.matrix[h][row][col] *
	    T1B.matrix[Isym][I][B] * T1B.matrix[Jsym][J][A];
      }
    }
    dpd_buf4_mat_irrep_close(&L2, h);
  }
  dpd_buf4_close(&L2);

  ST12AB = 0.0;

  if(params.ref == 0 || params.ref == 1)
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  else if(params.ref == 2)
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&L2, h);
    dpd_buf4_mat_irrep_rd(&L2, h);
    for(row=0; row < L2.params->rowtot[h]; row++) {
      i = L2.params->roworb[h][row][0];
      j = L2.params->roworb[h][row][1];
      I = T1A.params->rowidx[i]; Isym = T1A.params->psym[i];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      for(col=0; col < L2.params->coltot[h^L_irr]; col++) {
	a = L2.params->colorb[h^L_irr][col][0];
	b = L2.params->colorb[h^L_irr][col][1];
	A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
	B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
	if((Isym == Asym) && (Jsym == Bsym))
	  ST12AB += L2.matrix[h][row][col] *
	    T1A.matrix[Isym][I][A] * T1B.matrix[Jsym][J][B];
      }
    }
    dpd_buf4_mat_irrep_close(&L2, h);
  }
  dpd_buf4_close(&L2);


  dpd_file2_mat_close(&T1A);
  dpd_file2_close(&T1A);
  dpd_file2_mat_close(&T1B);
  dpd_file2_close(&T1B);

  /*
    fprintf(outfile, "\tST1A = %20.15f\n", ST1A);
    fprintf(outfile, "\tST1B = %20.15f\n", ST1B);
    fprintf(outfile, "\tST2AA = %20.15f\n", ST2AA);
    fprintf(outfile, "\tST2BB = %20.15f\n", ST2BB);
    fprintf(outfile, "\tST2AB = %20.15f\n", ST2AB);
    fprintf(outfile, "\tST12AA = %20.15f\n", ST12AA);
    fprintf(outfile, "\tST12BB = %20.15f\n", ST12BB);
    fprintf(outfile, "\tST12AB = %20.15f\n", ST12AB);
  */

  value = 1.0 - ST1A - ST1B - ST2AA - ST2BB - ST2AB + ST12AA + ST12BB + ST12AB;

  fprintf(outfile, "\tOverlap <L|e^T> = %20.11f\n", value);
}

}} // namespace psi::cclambda
