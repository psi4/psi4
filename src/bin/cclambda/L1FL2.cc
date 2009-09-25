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

/* L2FL2(): Computes the contributions of the Fme HBAR matrix elements
** to the Lambda doubles equations.  These contributions are given in
** spin orbitals as:
**
** L_ij^ab <-- P(ij) P(ab) L_i^a Fjb
**
** where Fjb = fjb + t_n^f <jn||bf>
**
** TDC, July 2002
*/

void L1FL2(int L_irr)
{
  int h, nirreps;
  int row,col;
  int i,j,a,b,I,J,A,B,Isym,Jsym,Asym,Bsym;
  dpdfile2 LIA, Lia, FJB, Fjb, L, F;
  dpdbuf4 newL2;

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&L, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_mat_init(&L);
    dpd_file2_mat_rd(&L);
    dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_mat_init(&F);
    dpd_file2_mat_rd(&F);

    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&newL2, h);
      dpd_buf4_mat_irrep_rd(&newL2, h);

      for(row=0; row < newL2.params->rowtot[h]; row++) {
	i = newL2.params->roworb[h][row][0];
	j = newL2.params->roworb[h][row][1];
	  
	for(col=0; col < newL2.params->coltot[h^L_irr]; col++) {
	  a = newL2.params->colorb[h^L_irr][col][0];
	  b = newL2.params->colorb[h^L_irr][col][1];

	  I = L.params->rowidx[i]; Isym = L.params->psym[i];
	  J = F.params->rowidx[j]; Jsym = F.params->psym[j];
	  A = L.params->colidx[a]; Asym = L.params->qsym[a];
	  B = F.params->colidx[b]; Bsym = F.params->qsym[b];
	  if(((Isym^Asym) == L_irr) && (Jsym == Bsym))
	    newL2.matrix[h][row][col] += (L.matrix[Isym][I][A] * F.matrix[Jsym][J][B]);

	  if((Isym == Asym) && ((Jsym^Bsym) == L_irr))
	    newL2.matrix[h][row][col] += (L.matrix[Jsym][J][B] * F.matrix[Isym][I][A]);
	}
      }

      dpd_buf4_mat_irrep_wrt(&newL2, h);
      dpd_buf4_mat_irrep_close(&newL2, h);
      
    }

    dpd_buf4_close(&newL2);

    dpd_file2_mat_close(&F);
    dpd_file2_close(&F);
    dpd_file2_mat_close(&L);
    dpd_file2_close(&L);

  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_mat_init(&LIA);
    dpd_file2_mat_rd(&LIA);
    dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "Lia");
    dpd_file2_mat_init(&Lia);
    dpd_file2_mat_rd(&Lia);
    dpd_file2_init(&FJB, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_mat_init(&FJB);
    dpd_file2_mat_rd(&FJB);
    dpd_file2_init(&Fjb, CC_OEI, 0, 0, 1, "Fme");
    dpd_file2_mat_init(&Fjb);
    dpd_file2_mat_rd(&Fjb);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_mat_init(&LIA);
    dpd_file2_mat_rd(&LIA);
    dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 2, 3, "Lia");
    dpd_file2_mat_init(&Lia);
    dpd_file2_mat_rd(&Lia);
    dpd_file2_init(&FJB, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_mat_init(&FJB);
    dpd_file2_mat_rd(&FJB);
    dpd_file2_init(&Fjb, CC_OEI, 0, 2, 3, "Fme");
    dpd_file2_mat_init(&Fjb);
    dpd_file2_mat_rd(&Fjb);
  
  }

  if(params.ref == 1) /** RHF/ROHF **/
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
  else if(params.ref == 2) /** UHF **/
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");

  if(params.ref == 1 || params.ref == 2) {
    /* loop over row irreps of LIJAB */
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&newL2, h);
      dpd_buf4_mat_irrep_rd(&newL2, h);

      /* loop over rows of irrep of LIJAB */
      for(row=0; row < newL2.params->rowtot[h]; row++) {
	i = newL2.params->roworb[h][row][0];
	j = newL2.params->roworb[h][row][1];
	  
	/* loop over cols of irrep of LIJAB */
	for(col=0; col < newL2.params->coltot[h^L_irr]; col++) {
	  a = newL2.params->colorb[h^L_irr][col][0];
	  b = newL2.params->colorb[h^L_irr][col][1];

	  I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
	  J = FJB.params->rowidx[j]; Jsym = FJB.params->psym[j];
	  A = LIA.params->colidx[a]; Asym = LIA.params->qsym[a];
	  B = FJB.params->colidx[b]; Bsym = FJB.params->qsym[b];

	  if( ((Isym^Asym) == L_irr) && (Jsym == Bsym) )
	    newL2.matrix[h][row][col] += (LIA.matrix[Isym][I][A] *
					  FJB.matrix[Jsym][J][B]);

	  J = LIA.params->rowidx[j]; Jsym = LIA.params->psym[j];
	  I = FJB.params->rowidx[i]; Isym = FJB.params->psym[i];

	  if( (Isym == Asym) && ((Jsym^Bsym) == L_irr) )
	    newL2.matrix[h][row][col] += (LIA.matrix[Jsym][J][B] *
					  FJB.matrix[Isym][I][A]);

	  I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
	  J = FJB.params->rowidx[j]; Jsym = FJB.params->psym[j];
	  B = LIA.params->colidx[b]; Bsym = LIA.params->qsym[b];
	  A = FJB.params->colidx[a]; Asym = FJB.params->qsym[a];

	  if( ((Jsym^Asym) == L_irr) && (Isym == Bsym))
	    newL2.matrix[h][row][col] -= (LIA.matrix[Jsym][J][A] *
					  FJB.matrix[Isym][I][B]);

	  J = LIA.params->rowidx[j]; Jsym = LIA.params->psym[j];
	  I = FJB.params->rowidx[i]; Isym = FJB.params->psym[i];

	  if( (Jsym == Asym) && ((Isym^Bsym) == L_irr) )
	    newL2.matrix[h][row][col] -= (LIA.matrix[Isym][I][B] *
					  FJB.matrix[Jsym][J][A]);
	}
      }

      dpd_buf4_mat_irrep_wrt(&newL2, h);
      dpd_buf4_mat_irrep_close(&newL2, h);
      
    }
    dpd_buf4_close(&newL2);
  }

  if(params.ref == 1) /** RHF/ROHF **/
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
  else if(params.ref == 2) /** UHF **/
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");

  if(params.ref == 1 || params.ref == 2) {
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&newL2, h);
      dpd_buf4_mat_irrep_rd(&newL2, h);

      for(row=0; row < newL2.params->rowtot[h]; row++) {
	i = newL2.params->roworb[h][row][0];
	j = newL2.params->roworb[h][row][1];
	  
	for(col=0; col < newL2.params->coltot[h^L_irr]; col++) {
	  a = newL2.params->colorb[h^L_irr][col][0];
	  b = newL2.params->colorb[h^L_irr][col][1];

	  I = Lia.params->rowidx[i]; Isym = Lia.params->psym[i];
	  J = Fjb.params->rowidx[j]; Jsym = Fjb.params->psym[j];
	  A = Lia.params->colidx[a]; Asym = Lia.params->qsym[a];
	  B = Fjb.params->colidx[b]; Bsym = Fjb.params->qsym[b];

	  if(((Isym^Asym) == L_irr) && (Jsym == Bsym))
	    newL2.matrix[h][row][col] += (Lia.matrix[Isym][I][A] *
					  Fjb.matrix[Jsym][J][B]);

	  J = Lia.params->rowidx[j]; Jsym = Lia.params->psym[j];
	  I = Fjb.params->rowidx[i]; Isym = Fjb.params->psym[i];

	  if((Isym == Asym) && ((Jsym^Bsym) == L_irr))
	    newL2.matrix[h][row][col] += (Lia.matrix[Jsym][J][B] *
					  Fjb.matrix[Isym][I][A]);

	  I = Lia.params->rowidx[i]; Isym = Lia.params->psym[i];
	  J = Fjb.params->rowidx[j]; Jsym = Fjb.params->psym[j];
	  B = Lia.params->colidx[b]; Bsym = Lia.params->qsym[b];
	  A = Fjb.params->colidx[a]; Asym = Fjb.params->qsym[a];

	  if(((Jsym^Asym) == L_irr) && (Isym == Bsym))
	    newL2.matrix[h][row][col] -= (Lia.matrix[Jsym][J][A] *
					  Fjb.matrix[Isym][I][B]);

	  J = Lia.params->rowidx[j]; Jsym = Lia.params->psym[j];
	  I = Fjb.params->rowidx[i]; Isym = Fjb.params->psym[i];

	  if((Jsym == Asym) && ((Isym^Bsym) == L_irr))
	    newL2.matrix[h][row][col] -= (Lia.matrix[Isym][I][B] *
					  Fjb.matrix[Jsym][J][A]);
	}
      }

      dpd_buf4_mat_irrep_wrt(&newL2, h);
      dpd_buf4_mat_irrep_close(&newL2, h);
      
    }
    dpd_buf4_close(&newL2);
  }

  if(params.ref == 1) /** RHF/ROHF **/
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
  else if(params.ref == 2) /** UHF **/
    dpd_buf4_init(&newL2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");

  if(params.ref == 1 || params.ref == 2) {
    for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&newL2, h);
      dpd_buf4_mat_irrep_rd(&newL2, h);

      for(row=0; row < newL2.params->rowtot[h]; row++) {
	i = newL2.params->roworb[h][row][0];
	j = newL2.params->roworb[h][row][1];
	  
	for(col=0; col < newL2.params->coltot[h^L_irr]; col++) {
	  a = newL2.params->colorb[h^L_irr][col][0];
	  b = newL2.params->colorb[h^L_irr][col][1];

	  I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
	  J = Fjb.params->rowidx[j]; Jsym = Fjb.params->psym[j];
	  A = LIA.params->colidx[a]; Asym = LIA.params->qsym[a];
	  B = Fjb.params->colidx[b]; Bsym = Fjb.params->qsym[b];

	  if(((Isym^Asym) == L_irr) && (Jsym == Bsym))
	    newL2.matrix[h][row][col] += (LIA.matrix[Isym][I][A] *
					  Fjb.matrix[Jsym][J][B]);

	  J = Lia.params->rowidx[j]; Jsym = Lia.params->psym[j];
	  I = FJB.params->rowidx[i]; Isym = FJB.params->psym[i];
	  B = Lia.params->colidx[b]; Bsym = Lia.params->qsym[b];
	  A = FJB.params->colidx[a]; Asym = FJB.params->qsym[a];

	  if((Isym == Asym) && ((Jsym^Bsym) == L_irr))
	    newL2.matrix[h][row][col] += (Lia.matrix[Jsym][J][B] *
					  FJB.matrix[Isym][I][A]);
	}
      }

      dpd_buf4_mat_irrep_wrt(&newL2, h);
      dpd_buf4_mat_irrep_close(&newL2, h);
      
    }
  }

  if(params.ref == 1 || params.ref == 2) {
    dpd_buf4_close(&newL2);

    dpd_file2_mat_close(&FJB);
    dpd_file2_close(&FJB);
    dpd_file2_mat_close(&Fjb);
    dpd_file2_close(&Fjb);
    dpd_file2_mat_close(&LIA);
    dpd_file2_close(&LIA);
    dpd_file2_mat_close(&Lia);
    dpd_file2_close(&Lia);
  }
}

}} // namespace psi::cclambda
