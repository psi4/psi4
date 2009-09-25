/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void cc3_l3l1(void)
{
  dpdfile2 L1, D1, L1new;
  dpdbuf4 Z, W;
  int nirreps, Gde, Gg, Gi, Ga;
  int de, ig, ag;
  int nrows, ncols, nlinks;

  nirreps = moinfo.nirreps;

  dpd_buf4_init(&Z, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZIGDE");
  dpd_buf4_sort(&Z, CC3_MISC, rspq, 5, 10, "CC3 ZIGDE (DE,IG)");
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZIgDe");
  dpd_buf4_sort(&Z, CC3_MISC, rspq, 5, 10, "CC3 ZIgDe (De,Ig)");
  dpd_buf4_close(&Z);

  dpd_file2_init(&L1, CC3_MISC, 0, 0, 1, "CC3 LIA");
  dpd_file2_mat_init(&L1);

  dpd_buf4_init(&W, CC3_HET1, 0, 5, 5, 7, 7, 0, "CC3 WABEF");
  dpd_buf4_init(&Z, CC3_MISC, 0, 5, 10, 5, 10, 0, "CC3 ZIGDE (DE,IG)");
  for(Gde=0; Gde < nirreps; Gde++) {
    if(Z.params->coltot[Gde] && W.params->coltot[Gde]) {
      Z.matrix[Gde] = dpd_block_matrix(1, Z.params->coltot[Gde]);
      W.matrix[Gde] = dpd_block_matrix(1, W.params->coltot[Gde]);
      for(de=0; de < Z.params->rowtot[Gde]; de++) {
	dpd_buf4_mat_irrep_rd_block(&W, Gde, de, 1);
	dpd_buf4_mat_irrep_rd_block(&Z, Gde, de, 1);

	for(Gg=0; Gg < nirreps; Gg++) {
	  Ga = Gi = Gg ^ Gde; /* totally symmetric */
	  nrows = L1.params->rowtot[Gi];
	  ncols = L1.params->coltot[Gi];
	  nlinks = Z.params->spi[Gg];

	  ig = Z.col_offset[Gde][Gi];
	  ag = W.col_offset[Gde][Ga];

	  if(nrows && ncols && nlinks)
	    C_DGEMM('n', 't', nrows, ncols, nlinks, 0.5, &(Z.matrix[Gde][0][ig]), nlinks,
		    &(W.matrix[Gde][0][ag]), nlinks, 1.0, L1.matrix[Gi][0], ncols);
	}
      }
      dpd_free_block(Z.matrix[Gde], 1, Z.params->coltot[Gde]);
      dpd_free_block(W.matrix[Gde], 1, W.params->coltot[Gde]);
    }
  }
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC3_HET1, 0, 5, 5, 5, 5, 0, "CC3 WAbEf");
  dpd_buf4_init(&Z, CC3_MISC, 0, 5, 10, 5, 10, 0, "CC3 ZIgDe (De,Ig)");
  for(Gde=0; Gde < nirreps; Gde++) {
    if(Z.params->coltot[Gde] && W.params->coltot[Gde]) {
      Z.matrix[Gde] = dpd_block_matrix(1, Z.params->coltot[Gde]);
      W.matrix[Gde] = dpd_block_matrix(1, W.params->coltot[Gde]);
      for(de=0; de < Z.params->rowtot[Gde]; de++) {
	dpd_buf4_mat_irrep_rd_block(&W, Gde, de, 1);
	dpd_buf4_mat_irrep_rd_block(&Z, Gde, de, 1);

	for(Gg=0; Gg < nirreps; Gg++) {
	  Ga = Gi = Gg ^ Gde; /* totally symmetric */
	  nrows = L1.params->rowtot[Gi];
	  ncols = L1.params->coltot[Gi];
	  nlinks = Z.params->spi[Gg];

	  ig = Z.col_offset[Gde][Gi];
	  ag = W.col_offset[Gde][Ga];

	  if(nrows && ncols && nlinks)
	    C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, &(Z.matrix[Gde][0][ig]), nlinks,
		    &(W.matrix[Gde][0][ag]), nlinks, 1.0, L1.matrix[Gi][0], ncols);
	}
      }
      dpd_free_block(Z.matrix[Gde], 1, Z.params->coltot[Gde]);
      dpd_free_block(W.matrix[Gde], 1, W.params->coltot[Gde]);
    }
  }
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  dpd_file2_mat_wrt(&L1);
  dpd_file2_mat_close(&L1);

  /* Wmbej --> L1 */

  dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMBEJ (ME,JB)");
  dpd_buf4_sort(&W, CC3_HET1, psrq, 10, 10, "CC3 WMBEJ (MB,JE)");
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
  dpd_buf4_sort(&W, CC3_HET1, psrq, 10, 10, "CC3 WMbEj (Mb,jE)");
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
  dpd_buf4_sort(&W, CC3_HET1, psrq, 10, 10, "CC3 WMbeJ (Mb,Je)");
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMBEJ (MB,JE)");
  dpd_buf4_init(&Z, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZDMAE (MD,AE)");
  dpd_contract442(&W, &Z, &L1, 0, 2, -1, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (Mb,jE)");
  dpd_buf4_init(&Z, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZDmAe (mD,Ae)");
  dpd_contract442(&W, &Z, &L1, 0, 2, -1, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Mb,Je)");
  dpd_buf4_init(&Z, CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZdMAe (Md,Ae)");
  dpd_contract442(&W, &Z, &L1, 0, 2, -1, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC3_HET1, 0, 10, 11, 10, 11, 0, "CC3 WMBEJ (MB,EJ)");
  dpd_buf4_init(&Z, CC3_MISC, 0, 0, 10, 0, 10, 0, "CC3 ZIMLE");
  dpd_contract442(&Z, &W, &L1, 0, 2, 1, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC3_HET1, 0, 10, 11, 10, 11, 0, "CC3 WMbEj (Mb,Ej)");
  dpd_buf4_init(&Z, CC3_MISC, 0, 0, 10, 0, 10, 0, "CC3 ZImLe");
  dpd_contract442(&Z, &W, &L1, 0, 2, 1, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC3_HET1, 0, 10, 11, 10, 11, 0, "CC3 WmBEj (mB,Ej)");
  dpd_buf4_init(&Z, CC3_MISC, 0, 0, 10, 0, 10, 0, "CC3 ZImlE");
  dpd_contract442(&Z, &W, &L1, 0, 2, 1, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  /* Wmnij -> L1 */

  dpd_buf4_init(&W, CC3_HET1, 0, 0, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
  dpd_buf4_init(&Z, CC3_MISC, 0, 0, 11, 0, 11, 0, "CC3 ZLMAO");
  dpd_contract442(&W, &Z, &L1, 0, 2, -0.5, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  dpd_buf4_init(&W, CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
  dpd_buf4_init(&Z, CC3_MISC, 0, 0, 11, 0, 11, 0, "CC3 ZLmAo");
  dpd_contract442(&W, &Z, &L1, 0, 2, -1, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);

  dpd_file2_init(&D1, CC_DENOM, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&D1, &L1);
  dpd_file2_close(&D1);
  dpd_file2_init(&L1new, CC_LAMBDA, 0, 0, 1, "New LIA");
  dpd_file2_axpy(&L1, &L1new, 1, 0);
  dpd_file2_copy(&L1new, CC_LAMBDA, "New Lia");
  dpd_file2_close(&L1new);
  dpd_file2_close(&L1);

}

}} // namespace psi::cclambda
