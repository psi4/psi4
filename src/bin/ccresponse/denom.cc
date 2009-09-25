/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void denom1(dpdfile2 *X1, double omega)
{
  int nirreps, h, irrep;
  int i, a;
  int *occpi, *virtpi;
  dpdfile2 FAE, FMI;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  virtpi = moinfo.virtpi;

  irrep = X1->my_irrep;

  if((!strcmp(params.wfn,"CC2")) || (!strcmp(params.wfn,"EOM_CC2"))) {
    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_mat_init(&FMI);
    dpd_file2_mat_rd(&FMI);

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_mat_init(&FAE);
    dpd_file2_mat_rd(&FAE);
  }
  else {
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_mat_init(&FAE);
    dpd_file2_mat_rd(&FAE);

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_mat_init(&FMI);
    dpd_file2_mat_rd(&FMI);
  }

  dpd_file2_mat_init(X1);
  dpd_file2_mat_rd(X1);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) 
      for(a=0; a < virtpi[h^irrep]; a++)
	X1->matrix[h][i][a] /= (FMI.matrix[h][i][i] - FAE.matrix[h^irrep][a][a] + omega);
  }
  dpd_file2_mat_wrt(X1);
  dpd_file2_mat_close(X1);

  dpd_file2_mat_close(&FAE);
  dpd_file2_mat_close(&FMI);
  dpd_file2_close(&FAE);
  dpd_file2_close(&FMI);
}

void denom2(dpdbuf4 *X2, double omega)
{
  int nirreps, h, row, col, irrep;
  int i, j, I, J, a, b, A, B, isym, jsym, asym, bsym;
  dpdfile2 FAE, FMI;

  nirreps = moinfo.nirreps;
  irrep = X2->file.my_irrep;

  if((!strcmp(params.wfn,"CC2")) || (!strcmp(params.wfn,"EOM_CC2"))) {
    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_mat_init(&FMI);
    dpd_file2_mat_rd(&FMI);

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_mat_init(&FAE);
    dpd_file2_mat_rd(&FAE);
  }
  else {
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_mat_init(&FAE);
    dpd_file2_mat_rd(&FAE);

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_mat_init(&FMI);
    dpd_file2_mat_rd(&FMI);
  }

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(X2, h);
    dpd_buf4_mat_irrep_rd(X2, h);

    for(row=0; row < X2->params->rowtot[h]; row++) {

      i = X2->params->roworb[h][row][0];
      j = X2->params->roworb[h][row][1];
      isym = X2->params->psym[i];
      jsym = X2->params->qsym[j];
      I = i - moinfo.occ_off[isym];
      J = j - moinfo.occ_off[jsym];

      for(col=0; col < X2->params->coltot[h^irrep]; col++) {

	a = X2->params->colorb[h^irrep][col][0];
	b = X2->params->colorb[h^irrep][col][1];
	asym = X2->params->rsym[a];
	bsym = X2->params->ssym[b];
	A = a - moinfo.vir_off[asym];
	B = b - moinfo.vir_off[bsym];

	X2->matrix[h][row][col] /= (FMI.matrix[isym][I][I] + FMI.matrix[jsym][J][J] - 
				    FAE.matrix[asym][A][A] - FAE.matrix[bsym][B][B] + omega);

      }
    }

    dpd_buf4_mat_irrep_wrt(X2, h);
    dpd_buf4_mat_irrep_close(X2, h);
  }

  dpd_file2_mat_close(&FAE);
  dpd_file2_mat_close(&FMI);
  dpd_file2_close(&FAE);
  dpd_file2_close(&FMI);
}

}} // namespace psi::ccresponse
