/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void denom1(dpdfile2 *X1, double omega);
void local_filter_T1(dpdfile2 *T1);

void X1_build(const char *pert, int irrep, double omega)
{
  dpdfile2 F, X1, X1new;
  dpdbuf4 W, X2;
  char lbl[32];
  int Gam, Gef, Gim, Gi, Ga, Gm, nrows, ncols, A, a, am;

  sprintf(lbl, "%sBAR_IA", pert);
  dpd_file2_init(&X1new, CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
  dpd_file2_copy(&X1new, CC_OEI, lbl);
  dpd_file2_close(&X1new);
  dpd_file2_init(&X1new, CC_OEI, irrep, 0, 1, lbl);

  /*** S-S ***/

  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);

  dpd_file2_axpy(&X1, &X1new, -omega, 0);

  dpd_file2_init(&F, CC_OEI, 0, 1, 1, "FAE");
  dpd_contract222(&X1, &F, &X1new, 0, 0, 1, 1);
  dpd_file2_close(&F);

  dpd_file2_init(&F, CC_OEI, 0, 0, 0, "FMI");
  dpd_contract222(&F, &X1, &X1new, 1, 1, -1, 1);
  dpd_file2_close(&F);

  dpd_buf4_init(&W, CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(jb,ME) + W(Jb,Me)");
  dpd_contract422(&W, &X1, &X1new, 0, 0, 1, 1);
  dpd_buf4_close(&W);

  dpd_file2_close(&X1);


  /*** S-D ***/

  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_dot24(&F, &X2, &X1new, 0, 0, 1, 1);
  dpd_buf4_close(&X2);
  dpd_file2_close(&F);

  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  /*  dpd_contract442(&X2, &W, &X1new, 0, 0, 1, 1); */
  /* ooc code below added 7/28/05, -TDC */
  dpd_file2_mat_init(&X1new);
  dpd_file2_mat_rd(&X1new);
  for(Gam=0; Gam < moinfo.nirreps; Gam++) {
    Gef = Gam; /* W is totally symmetric */
    Gim = Gef ^ irrep;

    dpd_buf4_mat_irrep_init(&X2, Gim);
    dpd_buf4_mat_irrep_rd(&X2, Gim);
    dpd_buf4_mat_irrep_shift13(&X2, Gim);

    for(Gi=0; Gi < moinfo.nirreps; Gi++) {
      Ga = Gi ^ irrep;
      Gm = Ga ^ Gam;

      W.matrix[Gam] = dpd_block_matrix(moinfo.occpi[Gm], W.params->coltot[Gef]);

      nrows = moinfo.occpi[Gi];
      ncols = moinfo.occpi[Gm] * W.params->coltot[Gef];

      for(A=0; A < moinfo.virtpi[Ga]; A++) {
	a = moinfo.vir_off[Ga] + A;
	am = W.row_offset[Gam][a];

	dpd_buf4_mat_irrep_rd_block(&W, Gam, am, moinfo.occpi[Gm]);

	if(nrows && ncols && moinfo.virtpi[Ga])
	  C_DGEMV('n',nrows,ncols,1,X2.shift.matrix[Gim][Gi][0],ncols,W.matrix[Gam][0],1,1,
		  &(X1new.matrix[Gi][0][A]), moinfo.virtpi[Ga]);
      }
      dpd_free_block(W.matrix[Gam], moinfo.occpi[Gm], W.params->coltot[Gef]);
    }

    dpd_buf4_mat_irrep_close(&X2, Gim);
  }
  dpd_file2_mat_wrt(&X1new);
  dpd_file2_mat_close(&X1new);
  dpd_buf4_close(&W);
  dpd_buf4_close(&X2);

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
  dpd_contract442(&W, &X2, &X1new, 3, 3, 1, 1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&X2);

  if(params.local && local.filter_singles) local_filter_T1(&X1new);
  else denom1(&X1new, omega);
  dpd_file2_close(&X1new);
}

}} // namespace psi::ccresponse
