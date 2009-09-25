/*! \file
    \ingroup CCRESPONSE
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

void cc2_X1_build(const char *pert, int irrep, double omega)
{
  int GX2, GX1, GW, Gam, Gim, Gef, Ga, Gi, Gm;
  int a, A, i, I, num_m, nlinks, length;
  dpdfile2 F, X1, X1new, Xme, Zia;
  dpdbuf4 W, X2, D, T2;
  char lbl[32];

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

  dpd_buf4_init(&W, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 2 W(jb,ME) + W(Jb,Me)");
  dpd_contract422(&W, &X1, &X1new, 0, 0, 1, 1);
  dpd_buf4_close(&W);

  sprintf(lbl, "X_%s_ME", pert);
  dpd_file2_init(&Xme, CC_OEI, irrep, 0, 1, lbl);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
  dpd_contract422(&D, &X1, &Xme, 0, 0, 1, 0);
  dpd_buf4_close(&D);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  dpd_contract422(&T2, &Xme, &X1new, 0, 0, 1, 1);
  dpd_buf4_close(&T2);
  dpd_file2_close(&Xme);

  dpd_file2_close(&X1);


  /*** S-D ***/

  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_dot24(&F, &X2, &X1new, 0, 0, 1, 1);
  dpd_file2_close(&F);

  /** begin out of core contract442 **/
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  GW = W.file.my_irrep;
  GX2 = X2.file.my_irrep;
  GX1 = X1new.my_irrep;

  dpd_file2_mat_init(&X1new);
  dpd_file2_mat_rd(&X1new);

  for(Gam=0; Gam < moinfo.nirreps; Gam++) {
    Gef = Gam^GW;
    Gim = Gef^GX2;

    for(Gi=0; Gi < moinfo.nirreps; Gi++) {
      Gm = Gef^Gi^GX2;
      Ga = Gi^GX1;

      num_m = W.params->qpi[Gm];
      dpd_buf4_mat_irrep_init_block(&W, Gam, num_m);
      dpd_buf4_mat_irrep_init_block(&X2, Gim, num_m);

      nlinks = W.params->coltot[Gam] * num_m;
      length = X1new.params->rowtot[Gi] * X1new.params->coltot[Ga];

      if(nlinks && length) {
	for(i=0; i < X1new.params->rowtot[Gi]; i++) {

	  I = X2.params->poff[Gi] + i;
	  dpd_buf4_mat_irrep_rd_block(&X2, Gim, X2.row_offset[Gim][I], num_m);

	  for(a=0; a < X1new.params->coltot[Ga]; a++) {

	    A = W.params->poff[Ga] + a;
	    dpd_buf4_mat_irrep_rd_block(&W, Gam, W.row_offset[Gam][A], num_m);

  	    X1new.matrix[Gi][i][a] += C_DDOT(nlinks, X2.matrix[Gim][0], 1,
  					     W.matrix[Gam][0], 1);
  	  }
  	}
      }
      dpd_buf4_mat_irrep_close_block(&X2, Gim, num_m);
      dpd_buf4_mat_irrep_close_block(&W, Gam, num_m);
    }
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
