/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double LHX1Y2(const char *pert_x, int irrep_x, double omega_x, 
	      const char *pert_y, int irrep_y, double omega_y)
{
  dpdfile2 z, z1, X1, l1, F;
  dpdbuf4 Z, Z1, Z2, I, Y2, L2, W;
  char lbl[32];
  double polar;
  int nirreps, Gbm, Gef, Gjf, Ge, Gf, Gj, bm, ef, jf;
  int Gfe, Gb, Gm, b, m, B, M, e, f, fe, nrows, ncols, nlinks;
  double *X;

  sprintf(lbl, "Z_%s_MI", pert_y);
  dpd_file2_init(&z1, CC_TMP0, irrep_y, 0, 0, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract442(&I, &Y2, &z1, 0, 0, 1, 0);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&I);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&z1, &X1, &z, 1, 1, -1, 0);
  dpd_file2_close(&X1);
  dpd_file2_close(&z1);
  dpd_file2_close(&z);

  sprintf(lbl, "Z_%s_AE", pert_y);
  dpd_file2_init(&z1, CC_TMP0, irrep_y, 1, 1, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract442(&Y2, &I, &z1, 3, 3, -1, 0);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&I);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract222(&X1, &z1, &z, 0, 0, 1, 1);
  dpd_file2_close(&X1);
  dpd_file2_close(&z1);
  dpd_file2_close(&z);


  sprintf(lbl, "Z_%s_ME", pert_x);
  dpd_file2_init(&z1, CC_TMP0, irrep_x, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_dot24(&X1, &I, &z1, 0, 0, 1, 0);
  dpd_buf4_close(&I);
  dpd_file2_close(&X1);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_dot24(&z1, &Y2, &z, 0, 0, 1, 1);
  dpd_buf4_close(&Y2);
  dpd_file2_close(&z1);
  dpd_file2_close(&z);

  dpd_file2_init(&z, CC_TMP0, 0, 0, 1, "Z(I,A) Final");
  dpd_file2_init(&l1, CC_LAMPS, 0, 0, 1, "LIA 0 -1");
  polar = 2.0 * dpd_file2_dot(&z, &l1);
  dpd_file2_close(&l1);
  dpd_file2_close(&z);

//  fprintf(outfile, "L(1)HX1Y2 = %20.12f\n", polar);


  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_scm(&Z, 0);

  sprintf(lbl, "Z_%s_MI", pert_x);
  dpd_file2_init(&z1, CC_TMP0, irrep_x, 0, 0, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  dpd_contract222(&F, &X1, &z1, 0, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_file2_close(&X1);

  dpd_buf4_init(&Z1, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract244(&z1, &Y2, &Z1, 0, 0, 0, 1, 0);
  dpd_buf4_close(&Y2);
  dpd_file2_close(&z1);
  dpd_buf4_axpy(&Z1, &Z, -1);
  dpd_buf4_sort(&Z1, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z1, &Z, -1);
  dpd_buf4_close(&Z1);

  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");

  sprintf(lbl, "Z_%s_AE", pert_x);
  dpd_file2_init(&z1, CC_TMP0, irrep_x, 1, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  dpd_contract222(&X1, &F, &z1, 1, 1, -1, 0);
  dpd_file2_close(&F);
  dpd_file2_close(&X1);

  dpd_buf4_init(&Z1, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract424(&Y2, &z1, &Z1, 3, 1, 0, 1, 0);
  dpd_buf4_close(&Y2);
  dpd_file2_close(&z1);
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_sort(&Z1, CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z1);

  dpd_buf4_close(&Z);

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);

  sprintf(lbl, "Z_%s_MbEj (bM,jE)", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 11, 10, 11, 10, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  /*  dpd_contract244(&X1, &W, &Z1, 1, 2, 1, 1, 0); */
  /* add out-of-core algorithm here */
  /* Z(bm,jf) <-- X1(j,e) * W(bm,ef) */
  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);
  for(Gbm=0; Gbm < moinfo.nirreps; Gbm++) {
    Gef = Gbm;  Gjf = Gbm^irrep_x;

    dpd_buf4_mat_irrep_init(&Z1, Gbm);
    dpd_buf4_mat_irrep_row_init(&W, Gbm);
    for(bm=0; bm < W.params->rowtot[Gbm]; bm++) {
      dpd_buf4_mat_irrep_row_rd(&W, Gbm, bm);

      for(Gj=0; Gj < moinfo.nirreps; Gj++) {
	Gf = Gjf^Gj;  Ge = Gef^Gf;

	ef = W.col_offset[Gbm][Ge];
	jf = Z1.col_offset[Gbm][Gj];

	nrows = moinfo.occpi[Gj];
	ncols = moinfo.virtpi[Gf];
	nlinks = moinfo.virtpi[Ge];

	if(nrows && ncols && nlinks)
	  C_DGEMM('n','n', nrows, ncols, nlinks, 1.0, &(X1.matrix[Gj][0][0]), nlinks,
		  &(W.matrix[Gbm][0][ef]), ncols, 0.0, &(Z1.matrix[Gbm][bm][jf]), ncols);
      }
    }
    dpd_buf4_mat_irrep_row_close(&W, Gbm);
    dpd_buf4_mat_irrep_wrt(&Z1, Gbm);
    dpd_buf4_mat_irrep_close(&Z1, Gbm);

  }
  dpd_file2_mat_close(&X1);
  /* out-of-core algorithm done */
  dpd_buf4_close(&W);

  sprintf(lbl, "Z_%s_MbEj (ME,jb)", pert_x);
  dpd_buf4_sort(&Z1, CC_TMP0, qsrp, 10, 10, lbl);
  dpd_buf4_close(&Z1);

  sprintf(lbl, "Z_%s_WMbEj (bM,Ej)", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 11, 11, 11, 11, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_contract244(&X1, &W, &Z1, 0, 0, 0, -1, 0);
  dpd_buf4_close(&W);
  sprintf(lbl, "Z_%s_MbEj (ME,jb)", pert_x);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qrsp, 10, 10, lbl, 1);
  dpd_buf4_close(&Z1);

  sprintf(lbl, "Z_%s_MbeJ (Mb,eJ)", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 10, 11, 10, 11, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_contract424(&W, &X1, &Z1, 1, 0, 1, 1, 0);
  dpd_buf4_close(&W);
  sprintf(lbl, "Z_%s_MbeJ (Me,Jb)", pert_x);
  dpd_buf4_sort(&Z1, CC_TMP0, prsq, 10, 10, lbl);
  dpd_buf4_close(&Z1);

  sprintf(lbl, "Z_%s_MbeJ (bM,eJ)", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 11, 11, 11, 11, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_contract424(&W, &X1, &Z1, 3, 1, 0, -1, 0);
  dpd_buf4_close(&W);
  sprintf(lbl, "Z_%s_MbeJ (Me,Jb)", pert_x);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qrsp, 10, 10, lbl, 1);
  dpd_buf4_close(&Z1);

  dpd_file2_close(&X1);

  sprintf(lbl, "Z_%s_MbEj (ME,jb)", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "Z_%s_(2MbEj+MbeJ) (ME,JB)", pert_x);
  dpd_buf4_scmcopy(&Z1, CC_TMP0, lbl, 2);
  dpd_buf4_close(&Z1);
  sprintf(lbl, "Z_%s_MbeJ (Me,Jb)", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "Z_%s_(2MbEj+MbeJ) (ME,JB)", pert_x);
  dpd_buf4_init(&Z2, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  dpd_buf4_axpy(&Z1, &Z2, 1);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) I");
  sprintf(lbl, "X_%s_(2IAjb-IbjA) (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "Z_%s_(2MbEj+MbeJ) (ME,JB)", pert_x);
  dpd_buf4_init(&Z, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  dpd_contract444(&Y2, &Z, &Z1, 0, 1, 0.5, 0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) Ib");
  sprintf(lbl, "Z_%s_MbeJ (Me,Jb)", pert_x);
  dpd_buf4_init(&Z, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
  dpd_contract444(&Y2, &Z, &Z1, 0, 1, 1, 0);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) I");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) Ib");
  dpd_buf4_axpy(&Z2, &Z1, 0.5);
  dpd_buf4_sort(&Z2, CC_TMP0, psrq, 10, 10, "Z(IA,jb) III");
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) I");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) III");
  dpd_buf4_axpy(&Z2, &Z1, 1);
  dpd_buf4_close(&Z2);
  dpd_buf4_sort(&Z1, CC_TMP0, prqs, 0, 5, "Z(Ij,Ab) I+III");
  dpd_buf4_close(&Z1);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) I+III");
  dpd_buf4_sort(&Z, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) II+IV");
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) I+III");
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) II+IV");
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&Z);

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);

  sprintf(lbl, "Z_%s_mj" , pert_x);
  dpd_file2_init(&z, CC_TMP0, irrep_x, 0, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
  dpd_dot23(&X1, &W, &z, 0, 0, 1, 0);
  dpd_buf4_close(&W);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  dpd_contract424(&Y2, &z, &Z1, 1, 0, 1, -1, 0);
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z1);
  dpd_file2_close(&z);

/*** marker ***/
/*** I suspect a bug in this term.  H2O2/STO-3G or DZ yields different results
 * with symmetry on and off here.  The discrepancy disappears for larger
 * basis sets, however. Note also that a similar bug may exist in the X1 or
 * X2 code, based on comparisons with a spin-adapted toy code.  
 * -TDC 5 May 2009 */
  sprintf(lbl, "Z_%s_AE", pert_x);
  dpd_file2_init(&z, CC_TMP0, irrep_x, 1, 1, lbl);

/*   dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); */
/*   dpd_dot24(&X1, &W, &z, 0, 0, 1, 0); */
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_file2_scm(&z, 0);
  dpd_file2_mat_init(&z);
  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);
  for(Gbm=0; Gbm < moinfo.nirreps; Gbm++) {
    Gfe = Gbm; /* W is totally symmetric */
    dpd_buf4_mat_irrep_row_init(&W, Gbm);
    X = init_array(W.params->coltot[Gfe]);
    for(bm=0; bm < W.params->rowtot[Gbm]; bm++) {
      dpd_buf4_mat_irrep_row_rd(&W, Gbm, bm);
      b = W.params->roworb[Gbm][bm][0];
      m = W.params->roworb[Gbm][bm][1];
      Gb = W.params->psym[b];
      Gm = Gbm ^ Gb;
      Ge = Gm ^ irrep_x;
      Gf = Gfe ^ Ge;
      B = b - moinfo.vir_off[Gb];
      M = m - moinfo.occ_off[Gm];
      zero_arr(X, W.params->coltot[Gfe]);
      for(fe=0; fe < W.params->coltot[Gfe]; fe++) {
	f = W.params->colorb[Gfe][fe][0];
	e = W.params->colorb[Gfe][fe][1];
	ef = W.params->colidx[e][f];
	X[fe] = 2.0 * W.matrix[Gbm][0][fe] - W.matrix[Gbm][0][ef];
      }
      nrows = moinfo.virtpi[Gf];
      ncols = moinfo.virtpi[Ge];
      if(nrows & ncols) {
	C_DGEMV('n',nrows,ncols,1,&X[W.col_offset[Gfe][Gf]],ncols,
		X1.matrix[Gm][M],1,1,z.matrix[Gb][B],1);
      }
    }
    free(X);
    dpd_buf4_mat_irrep_row_close(&W, Gbm);
  }
  dpd_file2_mat_close(&X1);
  dpd_file2_mat_wrt(&z);
  dpd_file2_mat_close(&z);
  dpd_buf4_close(&W);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) temp");
  dpd_contract424(&Y2, &z, &Z1, 3, 1, 0, 1, 0);
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z1);
  dpd_file2_close(&z);

  dpd_file2_close(&X1);
  dpd_buf4_close(&Y2);
/*** Marker ***/

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  sprintf(lbl, "Z_%s_MnjI", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  dpd_contract244(&X1, &W, &Z1, 1, 2, 1, 1, 0);
  dpd_buf4_close(&W);
  sprintf(lbl, "Z_%s_MnIj", pert_x);
  dpd_buf4_sort(&Z1, CC_TMP0, pqsr, 0, 0, lbl);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qprs, 0, 0, lbl, 1);
  dpd_buf4_close(&Z1);
  dpd_file2_close(&X1);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  sprintf(lbl, "Z_%s_MnIj", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_contract444(&Z1, &Y2, &Z, 1, 1, 1, 1);
  dpd_buf4_close(&Y2);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&Z);

  sprintf(lbl, "Z_%s_AmIj (Am,Ij)", pert_y);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_y, 11, 0, 11, 0, 0, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  dpd_buf4_init(&Y2, CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  dpd_contract444(&W, &Y2, &Z1, 0, 0, 1, 0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Y2);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_contract424(&Z1, &X1, &Z, 1, 0, 0, -2, 1);
  dpd_file2_close(&X1);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
  polar += dpd_buf4_dot(&L2, &Z);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&Z);

  return polar;
}

}} // namespace psi::ccresponse
