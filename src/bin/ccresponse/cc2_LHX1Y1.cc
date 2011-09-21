/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double cc2_LHX1Y1(const char *pert_x, int irrep_x, double omega_x,
		  const char *pert_y, int irrep_y, double omega_y)
{
	
  int am, a, A, m, M, fe, ef, f, e, E;
  int GW, GZae, Ga, Gm, Gf, Gam, Gef, Gab, Gei, GX;
  int hxbuf, hzbuf, Gi, Gj, Ge, GZ;
  int ab, mb, colx, colz, rowz, rowx;
  int ncols, nrows, nlinks;
  double *Xt;
  dpdfile2 F, X1, Y1, Zmi, Zae, ZIA, L1, t1;
  dpdbuf4 Z1, Z2, I, W1, ZIjAb, L2, Z, X, W, B;
  double polar;
  char lbl[32];
	
  /* The Lambda 1 contractions */
  dpd_file2_init(&ZIA, CC_TMP0, 0, 0, 1, "ZIA");
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
  dpd_file2_init(&Y1, CC_OEI, irrep_y, 0, 1, lbl);
	
  /* Contraction of FME, XIE, YMA */
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  sprintf(lbl, "Z_%s_MI" , pert_x);
  dpd_file2_init(&Zmi, CC_TMP0, irrep_x, 0, 0, lbl);
  dpd_contract222(&F, &X1, &Zmi, 0, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_contract222(&Zmi, &Y1, &ZIA, 1, 1, -1, 0);
	
  /* Contraction of FME, XMA, YIE */
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_init(&Zmi, CC_TMP0, irrep_x, 0, 0, lbl);
  dpd_contract222(&F, &Y1, &Zmi, 0, 0, 1, 0);
  dpd_file2_close(&F);
  dpd_contract222(&Zmi, &X1, &ZIA, 1, 1, -1, 1);
  dpd_file2_close(&Zmi);
	
  /* Contraction of WAMEF, XIE, YMF */
  /** Begin out-of-core dot24 contraction **/
  sprintf(lbl, "Z_%s_AE" , pert_x);
  dpd_file2_init(&Zae, CC_TMP0, irrep_y, 1, 1, lbl);
  dpd_file2_scm(&Zae, 0);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");

  GW = W.file.my_irrep; 
  GZae = Zae.my_irrep;

  dpd_file2_mat_init(&Zae);
  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);

  for(Gam=0; Gam < moinfo.nirreps; Gam++) {
    dpd_buf4_mat_irrep_row_init(&W, Gam);
    Xt = init_array(W.params->coltot[Gam]);

    for(am=0; am < W.params->rowtot[Gam]; am++) {
      dpd_buf4_mat_irrep_row_rd(&W, Gam, am);
      a = W.params->roworb[Gam][am][0];
      m = W.params->roworb[Gam][am][1];
      Ga = W.params->psym[a];
      Gm = W.params->qsym[m];
      Ge = Ga^GZae;  /* Zae is not totally symmetric */
      Gf = Gam^GW^Ge; /* X1 is not totally symmetric */
      Gef = Gam^GW;
      A = a - W.params->poff[Ga];
      M = m - W.params->qoff[Gm];

      zero_arr(Xt, W.params->coltot[Gam]);

      /* build spin-adapted W-integrals for current am */
      for(ef=0; ef < W.params->coltot[Gef]; ef++) {
	e = W.params->colorb[Gef][ef][0];
	f = W.params->colorb[Gef][ef][1];
	fe = W.params->colidx[f][e];
	Xt[ef] = 2.0 * W.matrix[Gam][0][ef] - W.matrix[Gam][0][fe];
      }
	
      nrows = moinfo.virtpi[Ge];
      ncols = moinfo.virtpi[Gf];
      if(nrows && ncols)
	C_DGEMV('n',nrows,ncols,1.0,&Xt[W.col_offset[Gam][Ge]],ncols,
		X1.matrix[Gm][M],1,1.0,
		Zae.matrix[Ga][A],1);
    }

    free(Xt);
    dpd_buf4_mat_irrep_row_close(&W, Gam);
  }
  dpd_buf4_close(&W);

  dpd_file2_mat_close(&X1);
  dpd_file2_mat_wrt(&Zae);
  dpd_file2_mat_close(&Zae);
  /** End out-of-core dot24 contraction **/

  dpd_contract222(&Y1, &Zae, &ZIA, 0, 0, 1, 1);
  dpd_file2_close(&Zae);
	
  /* Contraction of WAMEF, XMF, YIE */
  /** Begin out-of-core dot24 contraction **/
  sprintf(lbl, "Z_%s_AE" , pert_y);
  dpd_file2_init(&Zae, CC_TMP0, irrep_y, 1, 1, lbl);
  dpd_file2_scm(&Zae, 0);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");

  GW = W.file.my_irrep; 
  GZae = Zae.my_irrep;

  dpd_file2_mat_init(&Zae);
  dpd_file2_mat_init(&Y1);
  dpd_file2_mat_rd(&Y1);

  for(Gam=0; Gam < moinfo.nirreps; Gam++) {
    dpd_buf4_mat_irrep_row_init(&W, Gam);
    Xt = init_array(W.params->coltot[Gam]);

    for(am=0; am < W.params->rowtot[Gam]; am++) {
      dpd_buf4_mat_irrep_row_rd(&W, Gam, am);
      a = W.params->roworb[Gam][am][0];
      m = W.params->roworb[Gam][am][1];
      Ga = W.params->psym[a];
      Gm = W.params->qsym[m];
      Ge = Ga^GZae;  /* Zae is not totally symmetric */
      Gf = Gam^GW^Ge; /* Y1 is not totally symmetric */
      Gef = Gam^GW;
      A = a - W.params->poff[Ga];
      M = m - W.params->qoff[Gm];

      zero_arr(Xt, W.params->coltot[Gam]);

      /* build spin-adapted W-integrals for current am */
      for(ef=0; ef < W.params->coltot[Gef]; ef++) {
	e = W.params->colorb[Gef][ef][0];
	f = W.params->colorb[Gef][ef][1];
	fe = W.params->colidx[f][e];
	Xt[ef] = 2.0 * W.matrix[Gam][0][ef] - W.matrix[Gam][0][fe];
      }
	
      nrows = moinfo.virtpi[Ge];
      ncols = moinfo.virtpi[Gf];
      if(nrows && ncols)
	C_DGEMV('n',nrows,ncols,1.0,&Xt[W.col_offset[Gam][Ge]],ncols,
		Y1.matrix[Gm][M],1,1.0,
		Zae.matrix[Ga][A],1);
    }

    free(Xt);
    dpd_buf4_mat_irrep_row_close(&W, Gam);
  }
  dpd_buf4_close(&W);
  dpd_file2_mat_close(&Y1);

  dpd_file2_mat_wrt(&Zae);
  dpd_file2_mat_close(&Zae);
  /** End out-of-core dot24 contraction **/

  dpd_contract222(&X1, &Zae, &ZIA, 0, 0, 1, 1);
  dpd_file2_close(&Zae);
	
  /* Contraction of WAMEF, XMA, YNE */
  sprintf(lbl, "Z_%s_MI" , pert_y);
  dpd_file2_init(&Zmi, CC_TMP0, irrep_y, 0, 0, lbl);
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
  dpd_dot13(&Y1, &W1, &Zmi, 0, 0, 1, 0);
  dpd_buf4_close(&W1);
  dpd_contract222(&Zmi, &X1, &ZIA, 1, 1, 1, 1);
  dpd_file2_close(&Zmi);
	
  /* Contraction of WAMEF, XMA, YNE */
  sprintf(lbl, "Z_%s_MI" , pert_x);
  dpd_file2_init(&Zmi, CC_TMP0, irrep_x, 0, 0, lbl);
  dpd_buf4_init(&W1, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
  dpd_dot13(&X1, &W1, &Zmi, 0, 0, 1, 0);
  dpd_buf4_close(&W1);
  dpd_contract222(&Zmi, &Y1, &ZIA, 1, 1, 1, 1);
  dpd_file2_close(&Zmi);
	
  dpd_file2_close(&Y1);
  dpd_file2_close(&X1);
	
  /* Final contraction of ZIA intermediate with LIA */
  dpd_file2_init(&L1, CC_LAMPS, 0, 0, 1, "LIA 0 -1");
  polar = 2.0 * dpd_file2_dot(&ZIA, &L1);
  dpd_file2_close(&L1);
  dpd_file2_close(&ZIA);
	
  /*   fprintf(outfile, "L(1)HX1Y1 = %20.12f\n", polar); */
	
  /* The Lambda 2 contractions */
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
  dpd_file2_init(&Y1, CC_OEI, irrep_y, 0, 1, lbl);
	
	
  /* Contraction with Wmnij */
  sprintf(lbl, "Z_%s_MbIj", pert_y);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_y, 10, 0, 10, 0, 0, lbl);
  dpd_buf4_init(&W1, CC2_HET1, 0, 0, 0, 0, 0, 0, "CC2 WMnIj");
  dpd_contract424(&W1, &Y1, &Z1, 1, 0, 1, 1, 0);
  dpd_buf4_close(&W1);
	
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "ZIjAb");
  dpd_contract244(&X1, &Z1, &Z, 0, 0, 1, 1, 0);
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_scm(&ZIjAb, 0);
  dpd_buf4_axpy(&Z, &ZIjAb, 1);
  dpd_buf4_close(&ZIjAb);
  dpd_buf4_sort_axpy(&Z, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z);
	
  /* B -> Wabef */

  dpd_file2_mat_init(&Y1);
  dpd_file2_mat_rd(&Y1);

  /* Plus Combination */
  dpd_buf4_init(&B, CC_BINTS, 0, 5, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
  sprintf(lbl, "Z1_%s_(ei,a>=b)", pert_y);
  dpd_buf4_init(&Z1, CC_TMP8, irrep_y, 11, 8, 11, 8, 0, lbl);
  dpd_buf4_scm(&Z1, 0);

  for(Gef=0; Gef < moinfo.nirreps; Gef++) {
    Gab = Gef; /* B is totally symmetric */
    Gei = Gab ^ irrep_y; /* Z is not totally symmetrix */
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gf = Ge ^ Gef;
      Gi = Gf ^ irrep_y;  /* Y1 is not totally symmetric */
      B.matrix[Gef] = dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gef]);
      Z1.matrix[Gei] = dpd_block_matrix(moinfo.occpi[Gi],Z1.params->coltot[Gab]);
      nrows = moinfo.occpi[Gi];
      ncols = Z1.params->coltot[Gab];
      nlinks = moinfo.virtpi[Gf];
      if(nrows && ncols && nlinks) {
	for(E=0; E < moinfo.virtpi[Ge]; E++) {
	  e = moinfo.vir_off[Ge] + E;
	  dpd_buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
	  C_DGEMM('n','n',nrows,ncols,nlinks,0.5,Y1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
		  0.0,Z1.matrix[Gei][0],ncols);
	  dpd_buf4_mat_irrep_wrt_block(&Z1, Gei, Z1.row_offset[Gei][e], moinfo.occpi[Gi]);
	}
      }
      dpd_free_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gef]);
      dpd_free_block(Z1.matrix[Gei], moinfo.occpi[Gi], Z1.params->coltot[Gab]);
    }
  }

  dpd_buf4_close(&Z1);
  dpd_buf4_close(&B);

  /* Minus Combination */
  dpd_buf4_init(&B, CC_BINTS, 0, 5, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
  sprintf(lbl, "Z2_%s_(ei,a>=b)", pert_y);
  dpd_buf4_init(&Z2, CC_TMP8, irrep_y, 11, 9, 11, 9, 0, lbl);
  dpd_buf4_scm(&Z2, 0);

  for(Gef=0; Gef < moinfo.nirreps; Gef++) {
    Gab = Gef; /* B is totally symmetric */
    Gei = Gab ^ irrep_y; /* Z is not totally symmetrix */
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gf = Ge ^ Gef;
      Gi = Gf ^ irrep_y;  /* Y1 is not totally symmetric */
      B.matrix[Gef] = dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gef]);
      Z2.matrix[Gei] = dpd_block_matrix(moinfo.occpi[Gi],Z2.params->coltot[Gab]);
      nrows = moinfo.occpi[Gi];
      ncols = Z2.params->coltot[Gab];
      nlinks = moinfo.virtpi[Gf];
      if(nrows && ncols && nlinks) {
	for(E=0; E < moinfo.virtpi[Ge]; E++) {
	  e = moinfo.vir_off[Ge] + E;
	  dpd_buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
	  C_DGEMM('n','n',nrows,ncols,nlinks,0.5,Y1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
		  0.0,Z2.matrix[Gei][0],ncols);
	  dpd_buf4_mat_irrep_wrt_block(&Z2, Gei, Z2.row_offset[Gei][e], moinfo.occpi[Gi]);
	}
      }
      dpd_free_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gef]);
      dpd_free_block(Z2.matrix[Gei], moinfo.occpi[Gi], Z2.params->coltot[Gab]);
    }
  }

  dpd_buf4_close(&Z2);
  dpd_buf4_close(&B);
  dpd_file2_mat_close(&Y1);

  sprintf(lbl, "Z_%s_AbEj (Ej,Ab)", pert_y);
  dpd_buf4_init(&Z, CC_TMP8, irrep_y, 11, 5, 11, 5, 0, lbl);
  dpd_buf4_scm(&Z, 0);
  sprintf(lbl, "Z1_%s_(ei,a>=b)", pert_y);
  dpd_buf4_init(&Z1, CC_TMP8, irrep_y, 11, 5, 11, 8, 0, lbl);
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z1);
  sprintf(lbl, "Z2_%s_(ei,a>=b)", pert_y);
  dpd_buf4_init(&Z2, CC_TMP8, irrep_y, 11, 5, 11, 9, 0, lbl);
  dpd_buf4_axpy(&Z2, &Z, 1);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&Z);

  sprintf(lbl, "Z_%s_AbEj", pert_y);
  dpd_buf4_init(&Z, CC_TMP9, irrep_y, 5, 11, 5, 11, 0, lbl);
  dpd_buf4_scm(&Z, 0);
  dpd_buf4_close(&Z);

  sprintf(lbl, "Z_%s_AbEj (Ej,Ab)", pert_y);
  dpd_buf4_init(&Z, CC_TMP8, irrep_y, 11, 5, 11, 5, 0, lbl);
  sprintf(lbl, "Z_%s_AbEj", pert_y);
  dpd_buf4_sort_axpy(&Z, CC_TMP9, rspq, 5, 11, lbl, 1);
  dpd_buf4_close(&Z);

  /* Free some disk space */
  psio_close(CC_TMP8, 0);
  psio_open(CC_TMP8, 0);

  /*   sprintf(lbl, "Z_%s_AbEj", pert_y); */
  /*   dpd_buf4_init(&Z, CC_TMP9, irrep_y, 5, 11, 5, 11, 0, lbl); */
  /*   dpd_buf4_init(&I, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>"); */
  /*   dpd_contract424(&I, &Y1, &Z, 3, 1, 0, 1, 0); */
  /*   dpd_buf4_close(&I); */
  /*   dpd_buf4_close(&Z); */

  /** Begin out-of-core contract244 **/
  dpd_buf4_init(&Z, CC_TMP0, 0, 5, 0, 5, 0, 0, "ZIjAb (Ab,Ij)");
  sprintf(lbl, "Z_%s_AbEj", pert_y);
  dpd_buf4_init(&X, CC_TMP9, irrep_y, 5, 11, 5, 11, 0, lbl);

  /* Symmetry Info */
  GX = X.file.my_irrep; 
  GZ = Z.file.my_irrep;

  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);

  for(hxbuf=0; hxbuf < moinfo.nirreps; hxbuf++) {
    hzbuf = hxbuf;

    dpd_buf4_mat_irrep_row_init(&X, hxbuf);
    dpd_buf4_mat_irrep_row_init(&Z, hzbuf);

    /* Loop over rows of the X factor and the target */
    for(ab=0; ab < Z.params->rowtot[hzbuf]; ab++) {

      dpd_buf4_mat_irrep_row_zero(&X, hxbuf, ab);
      dpd_buf4_mat_irrep_row_rd(&X, hxbuf, ab);

      dpd_buf4_mat_irrep_row_zero(&Z, hzbuf, ab);

      for(Gj=0; Gj < moinfo.nirreps; Gj++) {
	Ge = Gj^hxbuf^GX;
	Gi = Gj^hzbuf^GZ;

	rowx = X.params->rpi[Ge];
	colx = X.params->spi[Gj];
	rowz = Z.params->rpi[Gi];
	colz = Z.params->spi[Gj];

	if(rowz && colz && rowx && colx) {
	  C_DGEMM('n','n',rowz,colz,rowx,1.0,
		  &(X1.matrix[Gi][0][0]),rowx,
		  &(X.matrix[hxbuf][0][X.col_offset[hxbuf][Ge]]),colx,0.0,
		  &(Z.matrix[hzbuf][0][Z.col_offset[hzbuf][Gi]]),colz);
	}
      }

      dpd_buf4_mat_irrep_row_wrt(&Z, hzbuf, ab);
    }

    dpd_buf4_mat_irrep_row_close(&X, hxbuf);
    dpd_buf4_mat_irrep_row_close(&Z, hzbuf);
  }
  dpd_file2_mat_close(&X1);
  dpd_buf4_close(&X);
  dpd_buf4_sort_axpy(&Z, CC_TMP0, rspq, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_sort_axpy(&Z, CC_TMP0, srqp, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z);
  /** End out-of-core contract244 **/
  psio_close(CC_TMP9, 0);
  psio_open(CC_TMP9, 0);
	
  /* D -> Wabef */
  sprintf(lbl, "Z_%s_MnIf", pert_x);
  dpd_buf4_init(&Z, CC_TMP0, irrep_x, 0, 10, 0, 10, 0, lbl);
  dpd_buf4_init(&I, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_contract244(&X1, &I, &Z, 1, 2, 1, 1, 0);
  dpd_buf4_close(&I);
	
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 0, 0, 0, 0, "ZMnIj temp");
  dpd_contract424(&Z, &Y1, &Z1, 3, 1, 0, 1, 0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 0, 0, 0, 0, "ZMnIj");
  dpd_buf4_scm(&Z, 0);
  dpd_buf4_axpy(&Z1, &Z, 1);
  dpd_buf4_close(&Z);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qpsr, 0, 0, "ZMnIj", 1);
  dpd_buf4_close(&Z1);
	
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "ZMbIj");
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 0, 0, 0, 0, "ZMnIj");
  dpd_contract424(&Z1, &t1, &Z, 1, 0, 1, 1, 0);
  dpd_buf4_close(&Z1);
	
  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "ZIjAb");
  dpd_contract244(&t1, &Z, &Z1, 0, 0, 1, 0.5, 0);
  dpd_buf4_close(&Z);
  dpd_file2_close(&t1);
  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_axpy(&Z1, &ZIjAb, 1);
  dpd_buf4_close(&ZIjAb);
  dpd_buf4_sort_axpy(&Z1, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z1);
	
  /* F -> Wabef */
  /** Begin out-of-core contract244 **/
  sprintf(lbl, "Z_%s_MbEj (Mb,jE)", pert_x);
  dpd_buf4_init(&Z, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  dpd_buf4_init(&X, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");

  /* Symmetry Info */
  GX = X.file.my_irrep; 
  GZ = Z.file.my_irrep;

  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);

  for(hxbuf=0; hxbuf < moinfo.nirreps; hxbuf++) {
    hzbuf = hxbuf;

    dpd_buf4_mat_irrep_row_init(&X, hxbuf);
    dpd_buf4_mat_irrep_row_init(&Z, hzbuf);

    /* Loop over rows of the X factor and the target */
    for(mb=0; mb < Z.params->rowtot[hzbuf]; mb++) {

      dpd_buf4_mat_irrep_row_zero(&X, hxbuf, mb);
      dpd_buf4_mat_irrep_row_rd(&X, hxbuf, mb);

      dpd_buf4_mat_irrep_row_zero(&Z, hzbuf, mb);

      for(Gj=0; Gj < moinfo.nirreps; Gj++) {
	Ge = Gj^hxbuf^GX;
	Gi = Gj^hzbuf^GZ;

	rowx = X.params->rpi[Ge];
	colx = X.params->spi[Gj];
	rowz = Z.params->rpi[Gi];
	colz = Z.params->spi[Gj];

	if(rowz && colz && rowx && colx) {
	  C_DGEMM('n','n',rowz,colz,rowx,1.0,
		  &(X1.matrix[Gi][0][0]),rowx,
		  &(X.matrix[hxbuf][0][X.col_offset[hxbuf][Ge]]),colx,0.0,
		  &(Z.matrix[hzbuf][0][Z.col_offset[hzbuf][Gi]]),colz);
	}
      }

      dpd_buf4_mat_irrep_row_wrt(&Z, hzbuf, mb);
    }

    dpd_buf4_mat_irrep_row_close(&X, hxbuf);
    dpd_buf4_mat_irrep_row_close(&Z, hzbuf);
  }
  dpd_file2_mat_close(&X1);
  dpd_buf4_close(&X);
  dpd_buf4_close(&Z);
  /** End out-of-core contract244 **/

  /*   sprintf(lbl, "Z_%s_MbEj (jE,Mb)", pert_x); */
  /*   dpd_buf4_init(&Z, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl); */
  /*   dpd_buf4_init(&I, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>"); */
  /*   dpd_contract244(&X1, &I, &Z, 1, 2, 0, 1, 0); */
  /*   dpd_buf4_close(&I); */
  /*   sprintf(lbl, "Z_%s_MbEj (Mb,jE)", pert_x); */
  /*   dpd_buf4_sort(&Z, CC_TMP0, rspq, 10, 10, lbl); */
  /*   dpd_buf4_close(&Z); */

  sprintf(lbl, "Z_%s_MbeJ (Mb,eJ)", pert_x);
  dpd_buf4_init(&Z, CC_TMP0, irrep_x, 10, 11, 10, 11, 0, lbl);
  dpd_buf4_init(&I, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_contract424(&I, &X1, &Z, 3, 1, 0, 1, 0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&Z);
	
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "ZMbIj");
  sprintf(lbl, "Z_%s_MbEj (Mb,jE)", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 10, 10, 10, 10, 0, lbl);
  dpd_contract424(&Z1, &Y1, &Z, 3, 1, 0, 1, 0);
  dpd_buf4_close(&Z1);
  sprintf(lbl, "Z_%s_MbeJ (Mb,eJ)", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 10, 11, 10, 11, 0, lbl);
  dpd_contract244(&Y1, &Z1, &Z, 1, 2, 1, 1, 1);
  dpd_buf4_close(&Z1);
  dpd_buf4_close(&Z);
	
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "ZIjAb");
  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 0, 10, 0, 0, "ZMbIj");
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&t1, &Z1, &Z, 0, 0, 1, -1, 0);
  dpd_file2_close(&t1);
  dpd_buf4_close(&Z1);
  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_axpy(&Z, &ZIjAb, 1);
  dpd_buf4_close(&ZIjAb);
  dpd_buf4_sort_axpy(&Z, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z);
	
  /* Contraction with Wmbej */	
  sprintf(lbl, "X_%s_jbMI", pert_x);
  dpd_buf4_init(&Z, CC_TMP0, irrep_x, 10, 0, 10, 0, 0, lbl);
  dpd_buf4_init(&W1, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj (ME,jb)");
  dpd_contract424(&W1, &X1, &Z, 1, 1, 0, -1, 0);
  dpd_buf4_close(&W1);
  sprintf(lbl, "X_%s_IjMb", pert_x);
  dpd_buf4_sort(&Z, CC_TMP0, sprq, 0, 10, lbl);
  dpd_buf4_close(&Z);
	
  sprintf(lbl, "X_%s_IbMj", pert_x);
  dpd_buf4_init(&Z, CC_TMP0, irrep_x, 10, 0, 10, 0, 0, lbl);
  dpd_buf4_init(&W1, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Me,Jb)");
  dpd_contract424(&W1, &X1, &Z, 1, 1, 0, 1, 0);
  dpd_buf4_close(&W1);
  sprintf(lbl, "X_%s_IjMb", pert_x);
  dpd_buf4_sort_axpy(&Z, CC_TMP0, psrq, 0, 10, lbl, 1);
  dpd_buf4_close(&Z);
	
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "ZIjAb");
  sprintf(lbl, "X_%s_IjMb", pert_x);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_x, 0, 10, 0, 10, 0, lbl);
  dpd_contract244(&Y1, &Z1, &Z, 0, 2, 1, 1, 0);
  dpd_buf4_close(&Z1);
		
  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_axpy(&Z, &ZIjAb, 1);
  dpd_buf4_close(&ZIjAb);
  dpd_buf4_sort_axpy(&Z, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z);
	
  sprintf(lbl, "Y_%s_jbMI", pert_y);
  dpd_buf4_init(&Z, CC_TMP0, irrep_y, 10, 0, 10, 0, 0, lbl);
  dpd_buf4_init(&W1, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbEj (ME,jb)");
  dpd_contract424(&W1, &Y1, &Z, 1, 1, 0, -1, 0);
  dpd_buf4_close(&W1);
  sprintf(lbl, "Y_%s_IjMb", pert_y);
  dpd_buf4_sort(&Z, CC_TMP0, sprq, 0, 10, lbl);
  dpd_buf4_close(&Z);
	
  sprintf(lbl, "Y_%s_IbMj", pert_y);
  dpd_buf4_init(&Z, CC_TMP0, irrep_y, 10, 0, 10, 0, 0, lbl);
  dpd_buf4_init(&W1, CC2_HET1, 0, 10, 10, 10, 10, 0, "CC2 WMbeJ (Me,Jb)");
  dpd_contract424(&W1, &Y1, &Z, 1, 1, 0, 1, 0);
  dpd_buf4_close(&W1);
  sprintf(lbl, "Y_%s_IjMb", pert_y);
  dpd_buf4_sort_axpy(&Z, CC_TMP0, psrq, 0, 10, lbl, 1);
  dpd_buf4_close(&Z);
	
  dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "ZIjAb");
  sprintf(lbl, "Y_%s_IjMb", pert_y);
  dpd_buf4_init(&Z1, CC_TMP0, irrep_y, 0, 10, 0, 10, 0, lbl);
  dpd_contract244(&X1, &Z1, &Z, 0, 2, 1, 1, 0);
  dpd_buf4_close(&Z1);
		
  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  dpd_buf4_axpy(&Z, &ZIjAb, 1);
  dpd_buf4_close(&ZIjAb);
  dpd_buf4_sort_axpy(&Z, CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) Final", 1);
  dpd_buf4_close(&Z);
	
  /* Close the X and Y matices */
  dpd_file2_close(&Y1);
  dpd_file2_close(&X1);
	 
  /* Final contraction with LIJAB */
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
  dpd_buf4_init(&ZIjAb, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  polar += dpd_buf4_dot(&L2, &ZIjAb);
  dpd_buf4_close(&ZIjAb);
  dpd_buf4_close(&L2);
	 
  return polar;
}

}} // namespace psi::ccresponse
