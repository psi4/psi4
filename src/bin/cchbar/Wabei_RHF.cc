/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

void build_Z1(void);
void ZFW(dpdbuf4 *Z, dpdbuf4 *F, dpdbuf4 *W, double alpha, double beta);

/* Wabei_RHF(): Builds the Wabei HBAR matrix elements for CCSD for
** spin-adapted, closed-shell cases.  (Numbering of individual terms
** is given in Wabei.c.)  This version produces a final storage of
** (Ei,Ab), uses only <ia|bc>-type F integrals, and avoids producing
** any intermediates of v^3 storage, beyond the final target.  In
** addition, all memory bottlenecks larger than v^3.
**
** -TDC, 7/05
*/

void Wabei_RHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, C, F1, F2, W1, W2, Tau;
  double value;
  int Gef, Gei, Gab, Ge, Gi, Gf, Gmi, Gm, nrows, ncols, nlinks, EE, e, row, Gnm;
  int Gma, ma, m, a, Ga, Gb, I, i, mi, E, ei, ab, ba, b, BB, fb, bf, fe, ef, mb, am;
  double ***WW1, ***WW2;
  int h, incore, core_total, rowtot, coltot, maxrows;

  /** Term I **/
  /** <Ei|Ab> **/
  if(params.print & 2) {
    fprintf(outfile, "\tF -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_sort(&F, CC_HBAR, qpsr, 11, 5, "WAbEi (Ei,Ab)");
  dpd_buf4_close(&F);
  if(params.print & 2) fprintf(outfile, "done.\n");

  /** Term II **/
  /** - F_ME t_Mi^Ab **/
  if(params.print & 2) {
    fprintf(outfile, "\tFME*T2 -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");

  for(Gei=0; Gei < moinfo.nirreps; Gei++) {
    Gmi = Gab = Gei;  /* W and T2 are totally symmetric */

    dpd_buf4_mat_irrep_init(&T2, Gmi);
    dpd_buf4_mat_irrep_rd(&T2, Gmi);

    row = 0;
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gm = Ge;  /* Fme is totally symmetric */
      Gi = Gm ^ Gmi;

      W.matrix[Gei] = dpd_block_matrix(moinfo.occpi[Gi], W.params->coltot[Gei]);

      nrows = moinfo.occpi[Gm];
      ncols = moinfo.occpi[Gi] * W.params->coltot[Gei];

      if(nrows && ncols) {
	for(E=0; E < moinfo.virtpi[Ge]; E++) {
	  e = moinfo.vir_off[Ge] + E;

	  dpd_buf4_mat_irrep_rd_block(&W, Gei, W.row_offset[Gei][e], moinfo.occpi[Gi]);

	  C_DGEMV('t',nrows,ncols,-1.0,&T2.matrix[Gmi][row][0],ncols,&Fme.matrix[Gm][0][E],
		  moinfo.virtpi[Ge],1.0,W.matrix[Gei][0],1);

	  dpd_buf4_mat_irrep_wrt_block(&W, Gei, W.row_offset[Gei][e], moinfo.occpi[Gi]);
	}
      }

      row += moinfo.occpi[Gm] * moinfo.occpi[Gi];
      dpd_free_block(W.matrix[Gei], moinfo.occpi[Gi], W.params->coltot[Gei]);
    }

    dpd_buf4_mat_irrep_close(&T2, Gmi);
  }

  dpd_buf4_close(&W);
  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);
  dpd_buf4_close(&T2);

  if(params.print & 2) fprintf(outfile, "done.\n");

  /** Term IIIa **/
  /** <Ab|Ef> t_i^f **/
  if(params.print & 2) {
    fprintf(outfile, "\tB*T1 -> Wabei...");
    fflush(outfile);
  }
  /* Term re-written to use only (Ei,Ab) ordering on the target, TDC, 5/11/05 */
  /* Code modified to use only symmetric and antisymmetry <ab|cd> ints, TDC, 9/25/05 */
  dpd_buf4_init(&B, CC_BINTS, 0, 5, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1);
  dpd_file2_mat_rd(&T1);
  dpd_buf4_init(&Z1, CC_TMP0, 0, 11, 8, 11, 8, 0, "Z1(ei,a>=b)");
  dpd_buf4_scm(&Z1, 0.0); /* this scm is necessary for cases with empty occpi or virtpi irreps */
  for(Gef=0; Gef < moinfo.nirreps; Gef++) {
    Gei = Gab = Gef; /* W and B are totally symmetric */
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gf = Ge ^ Gef; Gi = Gf;  /* T1 is totally symmetric */
      B.matrix[Gef] = dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gef]);
      Z1.matrix[Gef] = dpd_block_matrix(moinfo.occpi[Gi],Z1.params->coltot[Gei]);
      nrows = moinfo.occpi[Gi]; ncols = Z1.params->coltot[Gef]; nlinks = moinfo.virtpi[Gf];
      if(nrows && ncols && nlinks) {
	for(E=0; E < moinfo.virtpi[Ge]; E++) {
	  e = moinfo.vir_off[Ge] + E;
	  dpd_buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
	  C_DGEMM('n','n',nrows,ncols,nlinks,0.5,T1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
		  0.0,Z1.matrix[Gei][0],ncols);
	  dpd_buf4_mat_irrep_wrt_block(&Z1, Gei, Z1.row_offset[Gei][e], moinfo.occpi[Gi]);
	}
      }
      dpd_free_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gef]);
      dpd_free_block(Z1.matrix[Gef], moinfo.occpi[Gi], Z1.params->coltot[Gei]);
    }
  }
  dpd_buf4_close(&Z1);
  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&B);

  dpd_buf4_init(&B, CC_BINTS, 0, 5, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1);
  dpd_file2_mat_rd(&T1);
  dpd_buf4_init(&Z2, CC_TMP0, 0, 11, 9, 11, 9, 0, "Z2(ei,a>=b)");
  dpd_buf4_scm(&Z2, 0.0); /* this scm is necessary for cases with empty occpi or virtpi irreps */
  for(Gef=0; Gef < moinfo.nirreps; Gef++) {
    Gei = Gab = Gef; /* W and B are totally symmetric */
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gf = Ge ^ Gef; Gi = Gf;  /* T1 is totally symmetric */
      B.matrix[Gef] = dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gef]);
      Z2.matrix[Gef] = dpd_block_matrix(moinfo.occpi[Gi],Z2.params->coltot[Gei]);
      nrows = moinfo.occpi[Gi]; ncols = Z2.params->coltot[Gef]; nlinks = moinfo.virtpi[Gf];
      if(nrows && ncols && nlinks) {
	for(E=0; E < moinfo.virtpi[Ge]; E++) {
	  e = moinfo.vir_off[Ge] + E;
	  dpd_buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
	  C_DGEMM('n','n',nrows,ncols,nlinks,0.5,T1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
		  0.0,Z2.matrix[Gei][0],ncols);
	  dpd_buf4_mat_irrep_wrt_block(&Z2, Gei, Z2.row_offset[Gei][e], moinfo.occpi[Gi]);
	}
      }
      dpd_free_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gef]);
      dpd_free_block(Z2.matrix[Gef], moinfo.occpi[Gi], Z2.params->coltot[Gei]);
    }
  }
  dpd_buf4_close(&Z2);
  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&B);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 11, 5, 11, 8, 0, "Z1(ei,a>=b)");
  dpd_buf4_init(&Z2, CC_TMP0, 0, 11, 5, 11, 9, 0, "Z2(ei,a>=b)");
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
  dpd_buf4_axpy(&Z1, &W, 1);
  dpd_buf4_axpy(&Z2, &W, 1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&Z1);

  if(params.print & 2) fprintf(outfile, "done.\n");

  /** Terms IIIc + IIId + IV **/
  if(params.print & 2) {
    fprintf(outfile, "\tD*T1*T2 + E*T2 -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&Z, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (nM,eI)");
  dpd_buf4_sort(&Z, CC_HBAR, rspq, 11, 0, "WMnIe (eI,nM)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
  dpd_buf4_init(&Z, CC_HBAR, 0, 11, 0, 11, 0, 0, "WMnIe (eI,nM)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  /*   dpd_contract444(&Z, &T, &W, 1, 1, 1, 1); */
  for(Gei=0; Gei < moinfo.nirreps; Gei++) {
    Gab = Gnm = Gei; /* Everything is totally symmetric here */
    nrows = T.params->rowtot[Gnm];
    ncols = T.params->coltot[Gab];
    if(nrows && ncols) {
      dpd_buf4_mat_irrep_init(&Z, Gei);
      dpd_buf4_mat_irrep_rd(&Z, Gei);
      dpd_buf4_mat_irrep_init(&T, Gnm);
      dpd_buf4_mat_irrep_rd(&T, Gnm);
      dpd_buf4_mat_irrep_row_init(&W, Gei);
      for(ei=0; ei < W.params->rowtot[Gei]; ei++) {
	dpd_buf4_mat_irrep_row_rd(&W, Gei, ei);
	C_DGEMV('t',nrows,ncols,1,T.matrix[Gei][0],ncols,Z.matrix[Gei][ei],1,
		1,W.matrix[Gei][0],1);
	dpd_buf4_mat_irrep_row_wrt(&W, Gei, ei);
      }
      dpd_buf4_mat_irrep_row_close(&W, Gei);
      dpd_buf4_mat_irrep_close(&T, Gnm);
      dpd_buf4_mat_irrep_close(&Z, Gei);
    }
  }
  dpd_buf4_close(&T);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&W);
  if(params.print & 2) fprintf(outfile, "done.\n");

  /*** Terms IIIB + V ***/
  /* Wabei <-- Z1(mi,fb) <ma|fe> - t(mi,fb) <ma|ef> - tau(mi,af) <mb|ef> */
  if(params.print & 2) {
    fprintf(outfile, "\t(T2+T1*T1) * F -> Wabei...");
    fflush(outfile);
  }
  build_Z1(); /* Z1(ib,mf) = 2 t(mi,fb) - t(mi,bf) - t(m,b) t(i,f) */

  /* The default algorithm for terms IIIb+V requires storage of two additional
  ** ovvv quantities on disk besides the <ia|bc> integrals and the Wabei target. */
  if(!params.wabei_lowdisk) {

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_sort(&F, CC_TMP0, prsq, 10, 5, "F <ia|bc> (ib,ca)");
    dpd_buf4_close(&F);

    /* can we run these contractions fully in core? */
    incore = 1;
    core_total = 0;
    for(h=0; h < moinfo.nirreps; h++) {
      coltot = F.params->coltot[h];
      if(coltot)
	maxrows = DPD_BIGNUM/coltot;
      if(maxrows < 1) {
	fprintf(stderr, "\nWabei_RHF Error: A single row of ovvv > 2 GW.\n");
	exit(PSI_RETURN_FAILURE);
      }
      else maxrows = DPD_BIGNUM;
      rowtot = F.params->rowtot[h];
      for(; rowtot > maxrows; rowtot -= maxrows) {
	if(core_total > (core_total + 2*maxrows*coltot)) incore = 0;
	else core_total += 2*maxrows*coltot;
      }
      if(core_total > (core_total + 2*rowtot*coltot)) incore = 0;
      core_total += 2*rowtot*coltot;
    }
    if(core_total > dpd_memfree()) incore = 0;
    if(!incore && (params.print & 2)) 
      fprintf(outfile, "\nUsing out-of-core algorithm for (T2+T1*T1)*F -> Wabei.\n");

    dpd_buf4_init(&W, CC_TMP0, 0, 10, 5, 10, 5, 0, "W(ib,ea)");
    dpd_buf4_init(&F, CC_TMP0, 0, 10, 5, 10, 5, 0, "F <ia|bc> (ib,ca)");
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ib,mf)");
    if(incore) dpd_contract444(&Z, &F, &W, 0, 1, 1, 0);
    else ZFW(&Z, &F, &W, 1, 0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&F);
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&Z, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    if(incore) dpd_contract444(&Z, &F, &W, 0, 1, -1, 1);
    else ZFW(&Z, &F, &W, -1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&F);
    dpd_buf4_sort_axpy(&W, CC_HBAR, rpsq, 11, 5, "WAbEi (Ei,Ab)", 1);
    dpd_buf4_close(&W);
    psio_close(CC_TMP0, 0); /* delete the extra ovvv quantities on disk */
    psio_open(CC_TMP0, PSIO_OPEN_NEW);
    dpd_buf4_init(&W, CC_TMP0, 0, 10, 5, 10, 5, 0, "W(ia,eb)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&Z, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tauIjAb (Ib,jA)");
    if(incore) dpd_contract444(&Z, &F, &W, 0, 1, -1, 0);
    else ZFW(&Z, &F, &W, -1, 0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&F);
    dpd_buf4_sort_axpy(&W, CC_HBAR, rpqs, 11, 5, "WAbEi (Ei,Ab)", 1);
    dpd_buf4_close(&W);
    psio_close(CC_TMP0, 0); /* delete the extra ovvv quantity on disk */
    psio_open(CC_TMP0, PSIO_OPEN_NEW);
  }
  /* This is an alternative algorithm for terms IIIb+V that stores no additional
  ** ovvv terms beyond <ia|bc> integrals and the WAbEi target, but it's very slow. */
  else {
    if(params.print & 2) 
      fprintf(outfile, "\nUsing low-disk algorithm for (T2+T1*T1)*F -> Wabei.\n");

    dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
    /* prepare rows of W in advance */
    for(Gei=0; Gei < moinfo.nirreps; Gei++)
      dpd_buf4_mat_irrep_row_init(&W, Gei);

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ib,mf)");
    dpd_buf4_sort(&Z, CC_TMP0, rpsq, 0, 5, "Z1(mi,fb)");
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z1(mi,fb)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    WW1 = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    WW2 = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    for(Gma=0; Gma < moinfo.nirreps; Gma++) {
      dpd_buf4_mat_irrep_row_init(&F, Gma);
      for(ma=0; ma < F.params->rowtot[Gma]; ma++) {
	m = F.params->roworb[Gma][ma][0];
	a = F.params->roworb[Gma][ma][1];
	Gm = F.params->psym[m];
	Ga = F.params->qsym[a];
	dpd_buf4_mat_irrep_row_rd(&F, Gma, ma);
	for(Gi=0; Gi < moinfo.nirreps; Gi++) {
	  Gmi = Gm ^ Gi;
	  dpd_buf4_mat_irrep_row_init(&Z, Gmi);
	  dpd_buf4_mat_irrep_row_init(&T2, Gmi);
	  dpd_buf4_mat_irrep_row_init(&Tau, Gmi);
	  /* allocate space for WW1 target */
	  for(Gf=0; Gf < moinfo.nirreps; Gf++) {
	    Gb = Gf ^ Gmi; /* Z is totally symmetric */
	    Ge = Gf ^ Gma; /* F is totally symmetric */
	    WW1[Gb] = dpd_block_matrix(moinfo.virtpi[Gb], moinfo.virtpi[Ge]);
	    WW2[Gb] = dpd_block_matrix(moinfo.virtpi[Gb], moinfo.virtpi[Ge]);
	  }
	  for(I=0; I < moinfo.occpi[Gi]; I++) {
	    i = moinfo.occ_off[Gi] + I;
	    mi = Z.params->rowidx[m][i];
	    dpd_buf4_mat_irrep_row_rd(&Z, Gmi, mi);
	    dpd_buf4_mat_irrep_row_rd(&T2, Gmi, mi);
	    dpd_buf4_mat_irrep_row_rd(&Tau, Gmi, mi);
	    for(Gf=0; Gf < moinfo.nirreps; Gf++) {
	      Gb = Gf ^ Gmi; /* Z is totally symmetric */
	      Ge = Gf ^ Gma; /* F is totally symmetric */
	      Gei = Ge ^ Gi;
	      nrows = moinfo.virtpi[Gb];
	      ncols = moinfo.virtpi[Ge];
	      nlinks = moinfo.virtpi[Gf];
	      fb = Z.col_offset[Gmi][Gf];
	      bf = Z.col_offset[Gmi][Gb];
	      fe = F.col_offset[Gma][Gf];
	      ef = F.col_offset[Gma][Ge];
	      if(nrows && ncols && nlinks) {
		C_DGEMM('t','n',nrows,ncols,nlinks,1,&(Z.matrix[Gmi][0][fb]),nrows,
			&(F.matrix[Gma][0][fe]),ncols,0,WW1[Gb][0],ncols);
		C_DGEMM('t','t',nrows,ncols,nlinks,-1,&(T2.matrix[Gmi][0][fb]),nrows,
			&(F.matrix[Gma][0][ef]),nlinks,1,WW1[Gb][0],ncols);
		C_DGEMM('n','t',nrows,ncols,nlinks,-1,&(Tau.matrix[Gmi][0][bf]),nlinks,
			&(F.matrix[Gma][0][ef]),nlinks,0,WW2[Gb][0],ncols);
	      }

	      for(E=0; E < moinfo.virtpi[Ge]; E++) {
		e = moinfo.vir_off[Ge] + E;
		ei = W.params->rowidx[e][i];
		dpd_buf4_mat_irrep_row_rd(&W, Gei, ei);
		for(BB=0; BB < moinfo.virtpi[Gb]; BB++) {
		  b = moinfo.vir_off[Gb] + BB;
		  ab = W.params->colidx[a][b];
		  ba = W.params->colidx[b][a];
		  W.matrix[Gei][0][ab] += WW1[Gb][BB][E];
		  W.matrix[Gei][0][ba] += WW2[Gb][BB][E];
		}
		dpd_buf4_mat_irrep_row_wrt(&W, Gei, ei);
	      }
	    } /* Gf */
	  } /* I */
	  dpd_buf4_mat_irrep_row_close(&Z, Gmi);
	  dpd_buf4_mat_irrep_row_close(&T2, Gmi);
	  dpd_buf4_mat_irrep_row_close(&Tau, Gmi);
	} /* Gi */
	  /* free W1 target */
	for(Gf=0; Gf < moinfo.nirreps; Gf++) {
	  Gb = Gf ^ Gmi; /* Z is totally symmetric */
	  Ge = Gf ^ Gma; /* F is totally symmetric */
	  dpd_free_block(WW1[Gb],moinfo.virtpi[Gb], moinfo.virtpi[Ge]);
	  dpd_free_block(WW2[Gb],moinfo.virtpi[Gb], moinfo.virtpi[Ge]);
	}
      } /* ma */
      dpd_buf4_mat_irrep_row_close(&F, Gma);
    } /* Gma */
    free(WW1);
    free(WW2);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&F);
    for(Gei=0; Gei < moinfo.nirreps; Gei++)
      dpd_buf4_mat_irrep_row_close(&W, Gei);
    dpd_buf4_close(&W);
  }

  if(params.print & 2) fprintf(outfile, "done.\n");

  /** Terms VI and VII **/

  /** t_in^bf  <Mn|Ef> + t_iN^bF <MN||EF> --> Z1_MEib **/
  if(params.print & 2) {
    fprintf(outfile, "\tT1*(C+D*T2) -> Wabei...");
    fflush(outfile);
  }
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ME,ib)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  dpd_contract444(&D, &T2, &Z, 0, 0, -1, 1);
  dpd_buf4_close(&D);
  dpd_buf4_close(&T2);
  dpd_buf4_sort(&Z, CC_TMP0, qrps, 11, 10, "Z(Ei,Mb)");
  dpd_buf4_close(&Z);

  /** - t_M^A ( <Ei|Mb> + Z(Ei,Mb) ) --> W(Ei,Ab) **/
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort_axpy(&D, CC_TMP0, spqr, 11, 10, "Z(Ei,Mb)", 1);
  dpd_buf4_close(&D);
  dpd_buf4_init(&Z, CC_TMP0, 0, 11, 10, 11, 10, 0, "Z(Ei,Mb)");
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  /*  dpd_contract244(&T1, &Z, &W, 0, 2, 1, -1, 1); */
  dpd_file2_mat_init(&T1); 
  dpd_file2_mat_rd(&T1);
  for(Gei=0; Gei < moinfo.nirreps; Gei++) {
    dpd_buf4_mat_irrep_init(&Z, Gei);
    dpd_buf4_mat_irrep_rd(&Z, Gei);
    dpd_buf4_mat_irrep_row_init(&W, Gei);
    for(ei=0; ei < Z.params->rowtot[Gei]; ei++) {
      dpd_buf4_mat_irrep_row_rd(&W, Gei, ei);
      for(Gm=0; Gm < moinfo.nirreps; Gm++) {
	Ga = Gm; /* T1 is totally symmetric */
	Gb = Gm ^ Gei; /* Z is totally symmetric */
	nrows = moinfo.virtpi[Ga];
	ncols = moinfo.virtpi[Gb];
	nlinks = moinfo.occpi[Gm];
	mb = Z.col_offset[Gei][Gm];
	ab = W.col_offset[Gei][Ga];
	if(nrows && ncols && nlinks)
	  C_DGEMM('t','n',nrows,ncols,nlinks,-1,T1.matrix[Gm][0],nrows,
		  &(Z.matrix[Gei][ei][mb]),ncols,1,&(W.matrix[Gei][0][ab]),ncols);
      }
      dpd_buf4_mat_irrep_row_wrt(&W, Gei, ei);
    }
    dpd_buf4_mat_irrep_close(&Z, Gei);
  }
  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z);

  /** - t_Ni^Af <mN|fE> --> Z2_mEiA **/
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(mE,iA)");
  dpd_contract444(&D, &T2, &Z, 0, 0, 1, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);
  dpd_buf4_sort(&Z, CC_TMP0, qrsp, 11, 11, "Z(Ei,Am)");
  dpd_buf4_close(&Z);

  /** t_m^b ( - <mA|iE> + Z(mA,iE) ) --> Z2(Ab,Ei) **/
  dpd_buf4_init(&C, CC_CINTS, 0, 11, 11, 11, 11, 0, "C <ai|bj>");
  dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(Ei,Am)");
  dpd_buf4_axpy(&C, &Z, -1.0);
  dpd_buf4_close(&C);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  /*  dpd_contract424(&Z, &T1, &W, 3, 0, 0, 1, 1); */
  dpd_file2_mat_init(&T1);
  dpd_file2_mat_rd(&T1);
  for(Gei=0; Gei < moinfo.nirreps; Gei++) {
    dpd_buf4_mat_irrep_init(&Z, Gei);
    dpd_buf4_mat_irrep_rd(&Z, Gei);
    dpd_buf4_mat_irrep_row_init(&W, Gei);
    for(ei=0; ei < Z.params->rowtot[Gei]; ei++) {
      dpd_buf4_mat_irrep_row_rd(&W, Gei, ei);
      for(Gm=0; Gm < moinfo.nirreps; Gm++) {
	Gb = Gm; /* T1 is totally symmetric */
	Ga = Gm ^ Gei; /* Z is totally symmetric */
	nrows = moinfo.virtpi[Ga];
	ncols = moinfo.virtpi[Gb];
	nlinks = moinfo.occpi[Gm];
	am = Z.col_offset[Gei][Ga];
	ab = W.col_offset[Gei][Ga];
	if(nrows && ncols && nlinks)
	  C_DGEMM('n','n',nrows,ncols,nlinks,1,&(Z.matrix[Gei][ei][am]),nlinks,
		  T1.matrix[Gm][0],ncols,1,&(W.matrix[Gei][0][ab]),ncols);
      }
      dpd_buf4_mat_irrep_row_wrt(&W, Gei, ei);
    }
    dpd_buf4_mat_irrep_close(&Z, Gei);
  }
  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);
  dpd_buf4_close(&W);
  dpd_buf4_close(&Z);

  if(params.print & 2) fprintf(outfile, "done.\n");
}

/*
** Generate intermediate needed for efficient evaluation of
** <am||ef> contributions to Wabei HBAR elements:
**
**  Z1(ib,mf) = 2 t(mi,fb) - t(im,fb) - t(i,f) * t(m,b) 
**
** TDC, 7/10/05
*/
void build_Z1(void)
{
  dpdbuf4 T2, Z1;
  dpdfile2 T1;
  int h, row, col, p, q, r, s, P, Q, R, S, psym, qsym, rsym, ssym;

  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  dpd_buf4_copy(&T2, CC_TMP0, "Z1(ib,mf)");
  dpd_buf4_close(&T2);

  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1);
  dpd_file2_mat_rd(&T1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ib,mf)");
  for(h=0; h < moinfo.nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Z1, h);
    dpd_buf4_mat_irrep_rd(&Z1, h);

    for(row=0; row < Z1.params->rowtot[h]; row++) {
      p = Z1.params->roworb[h][row][0];
      q = Z1.params->roworb[h][row][1];
      P = T1.params->rowidx[p];
      Q = T1.params->colidx[q];
      psym = T1.params->psym[p];
      qsym = T1.params->qsym[q];
      for(col=0; col < Z1.params->coltot[h]; col++) {
	r = Z1.params->colorb[h][col][0];
	s = Z1.params->colorb[h][col][1];
	R = T1.params->rowidx[r];
	S = T1.params->colidx[s];
	rsym = T1.params->psym[r];
	ssym = T1.params->qsym[s];

	if(qsym == rsym && psym == ssym)
	  Z1.matrix[h][row][col] -= T1.matrix[rsym][R][Q] * T1.matrix[psym][P][S];
      }
    }
    dpd_buf4_mat_irrep_wrt(&Z1, h);
    dpd_buf4_mat_irrep_close(&Z1, h);
  }
  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);

  dpd_buf4_close(&Z1);
}

/*
** ZFW(): Compute the product:
** 
**  W(ib,ea) = alpha * Z(ib,mf) * F(mf,ea) + beta * W(ib,ea)
**
**  This algorithm uses all the available core by first storing all of
**  Z(ib,mf) and splitting the remainder between F and W.  It then
**  makes a single pass through the W(ib,ea) file, but multiple passes
**  through the F(mf,ea) file.
**
**  Note that in this code, I think beta can be only 0 or 1.
**
** TDC, 8/05
**
*/
void ZFW(dpdbuf4 *Z, dpdbuf4 *F, dpdbuf4 *W, double alpha, double beta)
{
  int Gib, Gea, Gmf;
  int m, n;
  int rows_per_bucket, nbuckets, rows_left;
  int W_row_start, F_row_start;
  int nrows, ncols, nlinks;

  for(Gib=0; Gib < moinfo.nirreps; Gib++) {
    Gea = Gmf = Gib;  /* F and W are totally symmetric */
    dpd_buf4_mat_irrep_init(Z, Gib);
    dpd_buf4_mat_irrep_rd(Z, Gib);

    /* how many rows of F/W can we store in half the remaining core? */
    rows_per_bucket = dpd_memfree()/(2 * F->params->coltot[Gea]);
    if(rows_per_bucket > F->params->rowtot[Gib])
      rows_per_bucket = F->params->rowtot[Gib];
    nbuckets = (int) ceil((double) F->params->rowtot[Gib]/(double) rows_per_bucket);
    rows_left = F->params->rowtot[Gib] % rows_per_bucket;

    /* allocate space for the full buckets */
    dpd_buf4_mat_irrep_init_block(F, Gib, rows_per_bucket);
    dpd_buf4_mat_irrep_init_block(W, Gib, rows_per_bucket);

    ncols = W->params->coltot[Gea];

    for(m=0; m < (rows_left ? nbuckets-1 : nbuckets); m++) {
      W_row_start = m*rows_per_bucket;
      zero_arr(W->matrix[Gib][0], rows_per_bucket*ncols);

      if(beta == 1.0)
	dpd_buf4_mat_irrep_rd_block(W, Gib, W_row_start, rows_per_bucket);

      for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {
	F_row_start = n*rows_per_bucket;
	dpd_buf4_mat_irrep_rd_block(F, Gib, F_row_start, rows_per_bucket);

	nrows = rows_per_bucket;
	nlinks = rows_per_bucket;

	if(nrows && ncols && nlinks)
	  C_DGEMM('n', 'n', nrows, ncols, nlinks, alpha,
		  &(Z->matrix[Gib][W_row_start][F_row_start]), Z->params->coltot[Gmf],
		  F->matrix[Gmf][0], ncols, 1, W->matrix[Gib][0], ncols);
      }
      if(rows_left) {
	F_row_start = n*rows_per_bucket;
	dpd_buf4_mat_irrep_rd_block(F, Gib, F_row_start, rows_left);

	nrows = rows_per_bucket;
	nlinks = rows_left;

	if(nrows && ncols && nlinks)
	  C_DGEMM('n', 'n', nrows, ncols, nlinks, alpha,
		  &(Z->matrix[Gib][W_row_start][F_row_start]), Z->params->coltot[Gmf],
		  F->matrix[Gmf][0], ncols, 1, W->matrix[Gib][0], ncols);
      }

      dpd_buf4_mat_irrep_wrt_block(W, Gib, W_row_start, rows_per_bucket);
    }
    if(rows_left) {
      W_row_start = m*rows_per_bucket;
      zero_arr(W->matrix[Gib][0], rows_per_bucket*ncols);

      if(beta == 1.0)
	dpd_buf4_mat_irrep_rd_block(W, Gib, W_row_start, rows_left);

      for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {
	F_row_start = n*rows_per_bucket;
	dpd_buf4_mat_irrep_rd_block(F, Gib, F_row_start, rows_per_bucket);

	nrows = rows_left;
	nlinks = rows_per_bucket;

	if(nrows && ncols && nlinks)
	  C_DGEMM('n', 'n', nrows, ncols, nlinks, alpha,
		  &(Z->matrix[Gib][W_row_start][F_row_start]), Z->params->coltot[Gmf],
		  F->matrix[Gmf][0], ncols, 1, W->matrix[Gib][0], ncols);
      }
      if(rows_left) {
	F_row_start = n*rows_per_bucket;
	dpd_buf4_mat_irrep_rd_block(F, Gib, F_row_start, rows_left);

	nrows = rows_left;
	nlinks = rows_left;

	if(nrows && ncols && nlinks)
	  C_DGEMM('n', 'n', nrows, ncols, nlinks, alpha,
		  &(Z->matrix[Gib][W_row_start][F_row_start]), Z->params->coltot[Gmf],
		  F->matrix[Gmf][0], ncols, 1, W->matrix[Gib][0], ncols);
      }
      dpd_buf4_mat_irrep_wrt_block(W, Gib, W_row_start, rows_left);
    } /* m */

    dpd_buf4_mat_irrep_close_block(F, Gib, rows_per_bucket);
    dpd_buf4_mat_irrep_close_block(W, Gib, rows_per_bucket);

    dpd_buf4_mat_irrep_close(Z, Gib);

  } /* Gib */
}

}} // namespace psi::cchbar
