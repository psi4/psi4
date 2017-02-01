/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
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
    outfile->Printf( "\tF -> Wabei...");

  }
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_HBAR, qpsr, 11, 5, "WAbEi (Ei,Ab)");
  global_dpd_->buf4_close(&F);
  if(params.print & 2) outfile->Printf( "done.\n");

  /** Term II **/
  /** - F_ME t_Mi^Ab **/
  if(params.print & 2) {
    outfile->Printf( "\tFME*T2 -> Wabei...");

  }
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->file2_mat_init(&Fme);
  global_dpd_->file2_mat_rd(&Fme);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");

  for(Gei=0; Gei < moinfo.nirreps; Gei++) {
    Gmi = Gab = Gei;  /* W and T2 are totally symmetric */

    global_dpd_->buf4_mat_irrep_init(&T2, Gmi);
    global_dpd_->buf4_mat_irrep_rd(&T2, Gmi);

    row = 0;
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gm = Ge;  /* Fme is totally symmetric */
      Gi = Gm ^ Gmi;

      W.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gi], W.params->coltot[Gei]);

      nrows = moinfo.occpi[Gm];
      ncols = moinfo.occpi[Gi] * W.params->coltot[Gei];

      if(nrows && ncols) {
	for(E=0; E < moinfo.virtpi[Ge]; E++) {
	  e = moinfo.vir_off[Ge] + E;

	  global_dpd_->buf4_mat_irrep_rd_block(&W, Gei, W.row_offset[Gei][e], moinfo.occpi[Gi]);

	  C_DGEMV('t',nrows,ncols,-1.0,&T2.matrix[Gmi][row][0],ncols,&Fme.matrix[Gm][0][E],
		  moinfo.virtpi[Ge],1.0,W.matrix[Gei][0],1);

	  global_dpd_->buf4_mat_irrep_wrt_block(&W, Gei, W.row_offset[Gei][e], moinfo.occpi[Gi]);
	}
      }

      row += moinfo.occpi[Gm] * moinfo.occpi[Gi];
      global_dpd_->free_dpd_block(W.matrix[Gei], moinfo.occpi[Gi], W.params->coltot[Gei]);
    }

    global_dpd_->buf4_mat_irrep_close(&T2, Gmi);
  }

  global_dpd_->buf4_close(&W);
  global_dpd_->file2_mat_close(&Fme);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  if(params.print & 2) outfile->Printf( "done.\n");

  /** Term IIIa **/
  /** <Ab|Ef> t_i^f **/
  if(params.print & 2) {
    outfile->Printf( "\tB*T1 -> Wabei...");

  }
  /* Term re-written to use only (Ei,Ab) ordering on the target, TDC, 5/11/05 */
  /* Code modified to use only symmetric and antisymmetry <ab|cd> ints, TDC, 9/25/05 */
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 8, 11, 8, 0, "Z1(ei,a>=b)");
  global_dpd_->buf4_scm(&Z1, 0.0); /* this scm is necessary for cases with empty occpi or virtpi irreps */
  for(Gef=0; Gef < moinfo.nirreps; Gef++) {
    Gei = Gab = Gef; /* W and B are totally symmetric */
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gf = Ge ^ Gef; Gi = Gf;  /* T1 is totally symmetric */
      B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gef]);
      Z1.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gi],Z1.params->coltot[Gei]);
      nrows = moinfo.occpi[Gi]; ncols = Z1.params->coltot[Gef]; nlinks = moinfo.virtpi[Gf];
      if(nrows && ncols && nlinks) {
	for(E=0; E < moinfo.virtpi[Ge]; E++) {
	  e = moinfo.vir_off[Ge] + E;
	  global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
	  C_DGEMM('n','n',nrows,ncols,nlinks,0.5,T1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
		  0.0,Z1.matrix[Gei][0],ncols);
	  global_dpd_->buf4_mat_irrep_wrt_block(&Z1, Gei, Z1.row_offset[Gei][e], moinfo.occpi[Gi]);
	}
      }
      global_dpd_->free_dpd_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gef]);
      global_dpd_->free_dpd_block(Z1.matrix[Gef], moinfo.occpi[Gi], Z1.params->coltot[Gei]);
    }
  }
  global_dpd_->buf4_close(&Z1);
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);

  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 9, 11, 9, 0, "Z2(ei,a>=b)");
  global_dpd_->buf4_scm(&Z2, 0.0); /* this scm is necessary for cases with empty occpi or virtpi irreps */
  for(Gef=0; Gef < moinfo.nirreps; Gef++) {
    Gei = Gab = Gef; /* W and B are totally symmetric */
    for(Ge=0; Ge < moinfo.nirreps; Ge++) {
      Gf = Ge ^ Gef; Gi = Gf;  /* T1 is totally symmetric */
      B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.virtpi[Gf],B.params->coltot[Gef]);
      Z2.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gi],Z2.params->coltot[Gei]);
      nrows = moinfo.occpi[Gi]; ncols = Z2.params->coltot[Gef]; nlinks = moinfo.virtpi[Gf];
      if(nrows && ncols && nlinks) {
	for(E=0; E < moinfo.virtpi[Ge]; E++) {
	  e = moinfo.vir_off[Ge] + E;
	  global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e], moinfo.virtpi[Gf]);
	  C_DGEMM('n','n',nrows,ncols,nlinks,0.5,T1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
		  0.0,Z2.matrix[Gei][0],ncols);
	  global_dpd_->buf4_mat_irrep_wrt_block(&Z2, Gei, Z2.row_offset[Gei][e], moinfo.occpi[Gi]);
	}
      }
      global_dpd_->free_dpd_block(B.matrix[Gef], moinfo.virtpi[Gf], B.params->coltot[Gef]);
      global_dpd_->free_dpd_block(Z2.matrix[Gef], moinfo.occpi[Gi], Z2.params->coltot[Gei]);
    }
  }
  global_dpd_->buf4_close(&Z2);
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 11, 5, 11, 8, 0, "Z1(ei,a>=b)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 5, 11, 9, 0, "Z2(ei,a>=b)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
  global_dpd_->buf4_axpy(&Z1, &W, 1);
  global_dpd_->buf4_axpy(&Z2, &W, 1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&Z1);

  if(params.print & 2) outfile->Printf( "done.\n");

  /** Terms IIIc + IIId + IV **/
  if(params.print & 2) {
    outfile->Printf( "\tD*T1*T2 + E*T2 -> Wabei...");

  }
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (nM,eI)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, rspq, 11, 0, "WMnIe (eI,nM)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 11, 0, 11, 0, 0, "WMnIe (eI,nM)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  /*   dpd_contract444(&Z, &T, &W, 1, 1, 1, 1); */
  for(Gei=0; Gei < moinfo.nirreps; Gei++) {
    Gab = Gnm = Gei; /* Everything is totally symmetric here */
    nrows = T.params->rowtot[Gnm];
    ncols = T.params->coltot[Gab];
    if(nrows && ncols) {
      global_dpd_->buf4_mat_irrep_init(&Z, Gei);
      global_dpd_->buf4_mat_irrep_rd(&Z, Gei);
      global_dpd_->buf4_mat_irrep_init(&T, Gnm);
      global_dpd_->buf4_mat_irrep_rd(&T, Gnm);
      global_dpd_->buf4_mat_irrep_row_init(&W, Gei);
      for(ei=0; ei < W.params->rowtot[Gei]; ei++) {
	global_dpd_->buf4_mat_irrep_row_rd(&W, Gei, ei);
	C_DGEMV('t',nrows,ncols,1,T.matrix[Gei][0],ncols,Z.matrix[Gei][ei],1,
		1,W.matrix[Gei][0],1);
	global_dpd_->buf4_mat_irrep_row_wrt(&W, Gei, ei);
      }
      global_dpd_->buf4_mat_irrep_row_close(&W, Gei);
      global_dpd_->buf4_mat_irrep_close(&T, Gnm);
      global_dpd_->buf4_mat_irrep_close(&Z, Gei);
    }
  }
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);
  if(params.print & 2) outfile->Printf( "done.\n");

  /*** Terms IIIB + V ***/
  /* Wabei <-- Z1(mi,fb) <ma|fe> - t(mi,fb) <ma|ef> - tau(mi,af) <mb|ef> */
  if(params.print & 2) {
    outfile->Printf( "\t(T2+T1*T1) * F -> Wabei...");

  }
  build_Z1(); /* Z1(ib,mf) = 2 t(mi,fb) - t(mi,bf) - t(m,b) t(i,f) */

  /* The default algorithm for terms IIIb+V requires storage of two additional
  ** ovvv quantities on disk besides the <ia|bc> integrals and the Wabei target. */
  if(!params.wabei_lowdisk) {

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, prsq, 10, 5, "F <ia|bc> (ib,ca)");
    global_dpd_->buf4_close(&F);

    /* can we run these contractions fully in core? */
    incore = 1;
    core_total = 0;
    for(h=0; h < moinfo.nirreps; h++) {
      coltot = F.params->coltot[h];
      if(coltot)
	maxrows = DPD_BIGNUM/coltot;
      if(maxrows < 1) {
	outfile->Printf( "\nWabei_RHF Error: A single row of ovvv > 2 GW.\n");
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
      outfile->Printf( "\nUsing out-of-core algorithm for (T2+T1*T1)*F -> Wabei.\n");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "W(ib,ea)");
    global_dpd_->buf4_init(&F, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "F <ia|bc> (ib,ca)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ib,mf)");
    if(incore) global_dpd_->contract444(&Z, &F, &W, 0, 1, 1, 0);
    else ZFW(&Z, &F, &W, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&Z, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    if(incore) global_dpd_->contract444(&Z, &F, &W, 0, 1, -1, 1);
    else ZFW(&Z, &F, &W, -1, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, rpsq, 11, 5, "WAbEi (Ei,Ab)", 1);
    global_dpd_->buf4_close(&W);
    psio_close(PSIF_CC_TMP0, 0); /* delete the extra ovvv quantities on disk */
    psio_open(PSIF_CC_TMP0, PSIO_OPEN_NEW);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 10, 5, 10, 5, 0, "W(ia,eb)");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&Z, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tauIjAb (Ib,jA)");
    if(incore) global_dpd_->contract444(&Z, &F, &W, 0, 1, -1, 0);
    else ZFW(&Z, &F, &W, -1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, rpqs, 11, 5, "WAbEi (Ei,Ab)", 1);
    global_dpd_->buf4_close(&W);
    psio_close(PSIF_CC_TMP0, 0); /* delete the extra ovvv quantity on disk */
    psio_open(PSIF_CC_TMP0, PSIO_OPEN_NEW);
  }
  /* This is an alternative algorithm for terms IIIb+V that stores no additional
  ** ovvv terms beyond <ia|bc> integrals and the WAbEi target, but it's very slow. */
  else {
    if(params.print & 2)
      outfile->Printf( "\nUsing low-disk algorithm for (T2+T1*T1)*F -> Wabei.\n");

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
    /* prepare rows of W in advance */
    for(Gei=0; Gei < moinfo.nirreps; Gei++)
      global_dpd_->buf4_mat_irrep_row_init(&W, Gei);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ib,mf)");
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, rpsq, 0, 5, "Z1(mi,fb)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z1(mi,fb)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    WW1 = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    WW2 = (double ***) malloc(moinfo.nirreps * sizeof(double **));
    for(Gma=0; Gma < moinfo.nirreps; Gma++) {
      global_dpd_->buf4_mat_irrep_row_init(&F, Gma);
      for(ma=0; ma < F.params->rowtot[Gma]; ma++) {
	m = F.params->roworb[Gma][ma][0];
	a = F.params->roworb[Gma][ma][1];
	Gm = F.params->psym[m];
	Ga = F.params->qsym[a];
	global_dpd_->buf4_mat_irrep_row_rd(&F, Gma, ma);
	for(Gi=0; Gi < moinfo.nirreps; Gi++) {
	  Gmi = Gm ^ Gi;
	  global_dpd_->buf4_mat_irrep_row_init(&Z, Gmi);
	  global_dpd_->buf4_mat_irrep_row_init(&T2, Gmi);
	  global_dpd_->buf4_mat_irrep_row_init(&Tau, Gmi);
	  /* allocate space for WW1 target */
	  for(Gf=0; Gf < moinfo.nirreps; Gf++) {
	    Gb = Gf ^ Gmi; /* Z is totally symmetric */
	    Ge = Gf ^ Gma; /* F is totally symmetric */
	    WW1[Gb] = global_dpd_->dpd_block_matrix(moinfo.virtpi[Gb], moinfo.virtpi[Ge]);
	    WW2[Gb] = global_dpd_->dpd_block_matrix(moinfo.virtpi[Gb], moinfo.virtpi[Ge]);
	  }
	  for(I=0; I < moinfo.occpi[Gi]; I++) {
	    i = moinfo.occ_off[Gi] + I;
	    mi = Z.params->rowidx[m][i];
	    global_dpd_->buf4_mat_irrep_row_rd(&Z, Gmi, mi);
	    global_dpd_->buf4_mat_irrep_row_rd(&T2, Gmi, mi);
	    global_dpd_->buf4_mat_irrep_row_rd(&Tau, Gmi, mi);
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
		global_dpd_->buf4_mat_irrep_row_rd(&W, Gei, ei);
		for(BB=0; BB < moinfo.virtpi[Gb]; BB++) {
		  b = moinfo.vir_off[Gb] + BB;
		  ab = W.params->colidx[a][b];
		  ba = W.params->colidx[b][a];
		  W.matrix[Gei][0][ab] += WW1[Gb][BB][E];
		  W.matrix[Gei][0][ba] += WW2[Gb][BB][E];
		}
		global_dpd_->buf4_mat_irrep_row_wrt(&W, Gei, ei);
	      }
	    } /* Gf */
	  } /* I */
	  global_dpd_->buf4_mat_irrep_row_close(&Z, Gmi);
	  global_dpd_->buf4_mat_irrep_row_close(&T2, Gmi);
	  global_dpd_->buf4_mat_irrep_row_close(&Tau, Gmi);
	} /* Gi */
	  /* free W1 target */
	for(Gf=0; Gf < moinfo.nirreps; Gf++) {
	  Gb = Gf ^ Gmi; /* Z is totally symmetric */
	  Ge = Gf ^ Gma; /* F is totally symmetric */
	  global_dpd_->free_dpd_block(WW1[Gb],moinfo.virtpi[Gb], moinfo.virtpi[Ge]);
	  global_dpd_->free_dpd_block(WW2[Gb],moinfo.virtpi[Gb], moinfo.virtpi[Ge]);
	}
      } /* ma */
      global_dpd_->buf4_mat_irrep_row_close(&F, Gma);
    } /* Gma */
    free(WW1);
    free(WW2);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&F);
    for(Gei=0; Gei < moinfo.nirreps; Gei++)
      global_dpd_->buf4_mat_irrep_row_close(&W, Gei);
    global_dpd_->buf4_close(&W);
  }

  if(params.print & 2) outfile->Printf( "done.\n");

  /** Terms VI and VII **/

  /** t_in^bf  <Mn|Ef> + t_iN^bF <MN||EF> --> Z1_MEib **/
  if(params.print & 2) {
    outfile->Printf( "\tT1*(C+D*T2) -> Wabei...");

  }
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(ME,ib)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, -1, 1);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qrps, 11, 10, "Z(Ei,Mb)");
  global_dpd_->buf4_close(&Z);

  /** - t_M^A ( <Ei|Mb> + Z(Ei,Mb) ) --> W(Ei,Ab) **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_sort_axpy(&D, PSIF_CC_TMP0, spqr, 11, 10, "Z(Ei,Mb)", 1);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "Z(Ei,Mb)");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  /*  dpd_contract244(&T1, &Z, &W, 0, 2, 1, -1, 1); */
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  for(Gei=0; Gei < moinfo.nirreps; Gei++) {
    global_dpd_->buf4_mat_irrep_init(&Z, Gei);
    global_dpd_->buf4_mat_irrep_rd(&Z, Gei);
    global_dpd_->buf4_mat_irrep_row_init(&W, Gei);
    for(ei=0; ei < Z.params->rowtot[Gei]; ei++) {
      global_dpd_->buf4_mat_irrep_row_rd(&W, Gei, ei);
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
      global_dpd_->buf4_mat_irrep_row_wrt(&W, Gei, ei);
    }
    global_dpd_->buf4_mat_irrep_close(&Z, Gei);
  }
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /** - t_Ni^Af <mN|fE> --> Z2_mEiA **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(mE,iA)");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qrsp, 11, 11, "Z(Ei,Am)");
  global_dpd_->buf4_close(&Z);

  /** t_m^b ( - <mA|iE> + Z(mA,iE) ) --> Z2(Ab,Ei) **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 11, 11, 11, 11, 0, "C <ai|bj>");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(Ei,Am)");
  global_dpd_->buf4_axpy(&C, &Z, -1.0);
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  /*  dpd_contract424(&Z, &T1, &W, 3, 0, 0, 1, 1); */
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  for(Gei=0; Gei < moinfo.nirreps; Gei++) {
    global_dpd_->buf4_mat_irrep_init(&Z, Gei);
    global_dpd_->buf4_mat_irrep_rd(&Z, Gei);
    global_dpd_->buf4_mat_irrep_row_init(&W, Gei);
    for(ei=0; ei < Z.params->rowtot[Gei]; ei++) {
      global_dpd_->buf4_mat_irrep_row_rd(&W, Gei, ei);
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
      global_dpd_->buf4_mat_irrep_row_wrt(&W, Gei, ei);
    }
    global_dpd_->buf4_mat_irrep_close(&Z, Gei);
  }
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  if(params.print & 2) outfile->Printf( "done.\n");
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

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  global_dpd_->buf4_copy(&T2, PSIF_CC_TMP0, "Z1(ib,mf)");
  global_dpd_->buf4_close(&T2);

  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z1(ib,mf)");
  for(h=0; h < moinfo.nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Z1, h);
    global_dpd_->buf4_mat_irrep_rd(&Z1, h);

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
    global_dpd_->buf4_mat_irrep_wrt(&Z1, h);
    global_dpd_->buf4_mat_irrep_close(&Z1, h);
  }
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);

  global_dpd_->buf4_close(&Z1);
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
    global_dpd_->buf4_mat_irrep_init(Z, Gib);
    global_dpd_->buf4_mat_irrep_rd(Z, Gib);

    /* how many rows of F/W can we store in half the remaining core? */
    rows_per_bucket = dpd_memfree()/(2 * F->params->coltot[Gea]);
    if(rows_per_bucket > F->params->rowtot[Gib])
      rows_per_bucket = F->params->rowtot[Gib];
    nbuckets = (int) ceil((double) F->params->rowtot[Gib]/(double) rows_per_bucket);
    rows_left = F->params->rowtot[Gib] % rows_per_bucket;

    /* allocate space for the full buckets */
    global_dpd_->buf4_mat_irrep_init_block(F, Gib, rows_per_bucket);
    global_dpd_->buf4_mat_irrep_init_block(W, Gib, rows_per_bucket);

    ncols = W->params->coltot[Gea];

    for(m=0; m < (rows_left ? nbuckets-1 : nbuckets); m++) {
      W_row_start = m*rows_per_bucket;
      zero_arr(W->matrix[Gib][0], rows_per_bucket*ncols);

      if(beta == 1.0)
	global_dpd_->buf4_mat_irrep_rd_block(W, Gib, W_row_start, rows_per_bucket);

      for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {
	F_row_start = n*rows_per_bucket;
	global_dpd_->buf4_mat_irrep_rd_block(F, Gib, F_row_start, rows_per_bucket);

	nrows = rows_per_bucket;
	nlinks = rows_per_bucket;

	if(nrows && ncols && nlinks)
	  C_DGEMM('n', 'n', nrows, ncols, nlinks, alpha,
		  &(Z->matrix[Gib][W_row_start][F_row_start]), Z->params->coltot[Gmf],
		  F->matrix[Gmf][0], ncols, 1, W->matrix[Gib][0], ncols);
      }
      if(rows_left) {
	F_row_start = n*rows_per_bucket;
	global_dpd_->buf4_mat_irrep_rd_block(F, Gib, F_row_start, rows_left);

	nrows = rows_per_bucket;
	nlinks = rows_left;

	if(nrows && ncols && nlinks)
	  C_DGEMM('n', 'n', nrows, ncols, nlinks, alpha,
		  &(Z->matrix[Gib][W_row_start][F_row_start]), Z->params->coltot[Gmf],
		  F->matrix[Gmf][0], ncols, 1, W->matrix[Gib][0], ncols);
      }

      global_dpd_->buf4_mat_irrep_wrt_block(W, Gib, W_row_start, rows_per_bucket);
    }
    if(rows_left) {
      W_row_start = m*rows_per_bucket;
      zero_arr(W->matrix[Gib][0], rows_per_bucket*ncols);

      if(beta == 1.0)
	global_dpd_->buf4_mat_irrep_rd_block(W, Gib, W_row_start, rows_left);

      for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {
	F_row_start = n*rows_per_bucket;
	global_dpd_->buf4_mat_irrep_rd_block(F, Gib, F_row_start, rows_per_bucket);

	nrows = rows_left;
	nlinks = rows_per_bucket;

	if(nrows && ncols && nlinks)
	  C_DGEMM('n', 'n', nrows, ncols, nlinks, alpha,
		  &(Z->matrix[Gib][W_row_start][F_row_start]), Z->params->coltot[Gmf],
		  F->matrix[Gmf][0], ncols, 1, W->matrix[Gib][0], ncols);
      }
      if(rows_left) {
	F_row_start = n*rows_per_bucket;
	global_dpd_->buf4_mat_irrep_rd_block(F, Gib, F_row_start, rows_left);

	nrows = rows_left;
	nlinks = rows_left;

	if(nrows && ncols && nlinks)
	  C_DGEMM('n', 'n', nrows, ncols, nlinks, alpha,
		  &(Z->matrix[Gib][W_row_start][F_row_start]), Z->params->coltot[Gmf],
		  F->matrix[Gmf][0], ncols, 1, W->matrix[Gib][0], ncols);
      }
      global_dpd_->buf4_mat_irrep_wrt_block(W, Gib, W_row_start, rows_left);
    } /* m */

    global_dpd_->buf4_mat_irrep_close_block(F, Gib, rows_per_bucket);
    global_dpd_->buf4_mat_irrep_close_block(W, Gib, rows_per_bucket);

    global_dpd_->buf4_mat_irrep_close(Z, Gib);

  } /* Gib */
}

}} // namespace psi::cchbar
