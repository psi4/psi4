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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/

/*! \defgroup CCLAMBDA cclambda: Coupled-Cluster Lambda Equations */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void halftrans(dpdbuf4 *Buf1, int dpdnum1, dpdbuf4 *Buf2, int dpdnum2, double ***C, int nirreps,
	       int **mo_row, int **so_row, int *mospi, int *sospi, int type, double alpha, double beta);

void AO_contribute(int p, int q, int r, int s, double value,
		   dpdbuf4 *tau1_AO, dpdbuf4 *tau2_AO, int anti);

void BL2_AO(int L_irr)
{
  int h, nirreps, i, Gc, Gd, Ga, Gb, ij;
  double ***C, **X;
  int *orbspi, *virtpi;
  int **T2_cd_row_start, **T2_pq_row_start, offset, cd, pq;
  dpdbuf4 tau, t2, tau1_AO, tau2_AO;
  psio_address next;
  struct iwlbuf InBuf;
  int idx, p, q, r, s, filenum;
  int lastbuf;
  double value, tolerance=1e-14;
  Value *valptr;
  Label *lblptr;

  nirreps = moinfo.nirreps;
  orbspi = moinfo.orbspi;
  virtpi = moinfo.virtpi;
  C = moinfo.C;

  T2_cd_row_start = init_int_matrix(nirreps,nirreps);
  for(h=0; h < nirreps; h++) {
    for(Gc=0,offset=0; Gc < nirreps; Gc++) {
      Gd = Gc ^ h;
      T2_cd_row_start[h][Gc] = offset;
      offset += virtpi[Gc] * virtpi[Gd];
    }
  }

  T2_pq_row_start = init_int_matrix(nirreps,nirreps);
  for(h=0; h < nirreps; h++) {
    for(Gc=0,offset=0; Gc < nirreps; Gc++) {
      Gd = Gc ^ h;
      T2_pq_row_start[h][Gc] = offset;
      offset += orbspi[Gc] * orbspi[Gd];
    }
  }

  /************************************* AA *****************************************/

  dpd_set_default(1);
  global_dpd_->buf4_init(&tau1_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIJPQ (1)");
  global_dpd_->buf4_scm(&tau1_AO, 0.0);

  dpd_set_default(0);
  global_dpd_->buf4_init(&tau, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "LIJAB");

  halftrans(&tau, 0, &tau1_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start,
	    virtpi, orbspi, 0, 1.0, 0.0);

  global_dpd_->buf4_close(&tau);
  global_dpd_->buf4_close(&tau1_AO);

  dpd_set_default(1);
  global_dpd_->buf4_init(&tau1_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIJPQ (1)");
  global_dpd_->buf4_init(&tau2_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIJPQ (2)");
  global_dpd_->buf4_scm(&tau2_AO, 0.0);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&tau1_AO, h);
    global_dpd_->buf4_mat_irrep_rd(&tau1_AO, h);
    global_dpd_->buf4_mat_irrep_init(&tau2_AO, h);
  }

  iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

  lblptr = InBuf.labels;
  valptr = InBuf.values;
  lastbuf = InBuf.lastbuf;

  for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
    p = abs((int) lblptr[idx++]);
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    value = (double) valptr[InBuf.idx];

    /*    outfile->Printf( "<%d %d %d %d = %20.10lf\n", p, q, r, s, value); */

    AO_contribute(p, q, r, s, value, &tau1_AO, &tau2_AO, 1);

  }
  while(!lastbuf) {
    iwl_buf_fetch(&InBuf);
    lastbuf = InBuf.lastbuf;
    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = abs((int) lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      /*      outfile->Printf( "<%d %d %d %d = %20.10lf\n", p, q, r, s, value); */

      AO_contribute(p, q, r, s, value, &tau1_AO, &tau2_AO, 1);

    }
  }

  iwl_buf_close(&InBuf, 1);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_wrt(&tau2_AO, h);
    global_dpd_->buf4_mat_irrep_close(&tau2_AO, h);
    global_dpd_->buf4_mat_irrep_close(&tau1_AO, h);
  }
  global_dpd_->buf4_close(&tau1_AO);
  global_dpd_->buf4_close(&tau2_AO);


  dpd_set_default(1);
  global_dpd_->buf4_init(&tau2_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIJPQ (2)");

  dpd_set_default(0);
  global_dpd_->buf4_init(&t2, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "New LIJAB");

  halftrans(&t2, 0, &tau2_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start,
	    virtpi, orbspi, 1, 0.5, 1.0);

  global_dpd_->buf4_close(&t2);
  global_dpd_->buf4_close(&tau2_AO);

  /************************************* BB *****************************************/

  dpd_set_default(1);
  global_dpd_->buf4_init(&tau1_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "Lijpq (1)");
  global_dpd_->buf4_scm(&tau1_AO, 0.0);

  dpd_set_default(0);
  global_dpd_->buf4_init(&tau, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "Lijab");

  halftrans(&tau, 0, &tau1_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start,
	    virtpi, orbspi, 0, 1.0, 0.0);

  global_dpd_->buf4_close(&tau);
  global_dpd_->buf4_close(&tau1_AO);

  dpd_set_default(1);
  global_dpd_->buf4_init(&tau1_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "Lijpq (1)");
  global_dpd_->buf4_init(&tau2_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "Lijpq (2)");
  global_dpd_->buf4_scm(&tau2_AO, 0.0);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&tau1_AO, h);
    global_dpd_->buf4_mat_irrep_rd(&tau1_AO, h);
    global_dpd_->buf4_mat_irrep_init(&tau2_AO, h);
  }

  iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

  lblptr = InBuf.labels;
  valptr = InBuf.values;
  lastbuf = InBuf.lastbuf;

  for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
    p = abs((int) lblptr[idx++]);
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    value = (double) valptr[InBuf.idx];

    /*    outfile->Printf( "<%d %d %d %d = %20.10lf\n", p, q, r, s, value); */

    AO_contribute(p, q, r, s, value, &tau1_AO, &tau2_AO, 1);

  }
  while(!lastbuf) {
    iwl_buf_fetch(&InBuf);
    lastbuf = InBuf.lastbuf;
    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = abs((int) lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      /*      outfile->Printf( "<%d %d %d %d = %20.10lf\n", p, q, r, s, value); */

      AO_contribute(p, q, r, s, value, &tau1_AO, &tau2_AO, 1);

    }
  }

  iwl_buf_close(&InBuf, 1);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_wrt(&tau2_AO, h);
    global_dpd_->buf4_mat_irrep_close(&tau2_AO, h);
    global_dpd_->buf4_mat_irrep_close(&tau1_AO, h);
  }
  global_dpd_->buf4_close(&tau1_AO);
  global_dpd_->buf4_close(&tau2_AO);


  dpd_set_default(1);
  global_dpd_->buf4_init(&tau2_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "Lijpq (2)");

  dpd_set_default(0);
  global_dpd_->buf4_init(&t2, PSIF_CC_LAMBDA, L_irr, 0, 5, 2, 7, 0, "New Lijab");

  halftrans(&t2, 0, &tau2_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start,
	    virtpi, orbspi, 1, 0.5, 1.0);

  global_dpd_->buf4_close(&t2);
  global_dpd_->buf4_close(&tau2_AO);

  /************************************* AB *****************************************/

  dpd_set_default(1);
  global_dpd_->buf4_init(&tau1_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjPq (1)");
  global_dpd_->buf4_scm(&tau1_AO, 0.0);

  dpd_set_default(0);
  global_dpd_->buf4_init(&tau, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");

  halftrans(&tau, 0, &tau1_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start,
	    virtpi, orbspi, 0, 1.0, 0.0);

  global_dpd_->buf4_close(&tau);
  global_dpd_->buf4_close(&tau1_AO);

  dpd_set_default(1);
  global_dpd_->buf4_init(&tau1_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjPq (1)");
  global_dpd_->buf4_init(&tau2_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjPq (2)");
  global_dpd_->buf4_scm(&tau2_AO, 0.0);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&tau1_AO, h);
    global_dpd_->buf4_mat_irrep_rd(&tau1_AO, h);
    global_dpd_->buf4_mat_irrep_init(&tau2_AO, h);
  }

  iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

  lblptr = InBuf.labels;
  valptr = InBuf.values;
  lastbuf = InBuf.lastbuf;

  for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
    p = abs((int) lblptr[idx++]);
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    value = (double) valptr[InBuf.idx];

    /*    outfile->Printf( "<%d %d %d %d = %20.10lf\n", p, q, r, s, value); */

    AO_contribute(p, q, r, s, value, &tau1_AO, &tau2_AO, 0);

  }
  while(!lastbuf) {
    iwl_buf_fetch(&InBuf);
    lastbuf = InBuf.lastbuf;
    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = abs((int) lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      /*      outfile->Printf( "<%d %d %d %d = %20.10lf\n", p, q, r, s, value); */

      AO_contribute(p, q, r, s, value, &tau1_AO, &tau2_AO, 0);

    }
  }

  iwl_buf_close(&InBuf, 1);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_wrt(&tau2_AO, h);
    global_dpd_->buf4_mat_irrep_close(&tau2_AO, h);
    global_dpd_->buf4_mat_irrep_close(&tau1_AO, h);
  }
  global_dpd_->buf4_close(&tau1_AO);
  global_dpd_->buf4_close(&tau2_AO);


  dpd_set_default(1);
  global_dpd_->buf4_init(&tau2_AO, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjPq (2)");

  dpd_set_default(0);
  global_dpd_->buf4_init(&t2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

  halftrans(&t2, 0, &tau2_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start,
	    virtpi, orbspi, 1, 1.0, 1.0);

  global_dpd_->buf4_close(&t2);
  global_dpd_->buf4_close(&tau2_AO);

  free(T2_cd_row_start);
  free(T2_pq_row_start);

  /* Reset the default dpd back to 0 --- this stuff gets really ugly */
  dpd_set_default(0);
}

void AO_contribute(int p, int q, int r, int s, double value, dpdbuf4
		   *tau1_AO, dpdbuf4 *tau2_AO, int anti)
{
  int Gp, Gq, Gr, Gs, Gpr, Gps, Gqr, Gqs, Grp, Gsp, Grq, Gsq;
  int pr, ps, qr, qs, rp, rq, sp, sq, pq, rs;
  int row;

  Gp = tau1_AO->params->rsym[p];
  Gq = tau1_AO->params->rsym[q];
  Gr = tau1_AO->params->rsym[r];
  Gs = tau1_AO->params->rsym[s];

  pq = tau1_AO->params->colidx[p][q];  rs = tau1_AO->params->colidx[r][s];

  if(p!=q && r!=s) {

    /* (pq|rs) */
    Gpr = Gp ^ Gr;
    pr = tau1_AO->params->colidx[p][r];
    qs = tau1_AO->params->colidx[q][s];
    sq = tau1_AO->params->colidx[s][q];

    for(row=0; row < tau1_AO->params->rowtot[Gpr]; row++) {
      tau2_AO->matrix[Gpr][row][pr] += value * tau1_AO->matrix[Gpr][row][qs];
      if(anti) tau2_AO->matrix[Gpr][row][pr] -= value * tau1_AO->matrix[Gpr][row][sq];
    }

    /* (pq|sr) */
    Gps = Gp ^ Gs;
    ps = tau1_AO->params->colidx[p][s];
    qr = tau1_AO->params->colidx[q][r];
    rq = tau1_AO->params->colidx[r][q];

    for(row=0; row < tau1_AO->params->rowtot[Gps]; row++) {
      tau2_AO->matrix[Gps][row][ps] += value * tau1_AO->matrix[Gps][row][qr];
      if(anti) tau2_AO->matrix[Gps][row][ps] -= value * tau1_AO->matrix[Gps][row][rq];
    }

    /* (qp|rs) */
    Gqr = Gq ^ Gr;
    qr = tau1_AO->params->colidx[q][r];
    ps = tau1_AO->params->colidx[p][s];
    sp = tau1_AO->params->colidx[s][p];

    for(row=0; row < tau1_AO->params->rowtot[Gqr]; row++) {
      tau2_AO->matrix[Gqr][row][qr] += value * tau1_AO->matrix[Gqr][row][ps];
      if(anti) tau2_AO->matrix[Gqr][row][qr] -= value * tau1_AO->matrix[Gqr][row][sp];
    }

    /* (qp|sr) */
    Gqs = Gq ^ Gs;
    qs = tau1_AO->params->colidx[q][s];
    pr = tau1_AO->params->colidx[p][r];
    rp = tau1_AO->params->colidx[r][p];

    for(row=0; row < tau1_AO->params->rowtot[Gqs]; row++) {
      tau2_AO->matrix[Gqs][row][qs] += value * tau1_AO->matrix[Gqs][row][pr];
      if(anti) tau2_AO->matrix[Gqs][row][qs] -= value * tau1_AO->matrix[Gqs][row][rp];
    }

    if(pq != rs) {
      /* (rs|pq) */
      Grp = Gp ^ Gr;
      rp = tau1_AO->params->colidx[r][p];
      sq = tau1_AO->params->colidx[s][q];
      qs = tau1_AO->params->colidx[q][s];

      for(row=0; row < tau1_AO->params->rowtot[Grp]; row++) {
	tau2_AO->matrix[Grp][row][rp] += value * tau1_AO->matrix[Grp][row][sq];
	if(anti) tau2_AO->matrix[Grp][row][rp] -= value * tau1_AO->matrix[Grp][row][qs];
      }

      /* (sr|pq) */
      Gsp = Gp ^ Gs;
      sp = tau1_AO->params->colidx[s][p];
      qr = tau1_AO->params->colidx[q][r];
      rq = tau1_AO->params->colidx[r][q];

      for(row=0; row < tau1_AO->params->rowtot[Gsp]; row++) {
	tau2_AO->matrix[Gsp][row][sp] += value * tau1_AO->matrix[Gsp][row][rq];
	if(anti) tau2_AO->matrix[Gsp][row][sp] -= value * tau1_AO->matrix[Gsp][row][qr];
      }

      /* (rs|qp) */
      Grq = Gq ^ Gr;
      rq = tau1_AO->params->colidx[r][q];
      ps = tau1_AO->params->colidx[p][s];
      sp = tau1_AO->params->colidx[s][p];

      for(row=0; row < tau1_AO->params->rowtot[Grq]; row++) {
	tau2_AO->matrix[Grq][row][rq] += value * tau1_AO->matrix[Grq][row][sp];
	if(anti) tau2_AO->matrix[Grq][row][rq] -= value * tau1_AO->matrix[Grq][row][ps];
      }

      /* (sr|qp) */
      Gsq = Gq ^ Gs;
      sq = tau1_AO->params->colidx[s][q];
      pr = tau1_AO->params->colidx[p][r];
      rp = tau1_AO->params->colidx[r][p];

      for(row=0; row < tau1_AO->params->rowtot[Gsq]; row++) {
	tau2_AO->matrix[Gsq][row][sq] += value * tau1_AO->matrix[Gsq][row][rp];
	if(anti) tau2_AO->matrix[Gsq][row][sq] -= value * tau1_AO->matrix[Gsq][row][pr];
      }
    }

  }
  else if(p!=q && r==s) {

    /* (pq|rs) */
    Gpr = Gp ^ Gr;
    pr = tau1_AO->params->colidx[p][r];
    qs = tau1_AO->params->colidx[q][s];
    sq = tau1_AO->params->colidx[s][q];

    for(row=0; row < tau1_AO->params->rowtot[Gpr]; row++) {
      tau2_AO->matrix[Gpr][row][pr] += value * tau1_AO->matrix[Gpr][row][qs];
      if(anti) tau2_AO->matrix[Gpr][row][pr] -= value * tau1_AO->matrix[Gpr][row][sq];
    }

    /* (qp|rs) */
    Gqr = Gq ^ Gr;
    qr = tau1_AO->params->colidx[q][r];
    ps = tau1_AO->params->colidx[p][s];
    sp = tau1_AO->params->colidx[s][p];

    for(row=0; row < tau1_AO->params->rowtot[Gqr]; row++) {
      tau2_AO->matrix[Gqr][row][qr] += value * tau1_AO->matrix[Gqr][row][ps];
      if(anti) tau2_AO->matrix[Gqr][row][qr] -= value * tau1_AO->matrix[Gqr][row][sp];
    }

    if(pq != rs) {

      /* (rs|pq) */
      Grp = Gp ^ Gr;
      rp = tau1_AO->params->colidx[r][p];
      sq = tau1_AO->params->colidx[s][q];
      qs = tau1_AO->params->colidx[q][s];

      for(row=0; row < tau1_AO->params->rowtot[Grp]; row++) {
	tau2_AO->matrix[Grp][row][rp] += value * tau1_AO->matrix[Grp][row][sq];
	if(anti) tau2_AO->matrix[Grp][row][rp] -= value * tau1_AO->matrix[Grp][row][qs];
      }

      /* (rs|qp) */
      Grq = Gq ^ Gr;
      rq = tau1_AO->params->colidx[r][q];
      ps = tau1_AO->params->colidx[p][s];
      sp = tau1_AO->params->colidx[s][p];

      for(row=0; row < tau1_AO->params->rowtot[Grq]; row++) {
	tau2_AO->matrix[Grq][row][rq] += value * tau1_AO->matrix[Grq][row][sp];
	if(anti) tau2_AO->matrix[Grq][row][rq] -= value * tau1_AO->matrix[Grq][row][ps];
      }
    }

  }

  else if(p==q && r!=s) {

    /* (pq|rs) */
    Gpr = Gp ^ Gr;
    pr = tau1_AO->params->colidx[p][r];
    qs = tau1_AO->params->colidx[q][s];
    sq = tau1_AO->params->colidx[s][q];

    for(row=0; row < tau1_AO->params->rowtot[Gpr]; row++) {
      tau2_AO->matrix[Gpr][row][pr] += value * tau1_AO->matrix[Gpr][row][qs];
      if(anti) tau2_AO->matrix[Gpr][row][pr] -= value * tau1_AO->matrix[Gpr][row][sq];
    }

    /* (pq|sr) */
    Gps = Gp ^ Gs;
    ps = tau1_AO->params->colidx[p][s];
    qr = tau1_AO->params->colidx[q][r];
    rq = tau1_AO->params->colidx[r][q];

    for(row=0; row < tau1_AO->params->rowtot[Gps]; row++) {
      tau2_AO->matrix[Gps][row][ps] += value * tau1_AO->matrix[Gps][row][qr];
      if(anti) tau2_AO->matrix[Gps][row][ps] -= value * tau1_AO->matrix[Gps][row][rq];
    }

    if(pq != rs) {

      /* (rs|pq) */
      Grp = Gp ^ Gr;
      rp = tau1_AO->params->colidx[r][p];
      sq = tau1_AO->params->colidx[s][q];
      qs = tau1_AO->params->colidx[q][s];

      for(row=0; row < tau1_AO->params->rowtot[Grp]; row++) {
	tau2_AO->matrix[Grp][row][rp] += value * tau1_AO->matrix[Grp][row][sq];
	if(anti) tau2_AO->matrix[Grp][row][rp] -= value * tau1_AO->matrix[Grp][row][qs];
      }

      /* (sr|pq) */
      Gsp = Gp ^ Gs;
      sp = tau1_AO->params->colidx[s][p];
      qr = tau1_AO->params->colidx[q][r];
      rq = tau1_AO->params->colidx[r][q];

      for(row=0; row < tau1_AO->params->rowtot[Gsp]; row++) {
	tau2_AO->matrix[Gsp][row][sp] += value * tau1_AO->matrix[Gsp][row][rq];
	if(anti) tau2_AO->matrix[Gsp][row][sp] -= value * tau1_AO->matrix[Gsp][row][qr];
      }
    }

  }

  else if(p==q && r==s) {

    /* (pq|rs) */
    Gpr = Gp ^ Gr;
    pr = tau1_AO->params->colidx[p][r];
    qs = tau1_AO->params->colidx[q][s];
    sq = tau1_AO->params->colidx[s][q];

    for(row=0; row < tau1_AO->params->rowtot[Gpr]; row++) {
      tau2_AO->matrix[Gpr][row][pr] += value * tau1_AO->matrix[Gpr][row][qs];
      if(anti) tau2_AO->matrix[Gpr][row][pr] -= value * tau1_AO->matrix[Gpr][row][sq];
    }

    if(pq != rs) {

      /* (rs|pq) */
      Grp = Gp ^ Gr;
      rp = tau1_AO->params->colidx[r][p];
      sq = tau1_AO->params->colidx[s][q];
      qs = tau1_AO->params->colidx[q][s];

      for(row=0; row < tau1_AO->params->rowtot[Grp]; row++) {
	tau2_AO->matrix[Grp][row][rp] += value * tau1_AO->matrix[Grp][row][sq];
	if(anti) tau2_AO->matrix[Grp][row][rp] -= value * tau1_AO->matrix[Grp][row][qs];
      }

    }

  }
  return;
}

}} // namespace psi::cclambda
