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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"

namespace psi {
namespace cctriples {

/* T3_UHF_AAA_abc(): Computes all connected and disconnected
 * T3(IJK,ABC) * amplitudes for a given A, B, C combination for input
 * C2, F, and E * intermediates.  This function will work for AAA or
 * BBB spin cases, * with either RHF/ROHF or UHF orbitals.  This is
 * analogous to the T3_UHF_AAA() * function *
 *
 * Arguments: *
 *
 *   double ***W: The target connected triples amplitudes in an *
 * nirreps x I x JK array.  The memory for this must be allocated *
 * externally. *
 *
 *   double ***V: The target disconnected triples amplitudes (if *
 * requested) in an nirreps x I x JK array.  The memory for this *
 * must be allocated externally. *
 *
 *   int disc: Boolean: 1 == computed disconnected triples; 0 == don't *
 *
 *   int nirreps: Number of irreps. *
 *
 *   int A: Absolute index of orbital A. *
 *
 *   int Ga: Irrep of A. *
 *
 *   int B: Absolute index of orbital B. *
 *
 *   int Gb: Irrep of B. *
 *
 *   int C: Absolute index of orbital C. *
 *
 *   int Gc: Irrep of C. *
 *
 *   dpdbuf4 *C2: Pointer to dpd buffer for double excitation amps, *
 * ordered (AB, IJ).  N.B. This is transposed relative to the
 * ordering *   used in the T3_UHF_AAA() routine. *
 *
 *   dpdbuf4 *F: Pointer to dpd buffer for three-virtual-index *
 * intermediate, ordered (BC,IA). N.B. This is transposed relative *
 * to the ordering used in the T3_UHF_AAA() routine *
 *
 *   dpdbuf4 *E: Pointer to dpd buffer for three-occupied-index *
 * intermediate, ordered (IJ,KA). *
 *
 *   dpdfile *C1: If disconnected T3's are requested, pointer to dpd *
 * buffer for single-excitation amps. *
 *
 *   dpdbuf4 *D: If disconnected T3's are requested, pointer to dpd *
 * buffer for <IJ||ab> integrals. *
 *
 *   dpdfile2 *fIA: Pointer to the dpd file2 for the occ-vir block of *
 * the Fock matrix (or other appropriate one-electron operator). *
 *
 *   dpdfile2 *fIJ: Pointer to the dpd file2 for the occ-occ block of *
 * the Fock matrix (or other appropriate one-electron operator). *   *
 * dpdfile2 *fAB: Pointer to the dpd file2 for the vir-vir block of *
 * the Fock matrix (or other appropriate one-electron operator). *   *
 * int *occpi: Number of occupied orbitals per irrep lookup array. *
 *
 *   int *occ_off: Offset lookup for translating between absolute and *
 * relative orbital indices for occupied space. *
 *
 *   int *virtpi: Number of virtual orbitals per irrep lookup array. *
 *
 *   int *vir_off: Offset lookup for translating between absolute and *
 * relative orbital indices for virtual space. *
 *
 *   double omega: a constant to add to the final denominators - *
 * needed for CC3 EOM *
 *
 * ACS, December 2008. */

void T3_UHF_AAA_abc(double ***W, double ***V, int disc, int nirreps, int A,
    int Ga, int B, int Gb, int C, int Gc, dpdbuf4 * C2, dpdbuf4 * F,
    dpdbuf4 * E, dpdfile2 * C1, dpdbuf4 * D, dpdfile2 * fIA, dpdfile2 * fIJ,
    dpdfile2 * fAB, int *occpi, int *occ_off, int *virtpi, int *vir_off,
    double omega) {
  int h;
  int i, j, k;
  int ij, ji, ik, ki, jk, kj;
  int Gij, Gji, Gik, Gki, Gjk, Gkj, Gijk;
  int Gba, Gbd, Gabc;
  int Gi, Gj, Gk;
  int Gd, Gl;
  int Gid, Gjd, Gkd;
  int Gab, Gcb, Gca, Gbc, Gac, Gdc, Gcd, Gad;
  int Gla, Glb, Glc;
  int Gil, Gjl, Gkl;
  int Gal, Gbl, Gcl;
  int a, b, c, I, J, K;
  int ab, ba, bc, ac, ca, cb;
  int cd, bd, ad;
  int id, jd, kd;
  int la, lb, lc;
  int il, jl, kl;
  int al, bl, cl;
  int nrows, ncols, nlinks;
  double dabc, denom;
  double ***W2;
  int GE, GF, GC, GX3;
  double t_ia, t_ib, t_ic, t_ja, t_jb, t_jc, t_ka, t_kb, t_kc;
  double f_ia, f_ib, f_ic, f_ja, f_jb, f_jc, f_ka, f_kb, f_kc;
  double t_jkbc, t_jkac, t_jkba, t_ikbc, t_ikac, t_ikba, t_jibc, t_jiac, t_jiba;
  double D_jkbc, D_jkac, D_jkba, D_ikbc, D_ikac, D_ikba, D_jibc, D_jiac, D_jiba;

  GC = C2->file.my_irrep;
  /* F and E are assumed to have same irrep */
  GF = GE = F->file.my_irrep;
  GX3 = GC ^ GF;

  global_dpd_->file2_mat_init(C1);
  global_dpd_->file2_mat_rd(C1);

  global_dpd_->file2_mat_init(fIJ);
  global_dpd_->file2_mat_init(fAB);
  global_dpd_->file2_mat_rd(fIJ);
  global_dpd_->file2_mat_rd(fAB);
  if (disc) {
    global_dpd_->file2_mat_init(fIA);
    global_dpd_->file2_mat_rd(fIA);
  }

  for (h = 0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(F, h);
    global_dpd_->buf4_mat_irrep_rd(F, h);
  }

  a = A - vir_off[Ga];
  b = B - vir_off[Gb];
  c = C - vir_off[Gc];

  Gab = Gba = Ga ^ Gb;
  Gac = Gca = Ga ^ Gc;
  Gbc = Gcb = Gb ^ Gc;
  Gabc = Ga ^ Gb ^ Gc;

  ab = F->params->rowidx[A][B];
  ac = F->params->rowidx[A][C];
  ba = F->params->rowidx[B][A];
  bc = F->params->rowidx[B][C];
  ca = F->params->rowidx[C][A];
  cb = F->params->rowidx[C][B];

  dabc = 0.0;
  if (fAB->params->rowtot[Ga])
    dabc -= fAB->matrix[Ga][a][a];
  if (fAB->params->rowtot[Gb])
    dabc -= fAB->matrix[Gb][b][b];
  if (fAB->params->rowtot[Gc])
    dabc -= fAB->matrix[Gc][c][c];

  W2 = (double ***) malloc(nirreps * sizeof(double **));

  for (Gij = 0; Gij < nirreps; ++Gij) {
    Gk = Gij ^ Gabc ^ GX3; /* changed */
    W2[Gij] = global_dpd_->dpd_block_matrix(C2->params->coltot[Gij], occpi[Gk]);
    if (C2->params->coltot[Gij] && occpi[Gk]) {
      memset(W[Gij][0], 0, C2->params->coltot[Gij] * occpi[Gk] * sizeof(double));
      memset(W2[Gij][0], 0, C2->params->coltot[Gij] * occpi[Gk]
          * sizeof(double));
      if (disc)
        memset(V[Gij][0], 0, C2->params->coltot[Gij] * occpi[Gk]
            * sizeof(double));
    }
  } /* Gij */

  /*=======================================*/
  /* The terms that contribute as W[ij][k] */
  /*=======================================*/

  for (Gd = 0; Gd < nirreps; Gd++) {
    /* -t_ijcd * F_kdab, rewritten as -t_cdij * F_abkd */
    Gcd = Gc ^ Gd;
    Gk = Gd ^ Gab ^ GF;
    Gij = Gcd ^ GC;

    kd = F->col_offset[Gab][Gk];
    cd = C2->row_offset[Gcd][C];

    ncols = occpi[Gk];
    nrows = C2->params->coltot[Gij];
    nlinks = virtpi[Gd];

    if (nrows && ncols && nlinks) {
      C2->matrix[Gcd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(C2, Gcd, cd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, C2->matrix[Gcd][0], nrows,
          &(F->matrix[Gab][ab][kd]), nlinks, 1.0, W[Gij][0], ncols);
      global_dpd_->free_dpd_block(C2->matrix[Gcd], nlinks, nrows);
    }

    /* -t_ijbd * F_kdca, rewritten as  -t_bdij * F_cakd */
    Gbd = Gb ^ Gd;
    Gk = Gd ^ Gca ^ GF;
    Gij = Gbd ^ GC;

    kd = F->col_offset[Gca][Gk];
    bd = C2->row_offset[Gbd][B];

    ncols = occpi[Gk];
    nrows = C2->params->coltot[Gij];
    nlinks = virtpi[Gd];

    if (nrows && ncols && nlinks) {
      C2->matrix[Gbd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(C2, Gbd, bd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, C2->matrix[Gbd][0], nrows,
          &(F->matrix[Gca][ca][kd]), nlinks, 1.0, W[Gij][0], ncols);
      global_dpd_->free_dpd_block(C2->matrix[Gbd], nlinks, nrows);
    }

    /* +t_ijad * F_kdcb, rewritten as t_adij * F_cbkd */
    Gad = Ga ^ Gd;
    Gk = Gd ^ Gcb ^ GF;
    Gij = Gad ^ GC;

    kd = F->col_offset[Gcb][Gk];
    ad = C2->row_offset[Gad][A];

    ncols = occpi[Gk];
    nrows = C2->params->coltot[Gij];
    nlinks = virtpi[Gd];

    if (nrows && ncols && nlinks) {
      C2->matrix[Gad] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(C2, Gad, ad, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, C2->matrix[Gad][0], nrows,
          &(F->matrix[Gcb][cb][kd]), nlinks, 1.0, W[Gij][0], ncols);
      global_dpd_->free_dpd_block(C2->matrix[Gad], nlinks, nrows);
    }
  } /* Gd */

  for (Gl = 0; Gl < nirreps; ++Gl) {
    /* +t_klab * E_jilc, rewritten as -t_abkl * E_clij */
    Gcl = Gc ^ Gl;
    Gk = Gab ^ Gl ^ GC;
    Gij = Gcl ^ GE;

    global_dpd_->buf4_mat_irrep_init(C2, Gab);
    global_dpd_->buf4_mat_irrep_rd(C2, Gab);

    cl = E->row_offset[Gcl][C];
    kl = C2->col_offset[Gab][Gk];

    ncols = occpi[Gk];
    nrows = E->params->coltot[Gij];
    nlinks = occpi[Gl];

    if (nrows && ncols && nlinks) {
      E->matrix[Gcl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(E, Gcl, cl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, E->matrix[Gcl][0], nrows,
          &(C2->matrix[Gab][ab][kl]), nlinks, 1.0, W[Gij][0], ncols);
      global_dpd_->free_dpd_block(E->matrix[Gcl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(C2, Gab);

    /* +t_klca * E_jilb, rewritten as -t_cakl * E_blij */
    Gbl = Gb ^ Gl;
    Gk = Gca ^ Gl ^ GC;
    Gij = Gbl ^ GE;

    global_dpd_->buf4_mat_irrep_init(C2, Gca);
    global_dpd_->buf4_mat_irrep_rd(C2, Gca);

    bl = E->row_offset[Gbl][B];
    kl = C2->col_offset[Gca][Gk];

    ncols = occpi[Gk];
    nrows = E->params->coltot[Gij];
    nlinks = occpi[Gl];

    if (nrows && ncols && nlinks) {
      E->matrix[Gbl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(E, Gbl, bl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, E->matrix[Gbl][0], nrows,
          &(C2->matrix[Gca][ca][kl]), nlinks, 1.0, W[Gij][0], ncols);
      global_dpd_->free_dpd_block(E->matrix[Gbl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(C2, Gca);

    /* -t_klcb * E_jila, rewritten as t_cbkl * E_alij */
    Gal = Ga ^ Gl;
    Gk = Gcb ^ Gl ^ GC;
    Gij = Gal ^ GE;

    global_dpd_->buf4_mat_irrep_init(C2, Gcb);
    global_dpd_->buf4_mat_irrep_rd(C2, Gcb);

    al = E->row_offset[Gal][A];
    kl = C2->col_offset[Gcb][Gk];

    ncols = occpi[Gk];
    nrows = E->params->coltot[Gij];
    nlinks = occpi[Gl];

    if (nrows && ncols && nlinks) {
      E->matrix[Gal] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(E, Gal, al, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, E->matrix[Gal][0], nrows,
          &(C2->matrix[Gcb][cb][kl]), nlinks, 1.0, W[Gij][0], ncols);
      global_dpd_->free_dpd_block(E->matrix[Gal], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(C2, Gcb);
  } /* Gl */

  /*=======================================*/
  /* The terms that contribute to W[ik][j] */
  /*=======================================*/
  for (Gd = 0; Gd < nirreps; Gd++) {
    /* +t_ikcd * F_jdab, rewritten as +t_cdik * F_abjd */
    Gcd = Gc ^ Gd;
    Gj = Gd ^ Gab ^ GF;
    Gik = Gcd ^ GC;

    cd = C2->row_offset[Gcd][C];
    jd = F->col_offset[Gab][Gj];

    ncols = occpi[Gj];
    nrows = C2->params->coltot[Gik];
    nlinks = virtpi[Gd];

    if (nrows && ncols && nlinks) {
      C2->matrix[Gcd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(C2, Gcd, cd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, C2->matrix[Gcd][0], nrows,
          &(F->matrix[Gab][ab][jd]), nlinks, 1.0, W2[Gik][0], ncols);
      global_dpd_->free_dpd_block(C2->matrix[Gcd], nlinks, nrows);
    }

    /* +t_ikbd * F_jdca, rewritten as +t_bdik * F_cajd */
    Gbd = Gb ^ Gd;
    Gj = Gd ^ Gca ^ GF;
    Gik = Gbd ^ GC;

    bd = C2->row_offset[Gbd][B];
    jd = F->col_offset[Gca][Gj];

    ncols = occpi[Gj];
    nrows = C2->params->coltot[Gik];
    nlinks = virtpi[Gd];

    if (nrows && ncols && nlinks) {
      C2->matrix[Gbd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(C2, Gbd, bd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, C2->matrix[Gbd][0], nrows,
          &(F->matrix[Gca][ca][jd]), nlinks, 1.0, W2[Gik][0], ncols);
      global_dpd_->free_dpd_block(C2->matrix[Gbd], nlinks, nrows);
    }

    /* -t_ikad * F_jdcb, rewritten as t_adik * F_bcjd */
    Gad = Ga ^ Gd;
    Gj = Gd ^ Gbc ^ GF;
    Gik = Gad ^ GC;

    ad = C2->row_offset[Gad][A];
    jd = F->col_offset[Gbc][Gj];

    ncols = occpi[Gj];
    nrows = C2->params->coltot[Gik];
    nlinks = virtpi[Gd];

    if (nrows && ncols && nlinks) {
      C2->matrix[Gad] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(C2, Gad, ad, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, C2->matrix[Gad][0], nrows,
          &(F->matrix[Gbc][bc][jd]), nlinks, 1.0, W2[Gik][0], ncols);
      global_dpd_->free_dpd_block(C2->matrix[Gad], nlinks, nrows);
    }
  } /* Gd */

  for (Gl = 0; Gl < nirreps; ++Gl) {
    /* +t_jlab * E_iklc, rewritten as +t_abjl * E_clik */
    Gcl = Gc ^ Gl;
    Gj = Gab ^ Gl ^ GC;
    Gik = Gcl ^ GE;

    global_dpd_->buf4_mat_irrep_init(C2, Gab);
    global_dpd_->buf4_mat_irrep_rd(C2, Gab);

    cl = E->row_offset[Gcl][C];
    jl = C2->col_offset[Gab][Gj];

    ncols = occpi[Gj];
    nrows = E->params->coltot[Gik];
    nlinks = occpi[Gl];

    if (nrows && ncols && nlinks) {
      E->matrix[Gcl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(E, Gcl, cl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, E->matrix[Gcl][0], nrows,
          &(C2->matrix[Gab][ab][jl]), nlinks, 1.0, W2[Gik][0], ncols);
      global_dpd_->free_dpd_block(E->matrix[Gcl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(C2, Gab);

    /* +t_jlca * E_iklb, rewritten as +t_cajl * E_blik */
    Gbl = Gb ^ Gl;
    Gj = Gca ^ Gl ^ GC;
    Gik = Gbl ^ GE;

    global_dpd_->buf4_mat_irrep_init(C2, Gca);
    global_dpd_->buf4_mat_irrep_rd(C2, Gca);

    bl = E->row_offset[Gbl][B];
    jl = C2->col_offset[Gca][Gj];

    ncols = occpi[Gj];
    nrows = E->params->coltot[Gik];
    nlinks = occpi[Gl];

    if (nrows && ncols && nlinks) {
      E->matrix[Gbl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(E, Gbl, bl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, E->matrix[Gbl][0], nrows,
          &(C2->matrix[Gca][ca][jl]), nlinks, 1.0, W2[Gik][0], ncols);
      global_dpd_->free_dpd_block(E->matrix[Gbl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(C2, Gca);

    /* -t_jlcb * E_ikla, rewritten as -t_cbjl * E_alik */
    Gal = Ga ^ Gl;
    Gj = Gcb ^ Gl ^ GC;
    Gik = Gal ^ GE;

    global_dpd_->buf4_mat_irrep_init(C2, Gcb);
    global_dpd_->buf4_mat_irrep_rd(C2, Gcb);

    al = E->row_offset[Gal][A];
    jl = C2->col_offset[Gcb][Gj];

    ncols = occpi[Gj];
    nrows = E->params->coltot[Gik];
    nlinks = occpi[Gl];

    if (nrows && ncols && nlinks) {
      E->matrix[Gal] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(E, Gal, al, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, E->matrix[Gal][0], nrows,
          &(C2->matrix[Gcb][cb][jl]), nlinks, 1.0, W2[Gik][0], ncols);
      global_dpd_->free_dpd_block(E->matrix[Gal], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(C2, Gcb);
  } /* Gl */

  //Sort W2[ik][j]->W[ij][k]
  global_dpd_->sort_3d(W2, W, nirreps, Gabc ^ GX3, C2->params->coltot,
      C2->params->colidx, C2->params->colorb, C2->params->rsym,
      C2->params->ssym, occ_off, occ_off, occpi, occ_off, C2->params->colidx,
      acb, 1);

  //Reset the W2 buffer for the next set of contractions
  for (Gij = 0; Gij < nirreps; ++Gij) {
    Gk = Gij ^ Gabc ^ GX3; /* changed */
    if (C2->params->coltot[Gij] && occpi[Gk]) {
      memset(W2[Gij][0], 0, C2->params->coltot[Gij] * occpi[Gk]
          * sizeof(double));
    }
  } /* Gij */

  /*=======================================*/
  /* The terms that contribute to W[jk][i] */
  /*=======================================*/
  for (Gd = 0; Gd < nirreps; Gd++) {
    /* +t_kjcd * F_idab, rewritten as -t_cdjk * F_abid */
    Gcd = Gc ^ Gd;
    Gi = Gd ^ Gab ^ GF;
    Gjk = Gcd ^ GC;

    cd = C2->row_offset[Gcd][C];
    id = F->col_offset[Gab][Gi];

    ncols = occpi[Gi];
    nrows = C2->params->coltot[Gjk];
    nlinks = virtpi[Gd];

    if (nrows && ncols && nlinks) {
      C2->matrix[Gcd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(C2, Gcd, cd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, C2->matrix[Gcd][0], nrows,
          &(F->matrix[Gab][ab][id]), nlinks, 1.0, W2[Gjk][0], ncols);
      global_dpd_->free_dpd_block(C2->matrix[Gcd], nlinks, nrows);
    }

    /* +t_kjbd * F_idca, rewritten as -t_bdjk * F_caid */
    Gbd = Gb ^ Gd;
    Gi = Gd ^ Gca ^ GF;
    Gjk = Gbd ^ GC;

    bd = C2->row_offset[Gbd][B];
    id = F->col_offset[Gca][Gi];

    ncols = occpi[Gi];
    nrows = C2->params->coltot[Gjk];
    nlinks = virtpi[Gd];

    if (nrows && ncols && nlinks) {
      C2->matrix[Gbd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(C2, Gbd, bd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, C2->matrix[Gbd][0], nrows,
          &(F->matrix[Gca][ca][id]), nlinks, 1.0, W2[Gjk][0], ncols);
      global_dpd_->free_dpd_block(C2->matrix[Gbd], nlinks, nrows);
    }

    /* -t_kjad * F_idcb, rewritten as t_adjk * F_cbid */
    Gad = Ga ^ Gd;
    Gi = Gd ^ Gcb ^ GF;
    Gjk = Gad ^ GC;

    ad = C2->row_offset[Gad][A];
    id = F->col_offset[Gcb][Gi];

    ncols = occpi[Gi];
    nrows = C2->params->coltot[Gjk];
    nlinks = virtpi[Gd];

    if (nrows && ncols && nlinks) {
      C2->matrix[Gad] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(C2, Gad, ad, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, C2->matrix[Gad][0], nrows,
          &(F->matrix[Gcb][cb][id]), nlinks, 1.0, W2[Gjk][0], ncols);
      global_dpd_->free_dpd_block(C2->matrix[Gad], nlinks, nrows);
    }
  } /* Gd */

  for (Gl = 0; Gl < nirreps; Gl++) {
    /* -t_ilab * E_jklc, rewritten as -t_abil * E_cljk */
    Gcl = Gc ^ Gl;
    Gi = Gab ^ Gl ^ GC;
    Gjk = Gcl ^ GE;

    global_dpd_->buf4_mat_irrep_init(C2, Gab);
    global_dpd_->buf4_mat_irrep_rd(C2, Gab);

    cl = E->row_offset[Gcl][C];
    il = C2->col_offset[Gab][Gi];

    ncols = occpi[Gi];
    nrows = E->params->coltot[Gjk];
    nlinks = occpi[Gl];

    if (nrows && ncols && nlinks) {
      E->matrix[Gcl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(E, Gcl, cl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, E->matrix[Gcl][0], nrows,
          &(C2->matrix[Gab][ab][il]), nlinks, 1.0, W2[Gjk][0], ncols);
      global_dpd_->free_dpd_block(E->matrix[Gcl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(C2, Gab);

    /* -t_ilca * E_jklb, rewritten as -t_cail * E_bljk */
    Gbl = Gb ^ Gl;
    Gi = Gca ^ Gl ^ GC;
    Gjk = Gbl ^ GE;

    global_dpd_->buf4_mat_irrep_init(C2, Gca);
    global_dpd_->buf4_mat_irrep_rd(C2, Gca);

    bl = E->row_offset[Gbl][B];
    il = C2->col_offset[Gca][Gi];

    ncols = occpi[Gi];
    nrows = E->params->coltot[Gjk];
    nlinks = occpi[Gl];

    if (nrows && ncols && nlinks) {
      E->matrix[Gbl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(E, Gbl, bl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, E->matrix[Gbl][0], nrows,
          &(C2->matrix[Gca][ca][il]), nlinks, 1.0, W2[Gjk][0], ncols);
      global_dpd_->free_dpd_block(E->matrix[Gbl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(C2, Gca);

    /* +t_ilcb * E_jkla, rewritten as t_cbil * E_aljk */
    Gal = Ga ^ Gl;
    Gi = Gcb ^ Gl ^ GC;
    Gjk = Gal ^ GE;

    global_dpd_->buf4_mat_irrep_init(C2, Gcb);
    global_dpd_->buf4_mat_irrep_rd(C2, Gcb);

    al = E->row_offset[Gal][A];
    il = C2->col_offset[Gcb][Gi];

    ncols = occpi[Gi];
    nrows = E->params->coltot[Gjk];
    nlinks = occpi[Gl];

    if (nrows && ncols && nlinks) {
      E->matrix[Gal] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(E, Gal, al, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, E->matrix[Gal][0], nrows,
          &(C2->matrix[Gcb][cb][il]), nlinks, 1.0, W2[Gjk][0], ncols);
      global_dpd_->free_dpd_block(E->matrix[Gal], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(C2, Gcb);
  } /* Gl */
  //Sort W2[jk][i]->W[ij][k]
  global_dpd_->sort_3d(W2, W, nirreps, Gabc ^ GX3, C2->params->coltot,
      C2->params->colidx, C2->params->colorb, C2->params->rsym,
      C2->params->ssym, occ_off, occ_off, occpi, occ_off, C2->params->colidx,
      cab, 1);
  /* The disconnected diagrams */
  if (disc) {
    for (h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(D, h);
      global_dpd_->buf4_mat_irrep_rd(D, h);
      global_dpd_->buf4_mat_irrep_init(C2, h);
      global_dpd_->buf4_mat_irrep_rd(C2, h);
    }

    for (Gij=0; Gij < nirreps; ++Gij) {
      Gk = Gij ^ Gabc;
      Gji = Gij;

      for (ij=0; ij < C2->params->coltot[Gij]; ++ij) {
        I = C2->params->colorb[Gij][ij][0];
        Gi = C2->params->rsym[I];
        i = I - occ_off[Gi];
        J = C2->params->colorb[Gij][ij][1];
        Gj = C2->params->ssym[J];
        j = J - occ_off[Gj];

        Gjk = Gkj = Gj ^ Gk;
        Gik = Gki = Gi ^ Gk;

        for (k=0; k < occpi[Gk]; ++k) {
          K = occ_off[Gk] + k;
          ij = D->params->rowidx[I][J];
          ik = D->params->rowidx[I][K];
          jk = D->params->rowidx[J][K];
          ji = D->params->rowidx[J][I];

          /* +t_ia * D_jkbc + f_ia * t_bcjk */
          if (Gi == Ga && Gjk == Gbc) {
            t_ia = D_jkbc = f_ia = t_jkbc = 0.0;

            if (C1->params->rowtot[Ga] && C1->params->coltot[Ga]) {
              t_ia = C1->matrix[Ga][i][a];
              f_ia = fIA->matrix[Ga][i][a];
            }

            if (D->params->rowtot[Gbc] && D->params->coltot[Gbc]) {
              D_jkbc = D->matrix[Gbc][jk][bc];
              t_jkbc = C2->matrix[Gbc][bc][jk];
            }
            V[Gij][ij][k] += t_ia * D_jkbc + f_ia * t_jkbc;
          }

          /* -t_ib * D_jkac - f_ib * t_acjk */
          if (Gi == Gb && Gjk == Gac) {
            t_ib = D_jkac = f_ib = t_jkac = 0.0;

            if (C1->params->rowtot[Gb] && C1->params->coltot[Gb]) {
              t_ib = C1->matrix[Gb][i][b];
              f_ib = fIA->matrix[Gb][i][b];
            }

            if (D->params->rowtot[Gac] && D->params->coltot[Gac]) {
              D_jkac = D->matrix[Gac][jk][ac];
              t_jkac = C2->matrix[Gac][ac][jk];
            }
            V[Gij][ij][k] -= t_ib * D_jkac + f_ib * t_jkac;
          }

          /* +t_ic * D_jkab + f_ic * t_abjk */
          if (Gi == Gc && Gjk == Gab) {
            t_ic = D_jkba = f_ic = t_jkba = 0.0;

            if (C1->params->rowtot[Gc] && C1->params->coltot[Gc]) {
              t_ic = C1->matrix[Gc][i][c];
              f_ic = fIA->matrix[Gc][i][c];
            }

            if (D->params->rowtot[Gab] && D->params->coltot[Gab]) {
              D_jkba = D->matrix[Gba][jk][ab];
              t_jkba = C2->matrix[Gab][ab][jk];
            }

            V[Gij][ij][k] += t_ic * D_jkba + f_ic * t_jkba;
          }

          /* -t_ja * D_ikbc - f_ja * t_bcik*/
          if (Gj == Ga && Gik == Gbc) {
            t_ja = D_ikbc = f_ja = t_ikbc = 0.0;

            if (C1->params->rowtot[Ga] && C1->params->coltot[Ga]) {
              t_ja = C1->matrix[Ga][j][a];
              f_ja = fIA->matrix[Ga][j][a];
            }

            if (D->params->rowtot[Gbc] && D->params->coltot[Gbc]) {
              D_ikbc = D->matrix[Gbc][ik][bc];
              t_ikbc = C2->matrix[Gbc][bc][ik];
            }

            V[Gij][ij][k] -= t_ja * D_ikbc + f_ja * t_ikbc;
          }

          /* +t_jb * D_ikac + f_jb * t_acik */
          if (Gj == Gb && Gik == Gac) {
            t_jb = D_ikac = f_jb = t_ikac = 0.0;

            if (C1->params->rowtot[Gb] && C1->params->coltot[Gb]) {
              t_jb = C1->matrix[Gb][j][b];
              f_jb = fIA->matrix[Gb][j][b];
            }

            if (D->params->rowtot[Gac] && D->params->coltot[Gac]) {
              D_ikac = D->matrix[Gac][ik][ac];
              t_ikac = C2->matrix[Gac][ac][ik];
            }

            V[Gij][ij][k] += t_jb * D_ikac + f_jb * t_ikac;
          }

          /* t_jc * D_ikba  f_jc * t_baik */
          if (Gj == Gc && Gik == Gba) {
            t_jc = D_ikba = f_jc = t_ikba = 0.0;

            if (C1->params->rowtot[Gc] && C1->params->coltot[Gc]) {
              t_jc = C1->matrix[Gc][j][c];
              f_jc = fIA->matrix[Gc][j][c];
            }

            if (D->params->rowtot[Gba] && D->params->coltot[Gba]) {
              D_ikba = D->matrix[Gba][ik][ba];
              t_ikba = C2->matrix[Gba][ba][ik];
            }

            V[Gij][ij][k] += t_jc * D_ikba + f_jc * t_ikba;
          }

          /* -t_ka * D_jibc - f_ka * t_jibc */
          if (Gk == Ga && Gji == Gbc) {
            t_ka = D_jibc = f_ka = t_jibc = 0.0;

            if (C1->params->rowtot[Ga] && C1->params->coltot[Ga]) {
              t_ka = C1->matrix[Ga][k][a];
              f_ka = fIA->matrix[Ga][k][a];
            }

            if (D->params->rowtot[Gbc] && D->params->coltot[Gbc]) {
              D_jibc = D->matrix[Gbc][ji][bc];
              t_jibc = C2->matrix[Gbc][bc][ji];
            }

            V[Gij][ij][k] -= t_ka * D_jibc + f_ka * t_jibc;
          }

          /* +t_kb * D_jiac + f_kb * t_jiac */
          if (Gk == Gb && Gji == Gac) {
            t_kb = D_jiac = f_kb = t_jiac = 0.0;

            if (C1->params->rowtot[Gb] && C1->params->coltot[Gb]) {
              t_kb = C1->matrix[Gb][k][b];
              f_kb = fIA->matrix[Gb][k][b];
            }

            if (D->params->rowtot[Gac] && D->params->coltot[Gac]) {
              D_jiac = D->matrix[Gac][ji][ac];
              t_jiac = C2->matrix[Gac][ac][ji];
            }

            V[Gij][ij][k] += t_kb * D_jiac + f_kb * t_jiac;
          }

          /* t_kc * D_jiba + f_kc * t_jiba*/
          if (Gk == Gc && Gji == Gba) {
            t_kc = D_jiba = f_kc = t_jiba = 0.0;

            if (C1->params->rowtot[Gc] && C1->params->coltot[Gc]) {
              t_kc = C1->matrix[Gc][k][c];
              f_kc = fIA->matrix[Gc][k][c];
            }

            if (D->params->rowtot[Gba] && D->params->coltot[Gba]) {
              D_jiba = D->matrix[Gba][ji][ba];
              t_jiba = C2->matrix[Gba][ba][ji];
            }

            V[Gij][ij][k] += t_kc * D_jiba + f_kc * t_jiba;
          }
        } /* k */
      } /* ij */
    } /* Gij */

    for (h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_close(C2, h);
      global_dpd_->buf4_mat_irrep_close(D, h);
    }
  } /**** Disconnected T3 complete ****/

  for (Gij=0; Gij < nirreps; ++Gij) {
    Gk = Gij ^ Gabc ^ GX3; /* assumes totally symmetric! */
    for (ij=0; ij < C2->params->coltot[Gij]; ++ij) {
      I = C2->params->colorb[Gij][ij][0];
      Gi = C2->params->rsym[I];
      i = I - occ_off[Gi];
      J = C2->params->colorb[Gij][ij][1];
      Gj = C2->params->ssym[J];
      j = J - occ_off[Gj];

      for (k=0; k < occpi[Gk]; ++k) {
        denom = dabc;
        if (fIJ->params->rowtot[Gi])
          denom += fIJ->matrix[Gi][i][i];
        if (fIJ->params->rowtot[Gj])
          denom += fIJ->matrix[Gj][j][j];
        if (fIJ->params->rowtot[Gk])
          denom += fIJ->matrix[Gk][k][k];
        W[Gij][ij][k] /= (omega + denom);
        if (disc)
          V[Gij][ij][k] /= (omega + denom);
      } /* k */
    } /* ij */
  } /* Gij */

  for (Gij = 0; Gij < nirreps; ++Gij) {
    Gk = Gij ^ Gabc ^ GX3; /* changed */
    global_dpd_->free_dpd_block(W2[Gij], C2->params->coltot[Gij], occpi[Gk]);
  } /* Gij */

  free(W2);

  global_dpd_->file2_mat_close(fIJ);
  global_dpd_->file2_mat_close(fAB);
  if (disc) {
    global_dpd_->file2_mat_close(fIA);
    global_dpd_->file2_mat_close(C1);
  }
  for (h = 0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_close(F, h);
  }

}

/* T3_UHF_AAB_abc(): Computes all T3(IJk,ABc) amplitudes for a given A, B,
 ** C combination for input T2, F, and E intermediates.  This function
 ** will work for AAB or BBA spin cases with either RHF/ROHF or UHF
 ** orbitals.  Based on the T3_UHF_AAB() routine.
 **
 ** Arguments:
 **
 **   double ***W: The target triples amplitudes in an nirreps x IJ x k
 **   array.  The memory for this must be allocated externally.
 **
 **   double ***V: The target disconnected triples amplitudes (if
 **   requested) in an nirreps x IJ x k array.  The memory for this
 **   must be allocated externally.
 **
 **   int disc: Boolean: 1 == computed disconnected triples; 0 == don't
 **
 **   int nirreps: Number of irreps.
 **
 **   int A: Absolute index of orbital A.
 **
 **   int Ga: Irrep of A.
 **
 **   int B: Absolute index of orbital B.
 **
 **   int Gb: Irrep of B.
 **
 **   int C: Absolute index of orbital C.
 **
 **   int Gc: Irrep of C.
 **
 **   dpdbuf4 *T2: Pointer to dpd buffer for double excitation amps,
 **   ordered (AB,IJ).
 **
 **   dpdbuf4 *F: Pointer to dpd buffer for three-virtual-index
 **   intermediate, ordered (BC,IA).
 **
 **   dpdbuf4 *E: Pointer to dpd buffer for three-occupied-index
 **   intermediate, ordered (AK,IJ).
 **
 **   dpdfile *C1: If disconnected T3's are requested, pointer to dpd
 **   buffer for single-excitation amps.
 **
 **   dpdbuf4 *D: If disconnected T3's are requested, pointer to dpd
 **   buffer for <IJ||ab> integrals.
 **
 **   dpdfile2 *fIA: Pointer to the dpd file2 for the occ-vir block of
 **   the Fock matrix (or other appropriate one-electron operator).
 **
 **   dpdfile2 *fIJ: Pointer to the dpd file2 for the occ-occ block of
 **   the Fock matrix (or other appropriate one-electron operator).
 **
 **   dpdfile2 *fAB: Pointer to the dpd file2 for the vir-vir block of
 **   the Fock matrix (or other appropriate one-electron operator).
 **
 **   int *occpi: Number of occupied orbitals per irrep lookup array.
 **
 **   int *occ_off: Offset lookup for translating between absolute and
 **   relative orbital indices for occupied space.
 **
 **   int *virtpi: Number of virtual orbitals per irrep lookup array.
 **
 **   int *vir_off: Offset lookup for translating between absolute and
 **   relative orbital indices for virtual space.
 **
 **   double omega: constant to add to denominators - needed for
 **   CC3 EOM
 **
 **   N.B. The ordering of the F and E integrals, and the T2 amplitude
 **   is different from that used in the original T3_UHF_AAB() function.
 ** ACS, December 2008
 */

void T3_UHF_AAB_abc(double ***W, double ***V, int disc, int nirreps, int A,
    int Ga, int B, int Gb, int C, int Gc, dpdbuf4 *T2AA, dpdbuf4 *T2AB,
    dpdbuf4 *T2BA, dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA, dpdbuf4 *EAA,
    dpdbuf4 *EAB, dpdbuf4 *EBA, dpdfile2 *T1A, dpdfile2 *T1B, dpdbuf4 *DAA,
    dpdbuf4 *DAB, dpdfile2 *fIA, dpdfile2 *fia, dpdfile2 *fIJ, dpdfile2 *fij,
    dpdfile2 *fAB, dpdfile2 *fab, int *aoccpi, int *aocc_off, int *boccpi,
    int *bocc_off, int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off,
    double omega) {
  int h;
  int a, b, c;
  int ca, cb;
  int ij, ji, ik, ki, jk, kj;
  int Gij, Gji, Gik, Gki, Gjk, Gkj, Gabc;
  int Gi, Gj, Gk;
  int Gd, Gl;
  int Gad, Gbd, Gcd;
  int Gab, Gcb, Gca, Gac, Gbc, Gba;
  int Gla, Glb, Glc;
  int Gal, Gbl, Gcl;
  int i, j, k, I, J, K;
  int ab, bc, ac;
  int dc, bd, ad, cd;
  int id, jd, kd;
  int il, jl, kl;
  int la, lb, lc;
  int al, bl, cl;
  int nrows, ncols, nlinks;
  double dabc, denom;
  double ***W2;
  int GE, GF, GC, GX3;
  double t_ia, t_ib, t_ja, t_jb, t_kc;
  double f_ia, f_ib, f_ja, f_jb, f_kc;
  double D_jkbc, D_jkac, D_ikbc, D_ikac, D_jiab;
  double t_jkbc, t_jkac, t_ikbc, t_ikac, t_jiab;

  GC = T2AA->file.my_irrep;
  /* F and E are assumed to have same irrep */
  GF = GE = FAA->file.my_irrep;
  GX3 = GC^GF;

  global_dpd_->file2_mat_init(fIJ);
  global_dpd_->file2_mat_init(fAB);
  global_dpd_->file2_mat_init(fij);
  global_dpd_->file2_mat_init(fab);
  if (disc) {
    global_dpd_->file2_mat_init(fIA);
    global_dpd_->file2_mat_init(fia);
    global_dpd_->file2_mat_init(T1A);
    global_dpd_->file2_mat_init(T1B);
  }

  global_dpd_->file2_mat_rd(fIJ);
  global_dpd_->file2_mat_rd(fAB);
  global_dpd_->file2_mat_rd(fij);
  global_dpd_->file2_mat_rd(fab);
  if (disc) {
    global_dpd_->file2_mat_rd(fIA);
    global_dpd_->file2_mat_rd(fia);
    global_dpd_->file2_mat_rd(T1A);
    global_dpd_->file2_mat_rd(T1B);
  }

  for (h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(FAA, h);
    global_dpd_->buf4_mat_irrep_rd(FAA, h);
    global_dpd_->buf4_mat_irrep_init(FAB, h);
    global_dpd_->buf4_mat_irrep_rd(FAB, h);
    global_dpd_->buf4_mat_irrep_init(FBA, h);
    global_dpd_->buf4_mat_irrep_rd(FBA, h);
  }

  a = A - avir_off[Ga];
  b = B - avir_off[Gb];
  c = C - bvir_off[Gc];

  Gab = Gba = Ga ^ Gb;
  Gac = Gca = Ga ^ Gc;
  Gbc = Gcb = Gb ^ Gc;
  Gabc = Ga ^ Gb ^ Gc;

  ab = T2AA->params->rowidx[A][B];
  cb = T2BA->params->rowidx[C][B];
  ca = T2BA->params->rowidx[C][A];
  ac = T2AB->params->rowidx[A][C];
  bc = T2AB->params->rowidx[B][C];

  dabc = 0.0;
  if (fAB->params->rowtot[Ga])
    dabc -= fAB->matrix[Ga][a][a];
  if (fAB->params->rowtot[Gb])
    dabc -= fAB->matrix[Gb][b][b];
  if (fab->params->rowtot[Gc])
    dabc -= fab->matrix[Gc][c][c];

  W2 = (double ***) malloc(nirreps * sizeof(double **)); /* alpha-beta-alpha */

  /* clear out the old W */
  for (Gij=0; Gij < nirreps; Gij++) {
    Gk = Gabc ^ Gij ^ GX3; /* assumes totally symmetric! */

    if (T2AA->params->coltot[Gij] && boccpi[Gk]) {
      memset(W[Gij][0], 0, T2AA->params->coltot[Gij]*boccpi[Gk]*sizeof(double));
      memset(V[Gij][0], 0, T2AA->params->coltot[Gij]*boccpi[Gk]*sizeof(double));
    }
  }

  for (Gd=0; Gd < nirreps; Gd++) {
    /* -t_JIAD * F_kDcB, rewritten as t_ADIJ * F_cBkD */
    Gad = Ga ^ Gd;
    Gk = Gd ^ Gcb ^ GF;
    Gij = Gad ^ GC;

    kd = FBA->col_offset[Gcb][Gk];
    ad = T2AA->row_offset[Gad][A];

    ncols = boccpi[Gk];
    nrows = T2AA->params->coltot[Gij];
    nlinks = avirtpi[Gd];

    if (nrows && ncols && nlinks) {
      T2AA->matrix[Gad] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(T2AA, Gad, ad, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, T2AA->matrix[Gad][0], nrows,
          &(FBA->matrix[Gcb][cb][kd]), nlinks, 1.0, W[Gij][0], ncols);
      global_dpd_->free_dpd_block(T2AA->matrix[Gad], nlinks, nrows);
    }

    /* +t_JIBD * F_kDcA, rewritten as -t_BDIJ * F_cAkD */
    Gbd = Gb ^ Gd;
    Gk = Gd ^ Gca ^ GF;
    Gij = Gbd ^ GC;

    kd = FBA->col_offset[Gca][Gk];
    bd = T2AA->row_offset[Gbd][B];

    ncols = boccpi[Gk];
    nrows = T2AA->params->coltot[Gij];
    nlinks = avirtpi[Gd];

    if (nrows && ncols && nlinks) {
      T2AA->matrix[Gbd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(T2AA, Gbd, bd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, T2AA->matrix[Gbd][0],
          nrows, &(FBA->matrix[Gca][ca][kd]), nlinks, 1.0, W[Gij][0], ncols);
      global_dpd_->free_dpd_block(T2AA->matrix[Gbd], nlinks, nrows);
    }
  } /* Gd */

  for (Gl=0; Gl < nirreps; ++Gl) {
    /* -t_kLcB * E_JILA, rewritten as t_cBkL * E_ALIJ */
    Gal = Ga ^ Gl;
    Gk = Gcb ^ Gl ^ GC;
    Gij = Gal ^ GE;

    global_dpd_->buf4_mat_irrep_init(T2BA, Gcb);
    global_dpd_->buf4_mat_irrep_rd(T2BA, Gcb);

    al = EAA->row_offset[Gal][A];
    kl = T2BA->col_offset[Gcb][Gk];

    ncols = boccpi[Gk];
    nrows = EAA->params->coltot[Gij];
    nlinks = aoccpi[Gl];

    if (nrows && ncols && nlinks) {
      EAA->matrix[Gal] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(EAA, Gal, al, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, EAA->matrix[Gal][0], nrows,
          &(T2BA->matrix[Gcb][cb][kl]), nlinks, 1.0, W[Gij][0], ncols);
      global_dpd_->free_dpd_block(EAA->matrix[Gal], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(T2BA, Gcb);

    /* +t_kLcA * E_JILB, rewritten as -t_cAkL * E_BLIJ */
    Gbl = Gb ^ Gl;
    Gk = Gca ^ Gl ^ GC;
    Gij = Gbl ^ GE;

    global_dpd_->buf4_mat_irrep_init(T2BA, Gca);
    global_dpd_->buf4_mat_irrep_rd(T2BA, Gca);

    bl = EAA->row_offset[Gbl][B];
    kl = T2BA->col_offset[Gca][Gk];

    ncols = boccpi[Gk];
    nrows = EAA->params->coltot[Gij];
    nlinks = aoccpi[Gl];

    if (nrows && ncols && nlinks) {
      EAA->matrix[Gbl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(EAA, Gbl, bl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, EAA->matrix[Gbl][0], nrows,
          &(T2BA->matrix[Gca][ca][kl]), nlinks, 1.0, W[Gij][0], ncols);
      global_dpd_->free_dpd_block(EAA->matrix[Gbl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(T2BA, Gca);
  } /* Gl */

  /* Allocate the alpha-beta-alpha array */
  for (Gij=0; Gij < nirreps; ++Gij) {
    Gk = Gabc ^ Gij ^GX3; /* assumes totally symmetric! */
    W2[Gij] = global_dpd_->dpd_block_matrix(T2AB->params->coltot[Gij], aoccpi[Gk]);
  }

  for (Gd=0; Gd < nirreps; ++Gd) {
    /* -t_IkBd * F_JdAc, rewritten as -tBdIk * F_AcJd */
    Gbd = Gb ^ Gd;
    Gj = Gd ^ Gac ^ GF;
    Gik = Gbd ^ GC;

    jd = FAB->col_offset[Gac][Gj];
    bd = T2AB->row_offset[Gbd][B];

    ncols = aoccpi[Gj];
    nrows = T2AB->params->coltot[Gik];
    nlinks = bvirtpi[Gd];

    if (nrows && ncols && nlinks) {
      T2AB->matrix[Gbd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(T2AB, Gbd, bd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, T2AB->matrix[Gbd][0],
          nrows, &(FAB->matrix[Gac][ac][jd]), nlinks, 1.0, W2[Gik][0], ncols);
      global_dpd_->free_dpd_block(T2AB->matrix[Gbd], nlinks, nrows);
    }

    /* +t_IkAd * F_JdBc, rewritten as t_AdIk * F_BcJd */
    Gad = Ga ^ Gd;
    Gj = Gd ^ Gbc ^ GF;
    Gik = Gad ^ GC;

    jd = FAB->col_offset[Gbc][Gj];
    ad = T2AB->row_offset[Gad][A];

    ncols = aoccpi[Gj];
    nrows = T2AB->params->coltot[Gik];
    nlinks = bvirtpi[Gd];

    if (nrows && ncols && nlinks) {
      T2AB->matrix[Gad] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(T2AB, Gad, ad, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, T2AB->matrix[Gad][0], nrows,
          &(FAB->matrix[Gbc][bc][jd]), nlinks, 1.0, W2[Gik][0], ncols);
      global_dpd_->free_dpd_block(T2AB->matrix[Gad], nlinks, nrows);
    }
  } /* Gd */

  for (Gl=0; Gl < nirreps; ++Gl) {
    /* +t_JLAB * E_IkLc, rewritten as t_ABJL * E_cLIk */
    Gcl = Gc ^ Gl;
    Gj = Gab ^ Gl ^ GC;
    Gik = Gcl ^ GE;

    global_dpd_->buf4_mat_irrep_init(T2AA, Gab);
    global_dpd_->buf4_mat_irrep_rd(T2AA, Gab);

    cl = EAB->row_offset[Gcl][C];
    jl = T2AA->col_offset[Gab][Gj];

    ncols = aoccpi[Gj];
    nrows = EAB->params->coltot[Gik];
    nlinks = aoccpi[Gl];

    if (nrows && ncols && nlinks) {
      EAB->matrix[Gcl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(EAB, Gcl, cl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, EAB->matrix[Gcl][0], nrows,
          &(T2AA->matrix[Gab][ab][jl]), nlinks, 1.0, W2[Gik][0], ncols);
      global_dpd_->free_dpd_block(EAB->matrix[Gcl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(T2AA, Gab);
  } /* Gl */

  /*=== Sort W2[Ik][J] -> W[IJ][k] ===*/
  global_dpd_->sort_3d(W2, W, nirreps, Gabc^GX3, T2AB->params->coltot,
      T2AB->params->colidx, T2AB->params->colorb, T2AB->params->rsym,
      T2AB->params->ssym, aocc_off, bocc_off, aoccpi, aocc_off,
      T2AA->params->colidx, acb, 1);

  /* Clear the W2 array */
  for (Gij=0; Gij < nirreps; ++Gij) {
    Gk = Gabc ^ Gij ^ GX3;
    if (T2BA->params->coltot[Gij] && aoccpi[Gk]) {
      memset(W2[Gij][0], 0, T2BA->params->coltot[Gij]*aoccpi[Gk]*sizeof(double));
    }
  }

  for (Gl=0; Gl < nirreps; ++Gl) {
    /* -t_IlAc * E_kJlB, rewritten as -t_AcIl * E_BlkJ */
    Gbl = Gb ^ Gl;
    Gi = Gac ^ Gl ^ GC;
    Gkj = Gbl ^ GE;

    global_dpd_->buf4_mat_irrep_init(T2AB, Gac);
    global_dpd_->buf4_mat_irrep_rd(T2AB, Gac);

    bl = EBA->row_offset[Gbl][B];
    il = T2AB->col_offset[Gac][Gi];

    ncols = aoccpi[Gi];
    nrows = EBA->params->coltot[Gkj];
    nlinks = boccpi[Gl];

    if (nrows && ncols && nlinks) {
      EBA->matrix[Gbl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(EBA, Gbl, bl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, EBA->matrix[Gbl][0], nrows,
          &(T2AB->matrix[Gac][ac][il]), nlinks, 1.0, W2[Gkj][0], ncols);
      global_dpd_->free_dpd_block(EBA->matrix[Gbl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(T2AB, Gac);

    /* +t_IlBc * E_kJlA, rewritten as t_BcIl * E_AlkJ */
    Gal = Ga ^ Gl;
    Gi = Gbc ^ Gl ^ GC;
    Gkj = Gal ^ GE;

    global_dpd_->buf4_mat_irrep_init(T2AB, Gbc);
    global_dpd_->buf4_mat_irrep_rd(T2AB, Gbc);

    al = EBA->row_offset[Gal][A];
    il = T2AB->col_offset[Gbc][Gi];

    ncols = aoccpi[Gi];
    nrows = EBA->params->coltot[Gkj];
    nlinks = boccpi[Gl];

    if (nrows && ncols && nlinks) {
      EBA->matrix[Gal] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(EBA, Gal, al, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, EBA->matrix[Gal][0], nrows,
          &(T2AB->matrix[Gbc][bc][il]), nlinks, 1.0, W2[Gkj][0], ncols);
      global_dpd_->free_dpd_block(EBA->matrix[Gal], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(T2AB, Gbc);
  } /* Gl */

  /*=== Sort W2[kJ][I] -> W[IJ][k] ===*/
  global_dpd_->sort_3d(W2, W, nirreps, Gabc^GX3, T2BA->params->coltot,
      T2BA->params->colidx, T2BA->params->colorb, T2BA->params->rsym,
      T2BA->params->ssym, bocc_off, aocc_off, aoccpi, aocc_off,
      T2AA->params->colidx, cba, 1);

  /* Allocate the beta-alpha-alpha array, after closing the alpha-beta-alpha */
  for (Gij=0; Gij < nirreps; ++Gij) {
    Gk = Gabc ^ Gij ^GX3; /* assumes totally symmetric! */
    global_dpd_->free_dpd_block(W2[Gij], T2AB->params->coltot[Gij], aoccpi[Gk]);
    W2[Gij] = global_dpd_->dpd_block_matrix(T2BA->params->coltot[Gij], aoccpi[Gk]);
  }

  for (Gd=0; Gd < nirreps; ++Gd) {
    /* -t_IkDc * F_JDAB, rewritten as -t_cDkI * F_ABJD */
    Gcd = Gc ^ Gd;
    Gj = Gd ^ Gab ^ GF;
    Gki = Gcd ^ GC;

    jd = FAA->col_offset[Gab][Gj];
    cd = T2BA->row_offset[Gcd][C];

    ncols = aoccpi[Gj];
    nrows = T2BA->params->coltot[Gki];
    nlinks = avirtpi[Gd];

    if (nrows && ncols && nlinks) {
      T2BA->matrix[Gcd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(T2BA, Gcd, cd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, T2BA->matrix[Gcd][0],
          nrows, &(FAA->matrix[Gab][ab][jd]), nlinks, 1.0, W2[Gki][0], ncols);
      global_dpd_->free_dpd_block(T2BA->matrix[Gcd], nlinks, nrows);
    }
  } /* Gd */

  for (Gl=0; Gl < nirreps; ++Gl) {
    /* -t_JlBc * E_kIlA, rewritten as -t_BcJl * E_AlkI */
    Gal = Ga ^ Gl;
    Gj = Gbc ^ Gl ^ GC;
    Gki = Gal ^ GE;

    global_dpd_->buf4_mat_irrep_init(T2AB, Gbc);
    global_dpd_->buf4_mat_irrep_rd(T2AB, Gbc);

    al = EBA->row_offset[Gal][A];
    jl = T2AB->col_offset[Gbc][Gj];

    ncols = aoccpi[Gj];
    nrows = EBA->params->coltot[Gki];
    nlinks = boccpi[Gl];

    if (nrows && ncols && nlinks) {
      EBA->matrix[Gal] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(EBA, Gal, al, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, EBA->matrix[Gal][0], nrows,
          &(T2AB->matrix[Gbc][bc][jl]), nlinks, 1.0, W2[Gki][0], ncols);
      global_dpd_->free_dpd_block(EBA->matrix[Gal], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(T2AB, Gbc);

    /* +t_JlAc * E_kIlB, rewritten as t_AcJl * E_BlkI */
    Gbl = Gb ^ Gl;
    Gj = Gac ^ Gl ^ GC;
    Gki = Gbl ^ GE;

    global_dpd_->buf4_mat_irrep_init(T2AB, Gac);
    global_dpd_->buf4_mat_irrep_rd(T2AB, Gac);

    bl = EBA->row_offset[Gbl][B];
    jl = T2AB->col_offset[Gac][Gj];

    ncols = aoccpi[Gj];
    nrows = EBA->params->coltot[Gki];
    nlinks = boccpi[Gl];

    if (nrows && ncols && nlinks) {
      EBA->matrix[Gbl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(EBA, Gbl, bl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, EBA->matrix[Gbl][0], nrows,
          &(T2AB->matrix[Gac][ac][jl]), nlinks, 1.0, W2[Gki][0], ncols);
      global_dpd_->free_dpd_block(EBA->matrix[Gbl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(T2AB, Gac);
  } /* Gl */

  /*=== Sort W2[kI][J] -> W[IJ][k] ===*/
  global_dpd_->sort_3d(W2, W, nirreps, Gabc^GX3, T2BA->params->coltot,
      T2BA->params->colidx, T2BA->params->colorb, T2BA->params->rsym,
      T2BA->params->ssym, bocc_off, aocc_off, aoccpi, aocc_off,
      T2AA->params->colidx, bca, 1);

  /* Clear the W2 array */
  for (Gij=0; Gij < nirreps; ++Gij) {
    Gk = Gabc ^ Gij ^ GX3;
    if (T2BA->params->coltot[Gij] && aoccpi[Gk]) {
      memset(W2[Gij][0], 0, T2BA->params->coltot[Gij]*aoccpi[Gk]*sizeof(double));
    }
  }

  for (Gd=0; Gd < nirreps; ++Gd) {
    /* +t_JkDc * F_IDAB, rewritten as t_cDkJ * FABID */
    Gcd = Gc ^ Gd;
    Gi = Gd ^ Gab ^ GF;
    Gkj = Gcd ^ GC;

    id = FAA->col_offset[Gab][Gi];
    cd = T2BA->row_offset[Gcd][C];

    ncols = aoccpi[Gi];
    nrows = T2BA->params->coltot[Gkj];
    nlinks = avirtpi[Gd];

    if (nrows && ncols && nlinks) {
      T2BA->matrix[Gcd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(T2BA, Gcd, cd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, T2BA->matrix[Gcd][0], nrows,
          &(FAA->matrix[Gab][ab][id]), nlinks, 1.0, W2[Gkj][0], ncols);
      global_dpd_->free_dpd_block(T2BA->matrix[Gcd], nlinks, nrows);
    }
  } /* Gd */

  /*=== Sort W2[kJ][I] -> W[IJ][k] ===*/
  global_dpd_->sort_3d(W2, W, nirreps, Gabc^GX3, T2BA->params->coltot,
      T2BA->params->colidx, T2BA->params->colorb, T2BA->params->rsym,
      T2BA->params->ssym, bocc_off, aocc_off, aoccpi, aocc_off,
      T2AA->params->colidx, cba, 1);

  /* Clear the W2 array */
  for (Gij=0; Gij < nirreps; ++Gij) {
    Gk = Gabc ^ Gij ^ GX3;
    if (T2AB->params->coltot[Gij] && aoccpi[Gk]) {
      memset(W2[Gij][0], 0, T2AB->params->coltot[Gij]*aoccpi[Gk]*sizeof(double));
    }
  }

  for (Gd=0; Gd < nirreps; ++Gd) {
    /* +t_JkBd * F_IdAc, rewritten as t_BdJk * F_AcId */
    Gbd = Gb ^ Gd;
    Gi = Gd ^ Gac ^ GF;
    Gjk = Gbd ^ GC;

    id = FAB->col_offset[Gac][Gi];
    bd = T2AB->row_offset[Gbd][B];

    ncols = aoccpi[Gi];
    nrows = T2AB->params->coltot[Gjk];
    nlinks = bvirtpi[Gd];

    if (nrows && ncols && nlinks) {
      T2AB->matrix[Gbd] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(T2AB, Gbd, bd, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, T2AB->matrix[Gbd][0], nrows,
          &(FAB->matrix[Gac][ac][id]), nlinks, 1.0, W2[Gjk][0], ncols);
      global_dpd_->free_dpd_block(T2AB->matrix[Gbd], nlinks, nrows);
    }

    /* -t_JkAd * F_IdBc, rewritten as -t_AdJk * F_BcId */
    Gad = Ga ^ Gd;
    Gi = Gd ^ Gbc ^ GF;
    Gjk = Gad ^ GC;

    id = FAB->col_offset[Gbc][Gi];
    ad = T2AB->row_offset[Gad][A];

    ncols = aoccpi[Gi];
    nrows = T2AB->params->coltot[Gjk];
    nlinks = bvirtpi[Gd];

    if (nrows && ncols && nlinks) {
      T2AB->matrix[Gad] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(T2AB, Gad, ad, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, T2AB->matrix[Gad][0],
          nrows, &(FAB->matrix[Gbc][bc][id]), nlinks, 1.0, W2[Gjk][0], ncols);
      global_dpd_->free_dpd_block(T2AB->matrix[Gad], nlinks, nrows);
    }
  } /* Gd */

  for (Gl=0; Gl < nirreps; ++Gl) {
    /* -t_ILAB * E_JkLc, rewritten as -t_ABIL * E_cLJk */
    Gcl = Gc ^ Gl;
    Gi = Gab ^ Gl ^ GC;
    Gjk = Gcl ^ GE;

    global_dpd_->buf4_mat_irrep_init(T2AA, Gab);
    global_dpd_->buf4_mat_irrep_rd(T2AA, Gab);

    cl = EAB->row_offset[Gcl][C];
    il = T2AA->col_offset[Gab][Gi];

    ncols = aoccpi[Gi];
    nrows = EAB->params->coltot[Gjk];
    nlinks = aoccpi[Gl];

    if (nrows && ncols && nlinks) {
      EAB->matrix[Gcl] = global_dpd_->dpd_block_matrix(nlinks, nrows);
      global_dpd_->buf4_mat_irrep_rd_block(EAB, Gcl, cl, nlinks);
      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, EAB->matrix[Gcl][0], nrows,
          &(T2AA->matrix[Gab][ab][il]), nlinks, 1.0, W2[Gjk][0], ncols);
      global_dpd_->free_dpd_block(EAB->matrix[Gcl], nlinks, nrows);
    }
    global_dpd_->buf4_mat_irrep_close(T2AA, Gab);
  } /* Gl */

  /*=== Sort W2[Jk][I] -> W[IJ][k] ===*/
  global_dpd_->sort_3d(W2, W, nirreps, Gabc^GX3, T2AB->params->coltot,
      T2AB->params->colidx, T2AB->params->colorb, T2AB->params->rsym,
      T2AB->params->ssym, aocc_off, bocc_off, aoccpi, aocc_off,
      T2AA->params->colidx, cab, 1);

  /* clean out the alpha-beta-alpha intermediate for next set of terms */
  for (Gij=0; Gij < nirreps; ++Gij) {
    Gk = Gabc ^ Gij ^ GX3;
    if (T2AB->params->coltot[Gij] && aoccpi[Gk]) {
      memset(W2[Gij][0], 0, T2AB->params->coltot[Gij]*aoccpi[Gk]*sizeof(double));
    }
  }

  /*** compute disconnected triples ***/
  if (disc) {
    for (h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(DAB, h);
      global_dpd_->buf4_mat_irrep_rd(DAB, h);
      global_dpd_->buf4_mat_irrep_init(T2AB, h);
      global_dpd_->buf4_mat_irrep_rd(T2AB, h);
      global_dpd_->buf4_mat_irrep_init(DAA, h);
      global_dpd_->buf4_mat_irrep_rd(DAA, h);
      global_dpd_->buf4_mat_irrep_init(T2AA, h);
      global_dpd_->buf4_mat_irrep_rd(T2AA, h);
    }

    for (Gij=0; Gij < nirreps; ++Gij) {
      Gk = Gij ^ Gabc;
      Gji = Gij;

      for (ij=0; ij < T2AA->params->coltot[Gij]; ++ij) {
        I = T2AA->params->colorb[Gij][ij][0];
        Gi = T2AA->params->rsym[I];
        i = I - aocc_off[Gi];
        J = T2AA->params->colorb[Gij][ij][1];
        Gj = T2AA->params->ssym[J];
        j = J - aocc_off[Gj];

        Gjk = Gkj = Gj ^ Gk;
        Gik = Gki = Gi ^ Gk;

        for (k=0; k < boccpi[Gk]; ++k) {
          K = bocc_off[Gk] + k;

          ik = DAB->params->rowidx[I][K];
          jk = DAB->params->rowidx[J][K];
          ji = DAA->params->rowidx[J][I];

          /* +t_IA * D_JkBc + f_IA * t_JkBc */
          if (Gi == Ga && Gjk == Gbc) {
            t_ia = D_jkbc = f_ia = t_jkbc = 0.0;

            if (T1A->params->rowtot[Ga] && T1A->params->coltot[Ga]) {
              t_ia = T1A->matrix[Ga][i][a];
              f_ia = fIA->matrix[Ga][i][a];
            }

            if (DAB->params->rowtot[Gbc] && DAB->params->coltot[Gbc]) {
              D_jkbc = DAB->matrix[Gbc][jk][bc];
              t_jkbc = T2AB->matrix[Gbc][bc][jk];
            }
            V[Gij][ij][k] += t_ia * D_jkbc + f_ia * t_jkbc;
          }

          /* -t_IB * D_JkAc - f_IB * t_JkAc */
          if (Gi == Gb && Gjk == Gac) {
            t_ib = D_jkac = f_ib = t_jkac = 0.0;

            if (T1A->params->rowtot[Gb] && T1A->params->coltot[Gb]) {
              t_ib = T1A->matrix[Gb][i][b];
              f_ib = fIA->matrix[Gb][i][b];
            }

            if (DAB->params->rowtot[Gac] && DAB->params->coltot[Gac]) {
              D_jkac = DAB->matrix[Gac][jk][ac];
              t_jkac = T2AB->matrix[Gac][ac][jk];
            }
            V[Gij][ij][k] -= t_ib * D_jkac + f_ib * t_jkac;
          }

          /* -t_JA * D_IkBc - f_JA * t_IkBc */
          if (Gj == Ga && Gik == Gbc) {
            t_ja = D_ikbc = f_ja = t_ikbc = 0.0;

            if (T1A->params->rowtot[Ga] && T1A->params->coltot[Ga]) {
              t_ja = T1A->matrix[Ga][j][a];
              f_ja = fIA->matrix[Ga][j][a];
            }

            if (DAB->params->rowtot[Gbc] && DAB->params->coltot[Gbc]) {
              D_ikbc = DAB->matrix[Gbc][ik][bc];
              t_ikbc = T2AB->matrix[Gbc][bc][ik];
            }
            V[Gij][ij][k] -= t_ja * D_ikbc + f_ja * t_ikbc;
          }

          /* +t_JB * D_IkAc + f_JB * t_IkAc */
          if (Gj == Gb && Gik == Gac) {
            t_jb = D_ikac = f_jb = t_ikac = 0.0;

            if (T1A->params->rowtot[Gb] && T1A->params->coltot[Gb]) {
              t_jb = T1A->matrix[Gb][j][b];
              f_jb = fIA->matrix[Gb][j][b];
            }

            if (DAB->params->rowtot[Gac] && DAB->params->coltot[Gac]) {
              D_ikac = DAB->matrix[Gac][ik][ac];
              t_ikac = T2AB->matrix[Gac][ac][ik];
            }
            V[Gij][ij][k] += t_jb * D_ikac + f_jb * t_ikac;
          }

          /* -t_kc * D_JIAB - f_kc * t_JIAB */
          if (Gk == Gc && Gji == Gab) {
            t_kc = D_jiab = f_kc = t_jiab = 0.0;

            if (T1B->params->rowtot[Gc] && T1B->params->coltot[Gc]) {
              t_kc = T1B->matrix[Gc][k][c];
              f_kc = fia->matrix[Gc][k][c];
            }

            if (DAA->params->rowtot[Gab] && DAA->params->coltot[Gab]) {
              D_jiab = DAA->matrix[Gab][ji][ab];
              t_jiab = T2AA->matrix[Gab][ab][ji];
            }
            V[Gij][ij][k] -= t_kc * D_jiab + f_kc * t_jiab;
          }
        } /* k */
      } /* ij */
    } /* Gij */

    for (h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_close(T2AB, h);
      global_dpd_->buf4_mat_irrep_close(DAB, h);
      global_dpd_->buf4_mat_irrep_close(T2AA, h);
      global_dpd_->buf4_mat_irrep_close(DAA, h);
    }
  } /**** Disconnected T3 complete ****/

  for (Gij=0; Gij < nirreps; ++Gij) {
    Gk = Gabc ^ Gij ^ GX3; /* assumes totally symmetric! */
    for (ij=0; ij < T2AA->params->coltot[Gij]; ++ij) {
      I = T2AA->params->colorb[Gij][ij][0];
      J = T2AA->params->colorb[Gij][ij][1];
      Gi = T2AA->params->rsym[I];
      Gj = T2AA->params->ssym[J];
      i = I - aocc_off[Gi];
      j = J - aocc_off[Gj];

      for (k=0; k < boccpi[Gk]; ++k) {
        K = bocc_off[Gk] + k;

        denom = dabc;
        if (fIJ->params->rowtot[Gi])
          denom += fIJ->matrix[Gi][i][i];
        if (fIJ->params->rowtot[Gj])
          denom += fIJ->matrix[Gj][j][j];
        if (fij->params->rowtot[Gk])
          denom += fij->matrix[Gk][k][k];

        W[Gij][ij][k] /= (denom + omega);
        if (disc)
          V[Gij][ij][k] /= (denom + omega);
      } /* k */
    } /* ij */
  } /* Gij */

  for (Gij=0; Gij < nirreps; Gij++) {
    Gk = Gabc ^ Gij ^ GX3;
    global_dpd_->free_dpd_block(W2[Gij], T2BA->params->coltot[Gij], aoccpi[Gk]);
  }

  free(W2);

  global_dpd_->file2_mat_close(fIJ);
  global_dpd_->file2_mat_close(fij);
  global_dpd_->file2_mat_close(fAB);
  global_dpd_->file2_mat_close(fab);
  if (disc) {
    global_dpd_->file2_mat_close(fIA);
    global_dpd_->file2_mat_close(fia);
    global_dpd_->file2_mat_close(T1A);
    global_dpd_->file2_mat_close(T1B);
  }

  for (h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_close(FAA, h);
    global_dpd_->buf4_mat_irrep_close(FAB, h);
    global_dpd_->buf4_mat_irrep_close(FBA, h);
  }
}

}
}
